from pathlib import Path
import re

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st


# ============================================================
# CONFIGURAÇÃO
# ============================================================

st.set_page_config(
    page_title="Orthoherpesviridae Dataset Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


HERE = Path(__file__).resolve().parent

# Execução no servidor:
# ROOT/DASHBOARD_DATASET_INICIAL/dashboard_dataset_inicial.py
#
# Execução no Streamlit Community Cloud:
# ROOT/streamlit_app.py
if (HERE / "DASHBOARD_DATASET_INICIAL" / "data").exists():
    DATA_DIR = HERE / "DASHBOARD_DATASET_INICIAL" / "data"
else:
    DATA_DIR = HERE / "data"


RECORDS_FILE = DATA_DIR / "records_dashboard.tsv.gz"
CDS_FILE = DATA_DIR / "cds_dashboard.tsv.gz"
QUALITY_FILE = DATA_DIR / "metadata_quality.tsv"
DPOL_FILE = DATA_DIR / "dpol_flow.tsv"


# Paleta Okabe-Ito / color-blind friendly.
# Evitamos depender apenas de vermelho x verde.
ACCESSIBLE_COLORS = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#D55E00",  # vermillion
    "#56B4E9",  # sky blue
    "#000000",  # black
    "#999999",  # gray
]

ANNOTATION_CLASS_COLORS = {
    "Annotated product": "#0072B2",
    "Hypothetical with annotation clue": "#009E73",
    "Hypothetical with identifier only": "#E69F00",
    "Hypothetical without additional annotation": "#D55E00",
    "No product annotation": "#999999",
}

INTEREST_ORDER = [
    "CalHV-3",
    "Cebus albifrons lymphocryptovirus 1",
    "Cebus apella lymphocryptovirus 1",
    "Callithrix penicillata lymphocryptovirus 1",
    "Alouatta macconnelli cytomegalovirus",
    "Alouatta palliata cytomegalovirus",
    "Alouatta seniculus cytomegalovirus",
    "HSV-1",
]

MISSING_TERMS = {
    "",
    "unknown",
    "unknown_country",
    "unknown year",
    "unknown_year",
    "na",
    "n/a",
    "none",
    "null",
    "missing",
    "not provided",
    "not_provided",
    "not available",
    "not_available",
    ".",
    "-",
}

PLOTLY_CONFIG = {
    "displaylogo": False,
    "responsive": True,
}


# ============================================================
# FUNÇÕES
# ============================================================

def available_columns(df, columns):
    return [column for column in columns if column in df.columns]


def normalized(series):
    return (
        series
        .fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
    )


def missing_mask(series):
    return normalized(series).isin(MISSING_TERMS)


def download_tsv(df):
    return df.to_csv(sep="\t", index=False).encode("utf-8")


def download_on_demand(label, df, file_name, key):
    """Evita serializar tabelas grandes em cada rerun do Streamlit."""
    prepare = st.checkbox(
        f"Preparar arquivo completo para download — {label}",
        value=False,
        key=f"{key}_prepare",
        help=(
            "O arquivo só é convertido para TSV quando esta opção é marcada. "
            "Isso reduz bastante o uso de memória no Community Cloud."
        ),
    )
    if prepare:
        with st.spinner("Preparando arquivo para download..."):
            payload = download_tsv(df)
        st.download_button(
            label,
            data=payload,
            file_name=file_name,
            mime="text/tab-separated-values",
            key=key,
        )


def style_fig(fig, height=None, legend_title=None):
    fig.update_layout(
        template="plotly_white",
        font=dict(size=14),
        hoverlabel=dict(font_size=13),
        margin=dict(l=20, r=20, t=60, b=30),
        legend=dict(
            title=legend_title,
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
        ),
    )

    if height is not None:
        fig.update_layout(height=height)

    fig.update_xaxes(
        showgrid=True,
        gridcolor="rgba(120,120,120,0.16)",
        zeroline=False,
        automargin=True,
    )
    fig.update_yaxes(
        showgrid=False,
        zeroline=False,
        automargin=True,
    )
    return fig


def categorical_height(n_items, minimum=440, per_item=27, maximum=12000):
    return max(minimum, min(maximum, int(n_items) * per_item + 160))


def count_limit_control(label, key, default="50", include_200=True):
    options = ["20", "50", "100"]
    if include_200:
        options.append("200")
    options.append("Todas")

    index = options.index(default) if default in options else 1
    choice = st.selectbox(label, options, index=index, key=key)

    if choice == "Todas":
        return None
    return int(choice)


def limit_frame(df, n):
    if n is None:
        return df.copy()
    return df.head(n).copy()


def filter_dataframe(
    df,
    countries,
    species,
    interest,
    year_range,
    include_undated,
):
    # Não copie a tabela inteira sem necessidade.
    # Isso é importante no Community Cloud, principalmente para as 446.993 CDS.
    out = df

    if countries and "country_final" in out.columns:
        out = out[out["country_final"].isin(countries)]

    if species and "species_final" in out.columns:
        out = out[out["species_final"].isin(species)]

    if interest and "interest_group" in out.columns:
        out = out[out["interest_group"].isin(interest)]

    if year_range is not None and "collection_year_num" in out.columns:
        years = pd.to_numeric(out["collection_year_num"], errors="coerce")
        mask = years.between(year_range[0], year_range[1], inclusive="both")
        if include_undated:
            mask |= years.isna()
        out = out[mask]

    return out


def safe_value_counts(series, name, value_name):
    return (
        series
        .fillna("")
        .astype(str)
        .str.strip()
        .loc[lambda x: x.ne("")]
        .value_counts()
        .rename_axis(name)
        .reset_index(name=value_name)
    )


def taxon_protein_summary(df, dimension):
    base = df.copy()
    base[dimension] = base[dimension].fillna("").astype(str).str.strip()
    base = base[base[dimension].ne("")]

    totals = base.groupby(dimension).size().rename("CDS_total")

    hyp = (
        base[base["is_hypothetical"].eq("YES")]
        .groupby(dimension)
        .size()
        .rename("CDS_hypothetical")
    )

    clue = (
        base[
            base["protein_annotation_class"].eq(
                "Hypothetical with annotation clue"
            )
        ]
        .groupby(dimension)
        .size()
        .rename("Hypothetical_com_pista")
    )

    summary = pd.concat([totals, hyp, clue], axis=1).fillna(0).reset_index()
    summary["CDS_hypothetical"] = summary["CDS_hypothetical"].astype(int)
    summary["Hypothetical_com_pista"] = summary["Hypothetical_com_pista"].astype(int)
    summary["Taxa_hypothetical_pct"] = np.where(
        summary["CDS_total"] > 0,
        100 * summary["CDS_hypothetical"] / summary["CDS_total"],
        0,
    )
    return summary


# ============================================================
# VALIDAÇÃO DOS DADOS
# ============================================================

for required in [RECORDS_FILE, CDS_FILE, QUALITY_FILE, DPOL_FILE]:
    if not required.exists():
        st.error(
            "Base do dashboard não encontrada:\n\n"
            f"`{required}`\n\n"
            "Execute primeiro `01_build_dashboard_data.py`."
        )
        st.stop()


# ============================================================
# CARREGAMENTO
# ============================================================

@st.cache_resource(show_spinner="Carregando registros...")
def load_records():
    df = pd.read_csv(
        RECORDS_FILE,
        sep="\t",
        compression="gzip",
        low_memory=False,
        dtype_backend="pyarrow",
    )

    if "collection_year_num" in df.columns:
        df["collection_year_num"] = pd.to_numeric(
            df["collection_year_num"],
            errors="coerce",
        )

    return df


@st.cache_resource(show_spinner="Carregando CDS...")
def load_cds():
    df = pd.read_csv(
        CDS_FILE,
        sep="\t",
        compression="gzip",
        low_memory=False,
        dtype_backend="pyarrow",
    )

    if "collection_year_num" in df.columns:
        df["collection_year_num"] = pd.to_numeric(
            df["collection_year_num"],
            errors="coerce",
        )

    return df


@st.cache_resource
def load_simple_table(path):
    return pd.read_csv(path, sep="\t", dtype_backend="pyarrow")


records = load_records()
cds = load_cds()
quality = load_simple_table(QUALITY_FILE)
dpol_flow = load_simple_table(DPOL_FILE)


# ============================================================
# CABEÇALHO
# ============================================================

st.title("🧬 Orthoherpesviridae Dataset Explorer")
st.caption(
    "Exploração do universo inicial do NCBI/GenBank, taxonomia normalizada, "
    "metadados, hospedeiros, genes, CDS, proteínas hipotéticas e fluxo de "
    "curadoria DPOL. Os gráficos possuem controles para exibir Top N ou todas "
    "as categorias quando isso for informativo."
)


# ============================================================
# SIDEBAR
# ============================================================

st.sidebar.header("Filtros globais")
st.sidebar.caption(
    "Seleções vazias significam 'todos'. Os filtros são aplicados tanto aos "
    "registros quanto às CDS quando os metadados correspondentes estão disponíveis."
)

country_options = sorted(
    records.loc[
        records["country_final"].fillna("").ne("Unknown"),
        "country_final",
    ]
    .dropna()
    .unique()
    .tolist()
)

species_options = sorted(
    records.loc[
        records["species_final"].fillna("").ne("Unresolved_species"),
        "species_final",
    ]
    .dropna()
    .astype(str)
    .loc[lambda x: x.str.strip().ne("")]
    .unique()
    .tolist()
)

interest_options = [
    value
    for value in INTEREST_ORDER
    if value in set(records["interest_group"].dropna())
]

selected_interest = st.sidebar.multiselect(
    "Espécies / vírus de interesse",
    interest_options,
)

selected_countries = st.sidebar.multiselect(
    "País",
    country_options,
)

selected_species = st.sidebar.multiselect(
    "species_final",
    species_options,
)

year_values = records["collection_year_num"].dropna()

if len(year_values):
    min_year = int(year_values.min())
    max_year = int(year_values.max())
    selected_year = st.sidebar.slider(
        "Ano de coleta",
        min_value=min_year,
        max_value=max_year,
        value=(min_year, max_year),
    )
else:
    selected_year = None

include_undated = st.sidebar.checkbox(
    "Incluir registros sem ano",
    value=True,
)

filtered_records = filter_dataframe(
    records,
    selected_countries,
    selected_species,
    selected_interest,
    selected_year,
    include_undated,
)

filtered_cds = filter_dataframe(
    cds,
    selected_countries,
    selected_species,
    selected_interest,
    selected_year,
    include_undated,
)


# ============================================================
# MÉTRICAS GLOBAIS
# ============================================================

m1, m2, m3, m4 = st.columns(4)
m1.metric("Registros no dataset", f"{len(records):,}")
m2.metric("Registros após filtros", f"{len(filtered_records):,}")
m3.metric("CDS no dataset", f"{len(cds):,}")
m4.metric("CDS após filtros", f"{len(filtered_cds):,}")


# ============================================================
# ABAS
# ============================================================

# IMPORTANTE: abas dinâmicas. Apenas a aba aberta é calculada.
# Isso evita executar simultaneamente gráficos/tabelas das 446.993 CDS.
tabs = st.tabs([
    "🏠 Visão geral",
    "⭐ Espécies de interesse",
    "🔎 Explorador taxonômico",
    "🌍 Geografia",
    "📅 Tempo",
    "🇧🇷 Brasil",
    "🐒 Hospedeiros",
    "🧫 Genes / CDS",
    "❓ Proteínas hipotéticas",
    "🧬 Comparação de proteínas",
    "✅ Qualidade",
    "🧪 DPOL / Curadoria",
    "📋 Registros",
], on_change="rerun", key="main_tabs")


# ============================================================
# VISÃO GERAL
# ============================================================

if tabs[0].open:
    with tabs[0]:
        st.header("Visão geral")

        species_nonblank = (
            records["species_final"]
            .fillna("")
            .astype(str)
            .str.strip()
        )

        c1, c2, c3, c4 = st.columns(4)
        c1.metric(
            "Nomes distintos em species_final",
            f"{species_nonblank[species_nonblank.ne('')].nunique():,}",
        )
        c2.metric(
            "Países conhecidos",
            f"{records.loc[records['country_final'] != 'Unknown', 'country_final'].nunique():,}",
        )
        c3.metric(
            "Registros brasileiros",
            f"{(records['country_final'] == 'Brazil').sum():,}",
        )
        c4.metric(
            "CDS hipotéticas",
            f"{(cds['is_hypothetical'] == 'YES').sum():,}",
        )

        st.subheader("Espécies mais representadas")
        species_all = safe_value_counts(
            filtered_records["species_final"],
            "Espécie",
            "Registros",
        )
        n_species = count_limit_control(
            "Quantidade de espécies no gráfico",
            "overview_species_n",
            default="50",
        )
        species_plot = limit_frame(species_all, n_species).sort_values("Registros")

        if len(species_plot):
            fig = px.bar(
                species_plot,
                x="Registros",
                y="Espécie",
                orientation="h",
                color_discrete_sequence=[ACCESSIBLE_COLORS[0]],
                title="Representação por species_final",
            )
            style_fig(fig, categorical_height(len(species_plot)))
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        with st.expander("Ver tabela completa de espécies"):
            st.dataframe(species_all, hide_index=True, width="stretch")
            st.download_button(
                "Baixar contagem completa de espécies",
                download_tsv(species_all),
                file_name="species_counts.tsv",
                mime="text/tab-separated-values",
                key="dl_species_counts",
            )

        st.subheader("Espécies / vírus destacados")
        highlighted = filtered_records[
            filtered_records["interest_group"].ne("Other")
        ]
        interest_counts = (
            highlighted["interest_group"]
            .value_counts()
            .rename_axis("Grupo")
            .reset_index(name="Registros")
            .sort_values("Registros")
        )

        if len(interest_counts):
            fig = px.bar(
                interest_counts,
                x="Registros",
                y="Grupo",
                orientation="h",
                color="Grupo",
                color_discrete_sequence=ACCESSIBLE_COLORS,
                text="Registros",
                title="Grupos de interesse",
            )
            style_fig(fig, categorical_height(len(interest_counts), minimum=440))
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)


    # ============================================================
    # ESPÉCIES DE INTERESSE
    # ============================================================

if tabs[1].open:
    with tabs[1]:
        st.header("Espécies e vírus de interesse")
        st.caption(
            "Esta categoria é apenas uma camada de visualização. A taxonomia "
            "normalizada permanece preservada em species_final."
        )

        focus = records[records["interest_group"].ne("Other")].copy()
        counts = focus["interest_group"].value_counts()

        metric_cols = st.columns(4)
        for index, label in enumerate(INTEREST_ORDER):
            metric_cols[index % 4].metric(label, f"{int(counts.get(label, 0)):,}")

        selected_compare = st.multiselect(
            "Comparar",
            INTEREST_ORDER,
            default=[label for label in INTEREST_ORDER if counts.get(label, 0) > 0],
            key="interest_compare",
        )

        compare = focus[focus["interest_group"].isin(selected_compare)]

        if len(compare):
            st.subheader("Distribuição geográfica")
            interest_country = (
                compare[compare["country_final"].ne("Unknown")]
                .groupby(["country_final", "interest_group"])
                .size()
                .reset_index(name="Registros")
            )

            if len(interest_country):
                country_totals = (
                    interest_country.groupby("country_final")["Registros"]
                    .sum()
                    .sort_values(ascending=False)
                )
                n_countries = count_limit_control(
                    "Países no gráfico",
                    "interest_country_n",
                    default="50",
                )
                keep = country_totals.index if n_countries is None else country_totals.head(n_countries).index
                tmp = interest_country[interest_country["country_final"].isin(keep)]

                fig = px.bar(
                    tmp,
                    x="country_final",
                    y="Registros",
                    color="interest_group",
                    barmode="group",
                    color_discrete_sequence=ACCESSIBLE_COLORS,
                    labels={
                        "country_final": "País",
                        "interest_group": "Espécie / vírus",
                    },
                    title="País × grupo de interesse",
                )
                style_fig(fig, height=560, legend_title="Espécie / vírus")
                fig.update_xaxes(tickangle=-45)
                st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

            st.subheader("Distribuição temporal")
            interest_year = (
                compare[compare["collection_year_num"].notna()]
                .groupby(["collection_year_num", "interest_group"])
                .size()
                .reset_index(name="Registros")
            )

            if len(interest_year):
                fig = px.line(
                    interest_year,
                    x="collection_year_num",
                    y="Registros",
                    color="interest_group",
                    markers=True,
                    color_discrete_sequence=ACCESSIBLE_COLORS,
                    labels={
                        "collection_year_num": "Ano",
                        "interest_group": "Espécie / vírus",
                    },
                    title="Registros por ano",
                )
                style_fig(fig, height=520, legend_title="Espécie / vírus")
                st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

            st.subheader("Hospedeiros dos grupos selecionados")
            focus_hosts = compare[~missing_mask(compare["host"])].copy()
            if len(focus_hosts):
                host_group = (
                    focus_hosts.groupby(["host", "interest_group"])
                    .size()
                    .reset_index(name="Registros")
                )
                host_totals = (
                    host_group.groupby("host")["Registros"]
                    .sum()
                    .sort_values(ascending=False)
                )
                n_hosts_interest = count_limit_control(
                    "Hospedeiros no gráfico",
                    "interest_hosts_n",
                    default="50",
                )
                keep_hosts = host_totals.index if n_hosts_interest is None else host_totals.head(n_hosts_interest).index
                host_group = host_group[host_group["host"].isin(keep_hosts)]

                fig = px.bar(
                    host_group,
                    x="Registros",
                    y="host",
                    color="interest_group",
                    orientation="h",
                    barmode="stack",
                    color_discrete_sequence=ACCESSIBLE_COLORS,
                    labels={"host": "Hospedeiro", "interest_group": "Espécie / vírus"},
                    title="Hospedeiro × grupo de interesse",
                )
                style_fig(
                    fig,
                    categorical_height(host_group["host"].nunique()),
                    legend_title="Espécie / vírus",
                )
                st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

            st.subheader("Nomes originais no GenBank")
            original_names = (
                compare.groupby(["interest_group", "organism"])
                .size()
                .reset_index(name="Registros")
                .sort_values(["interest_group", "Registros"], ascending=[True, False])
            )
            st.dataframe(original_names, hide_index=True, width="stretch")


    # ============================================================
    # EXPLORADOR TAXONÔMICO
    # ============================================================

if tabs[2].open:
    with tabs[2]:
        st.header("Explorador taxonômico")

        explorer_mode = st.radio(
            "Usar",
            ["species_final", "organism"],
            horizontal=True,
        )

        explorer_options = sorted(
            records[explorer_mode]
            .dropna()
            .astype(str)
            .loc[lambda x: x.str.strip().ne("")]
            .unique()
            .tolist()
        )

        explorer_selected = st.multiselect(
            "Selecione um ou mais táxons",
            explorer_options,
            key="explorer_taxa",
        )

        if explorer_selected:
            explorer = records[records[explorer_mode].isin(explorer_selected)]

            e1, e2, e3, e4 = st.columns(4)
            e1.metric("Registros", f"{len(explorer):,}")
            e2.metric(
                "Países",
                f"{explorer.loc[explorer['country_final'] != 'Unknown', 'country_final'].nunique():,}",
            )
            e3.metric(
                "Hospedeiros",
                f"{explorer.loc[~missing_mask(explorer['host']), 'host'].nunique():,}",
            )
            e4.metric(
                "Anos",
                f"{explorer['collection_year_num'].dropna().nunique():,}",
            )

            explorer_country = safe_value_counts(
                explorer.loc[explorer["country_final"].ne("Unknown"), "country_final"],
                "País",
                "Registros",
            )
            n_explorer_country = count_limit_control(
                "Países no gráfico",
                "explorer_country_n",
                default="50",
            )
            explorer_country_plot = limit_frame(explorer_country, n_explorer_country).sort_values("Registros")

            if len(explorer_country_plot):
                fig = px.bar(
                    explorer_country_plot,
                    x="Registros",
                    y="País",
                    orientation="h",
                    color_discrete_sequence=[ACCESSIBLE_COLORS[0]],
                    title="Distribuição geográfica dos táxons selecionados",
                )
                style_fig(fig, categorical_height(len(explorer_country_plot)))
                st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

            explorer_year = (
                explorer[explorer["collection_year_num"].notna()]
                .groupby(["collection_year_num", explorer_mode])
                .size()
                .reset_index(name="Registros")
            )
            if len(explorer_year):
                fig = px.line(
                    explorer_year,
                    x="collection_year_num",
                    y="Registros",
                    color=explorer_mode,
                    markers=True,
                    color_discrete_sequence=ACCESSIBLE_COLORS,
                    labels={"collection_year_num": "Ano"},
                    title="Distribuição temporal",
                )
                style_fig(fig, height=520, legend_title=explorer_mode)
                st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

            explorer_hosts = safe_value_counts(
                explorer.loc[~missing_mask(explorer["host"]), "host"],
                "Hospedeiro",
                "Registros",
            )
            n_explorer_hosts = count_limit_control(
                "Hospedeiros no gráfico",
                "explorer_hosts_n",
                default="50",
            )
            explorer_hosts_plot = limit_frame(explorer_hosts, n_explorer_hosts).sort_values("Registros")

            if len(explorer_hosts_plot):
                fig = px.bar(
                    explorer_hosts_plot,
                    x="Registros",
                    y="Hospedeiro",
                    orientation="h",
                    color_discrete_sequence=[ACCESSIBLE_COLORS[2]],
                    title="Hospedeiros dos táxons selecionados",
                )
                style_fig(fig, categorical_height(len(explorer_hosts_plot)))
                st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

            display = available_columns(
                explorer,
                [
                    "accession_version",
                    "taxid",
                    "organism",
                    "species_final",
                    "taxonomy_final_status",
                    "country_final",
                    "collection_date",
                    "collection_year",
                    "host",
                    "isolate",
                    "strain",
                    "definition",
                ],
            )
            st.dataframe(explorer[display], hide_index=True, width="stretch")
            st.download_button(
                "Baixar registros selecionados",
                download_tsv(explorer[display]),
                file_name="taxon_explorer.tsv",
                mime="text/tab-separated-values",
            )
        else:
            st.info("Selecione um ou mais táxons para iniciar a exploração.")


    # ============================================================
    # GEOGRAFIA
    # ============================================================

if tabs[3].open:
    with tabs[3]:
        st.header("Distribuição geográfica")
        geo = filtered_records[filtered_records["country_final"].ne("Unknown")]

        countries_all = (
            geo["country_final"]
            .value_counts()
            .rename_axis("País")
            .reset_index(name="Registros")
        )
        n_countries_geo = count_limit_control(
            "Países no gráfico de registros",
            "geo_country_n",
            default="50",
        )
        countries = limit_frame(countries_all, n_countries_geo).sort_values("Registros")

        if len(countries):
            fig = px.bar(
                countries,
                x="Registros",
                y="País",
                orientation="h",
                color_discrete_sequence=[ACCESSIBLE_COLORS[0]],
                title="Registros por país",
            )
            style_fig(fig, categorical_height(len(countries)))
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        if len(countries_all):
            fig = px.choropleth(
                countries_all,
                locations="País",
                locationmode="country names",
                color="Registros",
                hover_name="País",
                color_continuous_scale="Cividis",
                title="Mapa mundial de registros",
            )
            style_fig(fig, height=600)
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        diversity_all = (
            geo[geo["species_final"].ne("Unresolved_species")]
            .groupby("country_final")["species_final"]
            .nunique()
            .sort_values(ascending=False)
            .rename_axis("País")
            .reset_index(name="Espécies distintas")
        )

        st.subheader("Diversidade taxonômica")
        n_diversity = count_limit_control(
            "Países no gráfico de diversidade",
            "geo_diversity_n",
            default="50",
        )
        diversity = limit_frame(diversity_all, n_diversity).sort_values("Espécies distintas")

        if len(diversity):
            fig = px.bar(
                diversity,
                x="Espécies distintas",
                y="País",
                orientation="h",
                color_discrete_sequence=[ACCESSIBLE_COLORS[4]],
                title="Número de species_final distintos por país",
            )
            style_fig(fig, categorical_height(len(diversity)))
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)


    # ============================================================
    # TEMPO
    # ============================================================

if tabs[4].open:
    with tabs[4]:
        st.header("Distribuição temporal")
        dated = filtered_records[filtered_records["collection_year_num"].notna()]

        annual = (
            dated.groupby("collection_year_num")
            .agg(
                Registros=("accession_version", "count"),
                Espécies=("species_final", "nunique"),
                Países=("country_final", "nunique"),
            )
            .reset_index()
            .sort_values("collection_year_num")
        )

        if len(annual):
            fig = px.line(
                annual,
                x="collection_year_num",
                y="Registros",
                markers=True,
                color_discrete_sequence=[ACCESSIBLE_COLORS[0]],
                labels={"collection_year_num": "Ano"},
                title="Registros por ano de coleta",
            )
            style_fig(fig, height=500)
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

            fig = px.line(
                annual,
                x="collection_year_num",
                y="Espécies",
                markers=True,
                color_discrete_sequence=[ACCESSIBLE_COLORS[4]],
                labels={"collection_year_num": "Ano"},
                title="Diversidade de species_final por ano",
            )
            style_fig(fig, height=500)
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

            with st.expander("Ver série temporal completa"):
                st.dataframe(annual, hide_index=True, width="stretch")


    # ============================================================
    # BRASIL
    # ============================================================

if tabs[5].open:
    with tabs[5]:
        st.header("Brasil")
        brazil = records[records["country_final"].eq("Brazil")]

        b1, b2, b3, b4 = st.columns(4)
        b1.metric("Registros", f"{len(brazil):,}")
        b2.metric("species_final distintos", f"{brazil['species_final'].nunique():,}")
        b3.metric(
            "Hospedeiros distintos",
            f"{brazil.loc[~missing_mask(brazil['host']), 'host'].nunique():,}",
        )
        b4.metric(
            "% do dataset",
            f"{100 * len(brazil) / len(records):.2f}%" if len(records) else "0%",
        )

        br_species_all = safe_value_counts(
            brazil["species_final"],
            "Espécie",
            "Registros",
        )
        n_br_species = count_limit_control(
            "Espécies no gráfico",
            "br_species_n",
            default="50",
        )
        br_species = limit_frame(br_species_all, n_br_species).sort_values("Registros")

        if len(br_species):
            fig = px.bar(
                br_species,
                x="Registros",
                y="Espécie",
                orientation="h",
                color_discrete_sequence=[ACCESSIBLE_COLORS[0]],
                title="Espécies mais representadas no Brasil",
            )
            style_fig(fig, categorical_height(len(br_species)))
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        br_interest = (
            brazil[brazil["interest_group"].ne("Other")]["interest_group"]
            .value_counts()
            .rename_axis("Grupo")
            .reset_index(name="Registros")
            .sort_values("Registros")
        )

        if len(br_interest):
            st.subheader("Espécies de interesse no Brasil")
            fig = px.bar(
                br_interest,
                x="Registros",
                y="Grupo",
                orientation="h",
                color="Grupo",
                color_discrete_sequence=ACCESSIBLE_COLORS,
                text="Registros",
            )
            style_fig(fig, categorical_height(len(br_interest)))
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        br_hosts_all = safe_value_counts(
            brazil.loc[~missing_mask(brazil["host"]), "host"],
            "Hospedeiro",
            "Registros",
        )
        n_br_hosts = count_limit_control(
            "Hospedeiros no gráfico",
            "br_hosts_n",
            default="50",
        )
        br_hosts = limit_frame(br_hosts_all, n_br_hosts).sort_values("Registros")

        if len(br_hosts):
            st.subheader("Hospedeiros no Brasil")
            fig = px.bar(
                br_hosts,
                x="Registros",
                y="Hospedeiro",
                orientation="h",
                color_discrete_sequence=[ACCESSIBLE_COLORS[2]],
            )
            style_fig(fig, categorical_height(len(br_hosts)))
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)


    # ============================================================
    # HOSPEDEIROS — VERSÃO EXPANDIDA
    # ============================================================

if tabs[6].open:
    with tabs[6]:
        st.header("Hospedeiros")
        st.caption(
            "Os nomes de hospedeiro são os metadados disponíveis nos registros e "
            "podem conter variações de nomenclatura. A aba mostra tanto a frequência "
            "bruta quanto a associação hospedeiro × táxon."
        )

        host_records = filtered_records[
            ~missing_mask(filtered_records["host"])
        ].copy()

        h1, h2, h3, h4 = st.columns(4)
        h1.metric("Registros com hospedeiro", f"{len(host_records):,}")
        h2.metric("Hospedeiros distintos", f"{host_records['host'].nunique():,}")
        h3.metric(
            "species_final ligados a hospedeiro",
            f"{host_records['species_final'].nunique():,}",
        )
        h4.metric(
            "Completude",
            f"{100 * len(host_records) / len(filtered_records):.1f}%" if len(filtered_records) else "0%",
        )

        host_counts_all = safe_value_counts(
            host_records["host"],
            "Hospedeiro",
            "Registros",
        )
        n_hosts = count_limit_control(
            "Quantidade de hospedeiros no ranking",
            "hosts_ranking_n",
            default="50",
        )
        host_counts = limit_frame(host_counts_all, n_hosts)

        st.subheader("Ranking de hospedeiros")
        if len(host_counts):
            fig = px.bar(
                host_counts.sort_values("Registros"),
                x="Registros",
                y="Hospedeiro",
                orientation="h",
                color="Registros",
                color_continuous_scale="Cividis",
                title="Número de registros por hospedeiro",
            )
            style_fig(fig, categorical_height(len(host_counts)))
            fig.update_layout(coloraxis_colorbar=dict(title="Registros"))
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        st.subheader("Visão proporcional")
        n_treemap = count_limit_control(
            "Hospedeiros no treemap",
            "hosts_treemap_n",
            default="50",
        )
        treemap_data = limit_frame(host_counts_all, n_treemap)
        if len(treemap_data):
            fig = px.treemap(
                treemap_data,
                path=["Hospedeiro"],
                values="Registros",
                color="Registros",
                color_continuous_scale="Cividis",
                title="Proporção relativa dos hospedeiros",
            )
            style_fig(fig, height=650)
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        st.subheader("Hospedeiro × táxon")
        host_matrix_dimension = st.radio(
            "Agrupar vírus por",
            ["species_final", "interest_group"],
            horizontal=True,
            key="host_matrix_dimension",
        )

        host_selection_mode = st.radio(
            "Hospedeiros da matriz",
            ["Mais frequentes", "Escolher manualmente"],
            horizontal=True,
            key="host_selection_mode",
        )

        if host_selection_mode == "Mais frequentes":
            n_matrix_hosts = st.slider(
                "Número de hospedeiros",
                min_value=2,
                max_value=min(30, max(2, host_counts_all.shape[0])),
                value=min(10, max(2, host_counts_all.shape[0])),
                key="host_matrix_n",
            )
            selected_hosts_matrix = host_counts_all.head(n_matrix_hosts)["Hospedeiro"].tolist()
        else:
            selected_hosts_matrix = st.multiselect(
                "Selecione os hospedeiros",
                host_counts_all["Hospedeiro"].tolist(),
                max_selections=30,
                key="host_matrix_manual",
            )

        if selected_hosts_matrix:
            hm = host_records[host_records["host"].isin(selected_hosts_matrix)].copy()
            if host_matrix_dimension == "interest_group":
                hm = hm[hm["interest_group"].ne("Other")]

            if len(hm):
                taxon_totals = (
                    hm[host_matrix_dimension]
                    .fillna("")
                    .astype(str)
                    .str.strip()
                    .loc[lambda x: x.ne("")]
                    .value_counts()
                )
                max_taxa = min(25, len(taxon_totals))
                if max_taxa > 0:
                    n_taxa_matrix = st.slider(
                        "Táxons na matriz",
                        min_value=1,
                        max_value=max_taxa,
                        value=min(12, max_taxa),
                        key="host_matrix_taxa_n",
                    )
                    taxa_keep = taxon_totals.head(n_taxa_matrix).index
                    hm = hm[hm[host_matrix_dimension].isin(taxa_keep)]

                    matrix = pd.crosstab(
                        hm["host"],
                        hm[host_matrix_dimension],
                    )

                    fig = px.imshow(
                        matrix,
                        aspect="auto",
                        color_continuous_scale="Cividis",
                        labels={
                            "x": host_matrix_dimension,
                            "y": "Hospedeiro",
                            "color": "Registros",
                        },
                        title="Matriz de associação hospedeiro × táxon",
                    )
                    style_fig(fig, height=max(500, 42 * len(matrix) + 220))
                    st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

                    st.dataframe(
                        matrix.reset_index(),
                        hide_index=True,
                        width="stretch",
                    )

        with st.expander("Tabela completa de hospedeiros"):
            st.dataframe(host_counts_all, hide_index=True, width="stretch")
            st.download_button(
                "Baixar contagem completa de hospedeiros",
                download_tsv(host_counts_all),
                file_name="host_counts.tsv",
                mime="text/tab-separated-values",
                key="dl_host_counts",
            )


    # ============================================================
    # GENES / CDS
    # ============================================================

if tabs[7].open:
    with tabs[7]:
        st.header("Genes e CDS")

        annotation_counts = (
            filtered_cds["protein_annotation_class"]
            .value_counts()
            .rename_axis("Classe")
            .reset_index(name="CDS")
        )
        if len(annotation_counts):
            fig = px.bar(
                annotation_counts,
                x="Classe",
                y="CDS",
                color="Classe",
                color_discrete_map=ANNOTATION_CLASS_COLORS,
                text="CDS",
                title="Classes de anotação proteica",
            )
            style_fig(fig, height=500, legend_title="Classe")
            fig.update_xaxes(tickangle=-20)
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        if "gene" in filtered_cds.columns:
            gene_counts_all = safe_value_counts(
                filtered_cds.loc[~missing_mask(filtered_cds["gene"]), "gene"],
                "Gene",
                "CDS",
            )
            n_genes = count_limit_control(
                "Genes no gráfico",
                "genes_n",
                default="50",
            )
            genes = limit_frame(gene_counts_all, n_genes).sort_values("CDS")

            if len(genes):
                st.subheader("Genes mais representados")
                fig = px.bar(
                    genes,
                    x="CDS",
                    y="Gene",
                    orientation="h",
                    color_discrete_sequence=[ACCESSIBLE_COLORS[0]],
                )
                style_fig(fig, categorical_height(len(genes)))
                st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        if "product" in filtered_cds.columns:
            product_counts_all = safe_value_counts(
                filtered_cds.loc[~missing_mask(filtered_cds["product"]), "product"],
                "Produto",
                "CDS",
            )
            n_products = count_limit_control(
                "Produtos no gráfico",
                "products_n",
                default="50",
            )
            products = limit_frame(product_counts_all, n_products).sort_values("CDS")

            if len(products):
                st.subheader("Produtos mais representados")
                fig = px.bar(
                    products,
                    x="CDS",
                    y="Produto",
                    orientation="h",
                    color_discrete_sequence=[ACCESSIBLE_COLORS[2]],
                )
                style_fig(fig, categorical_height(len(products)))
                st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        st.subheader("Busca nas anotações")
        query = st.text_input(
            "Pesquisar",
            placeholder="polymerase, capsid, tegument, ORF, kinase...",
            key="cds_search",
        )

        cds_search_fields = available_columns(
            filtered_cds,
            [
                "gene",
                "gene_synonym",
                "product",
                "standard_name",
                "function",
                "note",
                "protein_id",
                "locus_tag",
                "inference",
                "organism",
                "species_final",
                "annotation_clue",
            ],
        )

        search_hits = filtered_cds.copy()
        if query:
            pattern = re.escape(query)
            mask = pd.Series(False, index=search_hits.index)
            for column in cds_search_fields:
                mask |= (
                    search_hits[column]
                    .fillna("")
                    .astype(str)
                    .str.contains(pattern, case=False, regex=True, na=False)
                )
            search_hits = search_hits[mask]

        st.write(f"CDS encontradas: **{len(search_hits):,}**")

        display = available_columns(
            search_hits,
            [
                "accession_version",
                "organism",
                "species_final",
                "country_final",
                "host",
                "gene",
                "gene_synonym",
                "product",
                "standard_name",
                "function",
                "note",
                "protein_id",
                "locus_tag",
                "inference",
                "annotation_clue",
                "protein_annotation_class",
            ],
        )

        preview_n = st.selectbox(
            "Linhas na prévia da tabela",
            [1000, 5000, 10000, 20000],
            index=0,
            key="cds_preview_n",
        )
        st.dataframe(
            search_hits[display].head(preview_n),
            hide_index=True,
            width="stretch",
        )
        if len(search_hits) > preview_n:
            st.caption(
                f"A prévia mostra {preview_n:,} linhas; o download contém todas as "
                f"{len(search_hits):,} CDS filtradas."
            )

        download_on_demand(
            "Baixar resultados completos",
            search_hits[display],
            "CDS_search.tsv",
            "dl_cds_search_full",
        )


    # ============================================================
    # PROTEÍNAS HIPOTÉTICAS — VERSÃO EXPANDIDA
    # ============================================================

if tabs[8].open:
    with tabs[8]:
        st.header("Proteínas hipotéticas")
        st.info(
            "'Hypothetical protein' indica uma CDS prevista sem atribuição funcional "
            "específica suficientemente estabelecida na anotação. O painel separa "
            "pistas de anotação funcional de identificadores como protein_id/locus_tag."
        )

        annotation_classes = (
            filtered_cds["protein_annotation_class"]
            .value_counts()
            .rename_axis("Classe")
            .reset_index(name="CDS")
        )

        if len(annotation_classes):
            fig = px.bar(
                annotation_classes,
                x="Classe",
                y="CDS",
                color="Classe",
                color_discrete_map=ANNOTATION_CLASS_COLORS,
                text="CDS",
                title="Estado da anotação proteica",
            )
            style_fig(fig, height=520, legend_title="Classe")
            fig.update_xaxes(tickangle=-20)
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        hyp = filtered_cds[filtered_cds["is_hypothetical"].eq("YES")].copy()

        h1, h2, h3, h4 = st.columns(4)
        h1.metric("Hypothetical", f"{len(hyp):,}")
        h2.metric(
            "Com pista de anotação",
            f"{(hyp['protein_annotation_class'] == 'Hypothetical with annotation clue').sum():,}",
        )
        h3.metric(
            "Somente identificador",
            f"{(hyp['protein_annotation_class'] == 'Hypothetical with identifier only').sum():,}",
        )
        h4.metric(
            "% das CDS filtradas",
            f"{100 * len(hyp) / len(filtered_cds):.2f}%" if len(filtered_cds) else "0%",
        )

        st.subheader("Onde estão as proteínas hipotéticas?")
        hyp_dimension = st.radio(
            "Agrupar por",
            ["species_final", "organism", "interest_group"],
            horizontal=True,
            key="hyp_dimension",
        )

        hyp_summary = taxon_protein_summary(filtered_cds, hyp_dimension)
        if hyp_dimension == "interest_group":
            hyp_summary = hyp_summary[hyp_summary[hyp_dimension].ne("Other")]

        min_cds = st.number_input(
            "Mínimo de CDS no grupo para calcular a taxa",
            min_value=1,
            value=1,
            step=1,
            key="hyp_min_cds",
        )
        hyp_summary = hyp_summary[hyp_summary["CDS_total"] >= min_cds]

        hyp_rank_mode = st.radio(
            "Ordenar por",
            ["Número de hypothetical", "% hypothetical"],
            horizontal=True,
            key="hyp_rank_mode",
        )
        sort_col = "CDS_hypothetical" if hyp_rank_mode == "Número de hypothetical" else "Taxa_hypothetical_pct"
        hyp_summary = hyp_summary.sort_values(sort_col, ascending=False)

        n_hyp_taxa = count_limit_control(
            "Táxons no gráfico",
            "hyp_taxa_n",
            default="50",
        )
        hyp_taxa_plot = limit_frame(hyp_summary, n_hyp_taxa).sort_values(sort_col)

        if len(hyp_taxa_plot):
            x_col = "CDS_hypothetical" if hyp_rank_mode == "Número de hypothetical" else "Taxa_hypothetical_pct"
            x_label = "CDS hypothetical" if hyp_rank_mode == "Número de hypothetical" else "% hypothetical"

            fig = px.bar(
                hyp_taxa_plot,
                x=x_col,
                y=hyp_dimension,
                orientation="h",
                color="Hypothetical_com_pista",
                color_continuous_scale="Cividis",
                hover_data={
                    "CDS_total": True,
                    "CDS_hypothetical": True,
                    "Hypothetical_com_pista": True,
                    "Taxa_hypothetical_pct": ":.2f",
                },
                labels={x_col: x_label, hyp_dimension: "Táxon/grupo"},
                title="Distribuição das CDS hipotéticas por táxon/grupo",
            )
            style_fig(fig, categorical_height(len(hyp_taxa_plot)))
            fig.update_layout(coloraxis_colorbar=dict(title="Com pista"))
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        st.subheader("Tipos de pista encontrados")
        clue_hyp = hyp[
            hyp["protein_annotation_class"].eq("Hypothetical with annotation clue")
        ].copy()

        clue_source_counts = safe_value_counts(
            clue_hyp["annotation_clue_source"],
            "Campo de origem",
            "CDS",
        )
        if len(clue_source_counts):
            fig = px.bar(
                clue_source_counts.sort_values("CDS"),
                x="CDS",
                y="Campo de origem",
                orientation="h",
                color="Campo de origem",
                color_discrete_sequence=ACCESSIBLE_COLORS,
                title="Campo que forneceu a primeira pista adicional",
            )
            style_fig(fig, categorical_height(len(clue_source_counts), minimum=420))
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        clue_counts_all = safe_value_counts(
            clue_hyp["annotation_clue"],
            "Pista de anotação",
            "CDS",
        )
        n_clues = count_limit_control(
            "Pistas no gráfico",
            "hyp_clues_n",
            default="50",
        )
        clue_counts = limit_frame(clue_counts_all, n_clues).sort_values("CDS")

        if len(clue_counts):
            fig = px.bar(
                clue_counts,
                x="CDS",
                y="Pista de anotação",
                orientation="h",
                color_discrete_sequence=[ACCESSIBLE_COLORS[2]],
                title="Pistas de anotação mais recorrentes entre hypothetical proteins",
            )
            style_fig(fig, categorical_height(len(clue_counts)))
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        st.subheader("Explorar CDS hipotéticas")
        class_filter = st.multiselect(
            "Categoria",
            sorted(hyp["protein_annotation_class"].dropna().unique().tolist()),
            key="hyp_class_filter",
        )
        if class_filter:
            hyp = hyp[hyp["protein_annotation_class"].isin(class_filter)]

        hyp_query = st.text_input(
            "Pesquisar nas proteínas hipotéticas",
            placeholder="ORF, tegument, capsid, locus tag, protein ID, espécie...",
            key="hyp_query",
        )

        if hyp_query:
            pattern = re.escape(hyp_query)
            mask = pd.Series(False, index=hyp.index)
            for column in available_columns(
                hyp,
                [
                    "organism",
                    "species_final",
                    "gene",
                    "gene_synonym",
                    "annotation_clue",
                    "note",
                    "function",
                    "protein_id",
                    "locus_tag",
                    "inference",
                ],
            ):
                mask |= (
                    hyp[column]
                    .fillna("")
                    .astype(str)
                    .str.contains(pattern, case=False, regex=True, na=False)
                )
            hyp = hyp[mask]

        hyp_display = available_columns(
            hyp,
            [
                "accession_version",
                "organism",
                "species_final",
                "interest_group",
                "country_final",
                "host",
                "gene",
                "gene_synonym",
                "product",
                "standard_name",
                "function",
                "note",
                "protein_id",
                "locus_tag",
                "inference",
                "annotation_clue",
                "annotation_clue_source",
                "identifier_clue",
                "identifier_clue_source",
                "protein_annotation_class",
            ],
        )

        st.write(f"CDS hipotéticas exibidas: **{len(hyp):,}**")
        hyp_preview_n = st.selectbox(
            "Linhas na prévia",
            [1000, 5000, 10000, 20000],
            index=0,
            key="hyp_preview_n",
        )
        st.dataframe(
            hyp[hyp_display].head(hyp_preview_n),
            hide_index=True,
            width="stretch",
        )
        if len(hyp) > hyp_preview_n:
            st.caption(
                f"A prévia mostra {hyp_preview_n:,} linhas; o download contém todas "
                f"as {len(hyp):,} CDS hipotéticas filtradas."
            )

        download_on_demand(
            "Baixar proteínas hipotéticas filtradas",
            hyp[hyp_display],
            "hypothetical_proteins.tsv",
            "dl_hypothetical_full",
        )


    # ============================================================
    # COMPARAÇÃO DE PROTEÍNAS — VERSÃO EXPANDIDA
    # ============================================================

if tabs[9].open:
    with tabs[9]:
        st.header("Comparação de proteínas")
        st.caption(
            "Compare a composição de produtos, genes ou pistas de anotação entre "
            "táxons. A matriz pode mostrar contagem, presença/ausência ou proporção "
            "das CDS dentro de cada táxon."
        )

        compare_dimension = st.radio(
            "Comparar táxons por",
            ["species_final", "organism", "interest_group"],
            horizontal=True,
            key="protein_compare_dimension",
        )

        comparison_base = filtered_cds.copy()
        comparison_base[compare_dimension] = (
            comparison_base[compare_dimension]
            .fillna("")
            .astype(str)
            .str.strip()
        )
        comparison_base = comparison_base[comparison_base[compare_dimension].ne("")]

        if compare_dimension == "interest_group":
            comparison_base = comparison_base[
                comparison_base["interest_group"].ne("Other")
            ]

        compare_options = sorted(
            comparison_base[compare_dimension].dropna().unique().tolist()
        )

        default_taxa = []
        if compare_dimension == "interest_group":
            default_taxa = [x for x in INTEREST_ORDER if x in compare_options]

        compare_taxa = st.multiselect(
            "Selecione os táxons",
            compare_options,
            default=default_taxa,
            max_selections=15,
            key="protein_compare_taxa",
        )

        feature_dimension = st.radio(
            "Comparar características por",
            ["product", "gene", "annotation_clue"],
            horizontal=True,
            key="feature_dimension",
        )

        metric_mode = st.radio(
            "Métrica da matriz",
            ["Contagem de CDS", "% das CDS do táxon", "Presença / ausência"],
            horizontal=True,
            key="matrix_metric",
        )

        if compare_taxa:
            selected_base = comparison_base[
                comparison_base[compare_dimension].isin(compare_taxa)
            ].copy()

            selected_base[feature_dimension] = (
                selected_base[feature_dimension]
                .fillna("")
                .astype(str)
                .str.strip()
            )
            selected_features_base = selected_base[
                selected_base[feature_dimension].ne("")
            ].copy()

            feature_mode = st.radio(
                "Escolha das características",
                ["Mais frequentes", "Selecionar manualmente"],
                horizontal=True,
                key="feature_select_mode",
            )

            feature_counts = selected_features_base[feature_dimension].value_counts()

            if feature_mode == "Mais frequentes":
                max_features = max(1, min(200, len(feature_counts)))
                default_features = min(30, max_features)
                n_features = st.slider(
                    "Número de características",
                    min_value=1,
                    max_value=max_features,
                    value=default_features,
                    key="protein_compare_feature_n",
                )
                features_keep = feature_counts.head(n_features).index.tolist()
            else:
                features_keep = st.multiselect(
                    "Selecione produtos/genes/pistas",
                    feature_counts.index.tolist(),
                    max_selections=100,
                    key="protein_compare_features_manual",
                )

            matrix_source = selected_features_base[
                selected_features_base[feature_dimension].isin(features_keep)
            ]

            if len(matrix_source) and features_keep:
                count_matrix = pd.crosstab(
                    matrix_source[feature_dimension],
                    matrix_source[compare_dimension],
                )
                count_matrix = count_matrix.reindex(index=features_keep).fillna(0)

                if metric_mode == "Contagem de CDS":
                    matrix = count_matrix.copy()
                    color_label = "CDS"
                    color_scale = "Cividis"
                elif metric_mode == "% das CDS do táxon":
                    totals = selected_base.groupby(compare_dimension).size()
                    matrix = count_matrix.copy().astype(float)
                    for col in matrix.columns:
                        denom = totals.get(col, 0)
                        matrix[col] = 100 * matrix[col] / denom if denom else 0
                    color_label = "% das CDS"
                    color_scale = "Cividis"
                else:
                    matrix = count_matrix.gt(0).astype(int)
                    color_label = "Presença"
                    color_scale = [
                        [0.0, "#F2F2F2"],
                        [0.49, "#F2F2F2"],
                        [0.5, "#0072B2"],
                        [1.0, "#0072B2"],
                    ]

                fig = px.imshow(
                    matrix,
                    aspect="auto",
                    color_continuous_scale=color_scale,
                    labels={
                        "x": compare_dimension,
                        "y": feature_dimension,
                        "color": color_label,
                    },
                    title=f"{feature_dimension} × {compare_dimension}",
                )
                style_fig(
                    fig,
                    height=max(560, min(12000, 32 * len(matrix) + 260)),
                )
                fig.update_xaxes(side="top")
                st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

                matrix_export = matrix.reset_index()
                st.dataframe(matrix_export, hide_index=True, width="stretch")
                st.download_button(
                    "Baixar matriz",
                    download_tsv(matrix_export),
                    file_name="taxa_protein_matrix.tsv",
                    mime="text/tab-separated-values",
                )

            st.subheader("Perfil de anotação por táxon")
            class_by_taxon = (
                selected_base.groupby([compare_dimension, "protein_annotation_class"])
                .size()
                .reset_index(name="CDS")
            )
            if len(class_by_taxon):
                fig = px.bar(
                    class_by_taxon,
                    x=compare_dimension,
                    y="CDS",
                    color="protein_annotation_class",
                    barmode="stack",
                    color_discrete_map=ANNOTATION_CLASS_COLORS,
                    labels={"protein_annotation_class": "Classe de anotação"},
                    title="Composição das classes de anotação",
                )
                style_fig(fig, height=560, legend_title="Classe de anotação")
                fig.update_xaxes(tickangle=-30)
                st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

                totals = class_by_taxon.groupby(compare_dimension)["CDS"].transform("sum")
                class_pct = class_by_taxon.copy()
                class_pct["Percentual"] = 100 * class_pct["CDS"] / totals
                fig = px.bar(
                    class_pct,
                    x=compare_dimension,
                    y="Percentual",
                    color="protein_annotation_class",
                    barmode="stack",
                    color_discrete_map=ANNOTATION_CLASS_COLORS,
                    labels={
                        "protein_annotation_class": "Classe de anotação",
                        "Percentual": "% das CDS",
                    },
                    title="Composição percentual das classes de anotação",
                )
                style_fig(fig, height=560, legend_title="Classe de anotação")
                fig.update_yaxes(range=[0, 100])
                fig.update_xaxes(tickangle=-30)
                st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)
        else:
            st.info("Selecione pelo menos um táxon para iniciar a comparação.")


    # ============================================================
    # QUALIDADE
    # ============================================================

if tabs[10].open:
    with tabs[10]:
        st.header("Completude dos metadados")

        fig = px.bar(
            quality,
            x="field",
            y="available_pct",
            text="available_pct",
            color="available_pct",
            color_continuous_scale="Cividis",
            labels={
                "field": "Campo",
                "available_pct": "% disponível",
            },
            title="Completude dos principais campos",
        )
        style_fig(fig, height=520)
        fig.update_layout(coloraxis_colorbar=dict(title="% disponível"))
        st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        st.dataframe(quality, hide_index=True, width="stretch")

        if "taxonomy_final_status" in records.columns:
            st.subheader("Status taxonômico")
            tax_status = (
                records["taxonomy_final_status"]
                .fillna("")
                .astype(str)
                .str.strip()
                .replace("", "Missing")
                .value_counts()
                .rename_axis("Status")
                .reset_index(name="Registros")
            )
            fig = px.bar(
                tax_status,
                x="Registros",
                y="Status",
                orientation="h",
                color="Status",
                color_discrete_sequence=ACCESSIBLE_COLORS,
            )
            style_fig(fig, categorical_height(len(tax_status), minimum=420))
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)
            st.dataframe(tax_status, hide_index=True, width="stretch")


    # ============================================================
    # DPOL
    # ============================================================

if tabs[11].open:
    with tabs[11]:
        st.header("Fluxo DPOL / curadoria")

        flow_clean = dpol_flow.dropna(subset=["n"]).copy()
        flow_clean["n"] = pd.to_numeric(flow_clean["n"], errors="coerce")
        flow_clean = flow_clean.dropna(subset=["n"])

        if len(flow_clean):
            funnel_colors = [
                ACCESSIBLE_COLORS[i % len(ACCESSIBLE_COLORS)]
                for i in range(len(flow_clean))
            ]
            fig = go.Figure(
                go.Funnel(
                    y=flow_clean["stage"],
                    x=flow_clean["n"],
                    textinfo="value",
                    marker=dict(color=funnel_colors),
                )
            )
            style_fig(fig, height=600)
            st.plotly_chart(fig, width="stretch", config=PLOTLY_CONFIG)

        st.dataframe(dpol_flow, hide_index=True, width="stretch")


    # ============================================================
    # REGISTROS
    # ============================================================

if tabs[12].open:
    with tabs[12]:
        st.header("Tabela geral")

        general_query = st.text_input(
            "Pesquisar nos registros",
            placeholder="accession, organismo, espécie, país, hospedeiro...",
            key="general_query",
        )

        table = filtered_records.copy()

        if general_query:
            pattern = re.escape(general_query)
            mask = pd.Series(False, index=table.index)
            for column in available_columns(
                table,
                [
                    "accession_version",
                    "definition",
                    "organism",
                    "species_final",
                    "interest_group",
                    "country_final",
                    "host",
                    "isolate",
                    "strain",
                ],
            ):
                mask |= (
                    table[column]
                    .fillna("")
                    .astype(str)
                    .str.contains(pattern, case=False, regex=True, na=False)
                )
            table = table[mask]

        display = available_columns(
            table,
            [
                "accession_version",
                "taxid",
                "definition",
                "length_nt",
                "sequence_status",
                "organism",
                "species_final",
                "interest_group",
                "taxonomy_final_status",
                "country_final",
                "country_source",
                "collection_date",
                "collection_year",
                "geo_loc_name",
                "host",
                "lab_host",
                "isolate",
                "strain",
            ],
        )

        st.write(f"Registros encontrados: **{len(table):,}**")
        record_preview_n = st.selectbox(
            "Linhas na prévia",
            [1000, 5000, 10000, 20000, 50000],
            index=0,
            key="record_preview_n",
        )
        st.dataframe(
            table[display].head(record_preview_n),
            hide_index=True,
            width="stretch",
        )

        if len(table) > record_preview_n:
            st.caption(
                f"A tela mostra {record_preview_n:,} linhas; o download contém todas "
                f"as {len(table):,} linhas filtradas."
            )

        download_on_demand(
            "Baixar tabela filtrada completa",
            table[display],
            "Orthoherpes_records_filtered.tsv",
            "dl_records_full",
        )


# ============================================================
# RODAPÉ
# ============================================================

st.divider()
st.caption(
    "Representação no dataset corresponde ao número de registros NCBI/GenBank "
    "e não deve ser interpretada como prevalência biológica. species_final "
    "preserva a classificação normalizada; organism preserva a nomenclatura "
    "original do registro. Paleta de cores escolhida para melhor distinção em "
    "públicos com deficiência de visão de cores."
)
