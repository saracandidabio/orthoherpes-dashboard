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


# Execução local:
# ROOT/DASHBOARD_DATASET_INICIAL/dashboard_dataset_inicial.py
#
# Community Cloud:
# ROOT/streamlit_app.py

if (
    HERE
    / "DASHBOARD_DATASET_INICIAL"
    / "data"
).exists():

    PROJECT_ROOT = HERE

    DATA_DIR = (
        HERE
        / "DASHBOARD_DATASET_INICIAL"
        / "data"
    )

else:

    PROJECT_ROOT = HERE.parent

    DATA_DIR = (
        HERE
        / "data"
    )


RECORDS_FILE = (
    DATA_DIR
    / "records_dashboard.tsv.gz"
)

CDS_FILE = (
    DATA_DIR
    / "cds_dashboard.tsv.gz"
)

INTEREST_FILE = (
    DATA_DIR
    / "interest_summary.tsv"
)

QUALITY_FILE = (
    DATA_DIR
    / "metadata_quality.tsv"
)

DPOL_FILE = (
    DATA_DIR
    / "dpol_flow.tsv"
)


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


# ============================================================
# FUNÇÕES
# ============================================================

def available_columns(
    df,
    columns,
):

    return [
        column
        for column in columns
        if column in df.columns
    ]


def normalized(series):

    return (
        series
        .fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
    )


def missing_mask(series):

    return normalized(
        series
    ).isin(
        MISSING_TERMS
    )


def download_tsv(df):

    return (
        df.to_csv(
            sep="\t",
            index=False,
        )
        .encode(
            "utf-8"
        )
    )


def filter_dataframe(
    df,
    countries,
    species,
    interest,
    year_range,
    include_undated,
):

    out = df.copy()


    if (
        countries
        and
        "country_final"
        in out.columns
    ):

        out = out[
            out[
                "country_final"
            ]
            .isin(
                countries
            )
        ]


    if (
        species
        and
        "species_final"
        in out.columns
    ):

        out = out[
            out[
                "species_final"
            ]
            .isin(
                species
            )
        ]


    if (
        interest
        and
        "interest_group"
        in out.columns
    ):

        out = out[
            out[
                "interest_group"
            ]
            .isin(
                interest
            )
        ]


    if (
        year_range is not None
        and
        "collection_year_num"
        in out.columns
    ):

        years = pd.to_numeric(
            out[
                "collection_year_num"
            ],
            errors="coerce",
        )


        mask = years.between(
            year_range[0],
            year_range[1],
            inclusive="both",
        )


        if include_undated:

            mask |= years.isna()


        out = out[
            mask
        ]


    return out


# ============================================================
# VALIDAR DADOS
# ============================================================

for required in [
    RECORDS_FILE,
    CDS_FILE,
    QUALITY_FILE,
    DPOL_FILE,
]:

    if not required.exists():

        st.error(
            "Base do dashboard não encontrada:\n\n"
            f"`{required}`\n\n"
            "Execute primeiro "
            "`01_build_dashboard_data.py`."
        )

        st.stop()


# ============================================================
# CARREGAMENTO
# ============================================================

@st.cache_data(
    show_spinner=(
        "Carregando registros..."
    )
)
def load_records():

    df = pd.read_csv(
        RECORDS_FILE,
        sep="\t",
        compression="gzip",
        low_memory=False,
    )

    if (
        "collection_year_num"
        in df.columns
    ):

        df[
            "collection_year_num"
        ] = pd.to_numeric(
            df[
                "collection_year_num"
            ],
            errors="coerce",
        )

    return df


@st.cache_data(
    show_spinner=(
        "Carregando CDS..."
    )
)
def load_cds():

    df = pd.read_csv(
        CDS_FILE,
        sep="\t",
        compression="gzip",
        low_memory=False,
    )

    if (
        "collection_year_num"
        in df.columns
    ):

        df[
            "collection_year_num"
        ] = pd.to_numeric(
            df[
                "collection_year_num"
            ],
            errors="coerce",
        )

    return df


@st.cache_data
def load_simple_table(
    path
):

    return pd.read_csv(
        path,
        sep="\t",
    )


records = load_records()
cds = load_cds()

quality = load_simple_table(
    QUALITY_FILE
)

dpol_flow = load_simple_table(
    DPOL_FILE
)


# ============================================================
# CABEÇALHO
# ============================================================

st.title(
    "🧬 Orthoherpesviridae Dataset Explorer"
)

st.caption(
    "Exploração do universo inicial do NCBI/GenBank, "
    "taxonomia normalizada, metadados, hospedeiros, "
    "genes, CDS, proteínas hipotéticas e fluxo de curadoria DPOL."
)


# ============================================================
# SIDEBAR
# ============================================================

st.sidebar.header(
    "Filtros globais"
)


country_options = sorted(
    records.loc[
        records[
            "country_final"
        ]
        .fillna("")
        .ne("Unknown"),
        "country_final"
    ]
    .dropna()
    .unique()
    .tolist()
)


species_options = sorted(
    records.loc[
        records[
            "species_final"
        ]
        .fillna("")
        .ne(
            "Unresolved_species"
        ),
        "species_final"
    ]
    .dropna()
    .unique()
    .tolist()
)


interest_options = [
    value
    for value
    in INTEREST_ORDER
    if value
    in set(
        records[
            "interest_group"
        ]
        .dropna()
    )
]


selected_interest = (
    st.sidebar.multiselect(
        "Espécies / vírus de interesse",
        interest_options,
    )
)


selected_countries = (
    st.sidebar.multiselect(
        "País",
        country_options,
    )
)


selected_species = (
    st.sidebar.multiselect(
        "species_final",
        species_options,
    )
)


year_values = (
    records[
        "collection_year_num"
    ]
    .dropna()
)


if len(year_values):

    min_year = int(
        year_values.min()
    )

    max_year = int(
        year_values.max()
    )

    selected_year = (
        st.sidebar.slider(
            "Ano de coleta",
            min_value=min_year,
            max_value=max_year,
            value=(
                min_year,
                max_year,
            ),
        )
    )

else:

    selected_year = None


include_undated = (
    st.sidebar.checkbox(
        "Incluir registros sem ano",
        value=True,
    )
)


filtered_records = (
    filter_dataframe(
        records,
        selected_countries,
        selected_species,
        selected_interest,
        selected_year,
        include_undated,
    )
)


filtered_cds = (
    filter_dataframe(
        cds,
        selected_countries,
        selected_species,
        selected_interest,
        selected_year,
        include_undated,
    )
)


# ============================================================
# MÉTRICAS GLOBAIS
# ============================================================

m1, m2, m3, m4 = (
    st.columns(4)
)


m1.metric(
    "Registros no dataset",
    f"{len(records):,}"
)


m2.metric(
    "Registros após filtros",
    f"{len(filtered_records):,}"
)


m3.metric(
    "CDS no dataset",
    f"{len(cds):,}"
)


m4.metric(
    "CDS após filtros",
    f"{len(filtered_cds):,}"
)


# ============================================================
# ABAS
# ============================================================

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
])


# ============================================================
# VISÃO GERAL
# ============================================================

with tabs[0]:

    st.header(
        "Visão geral"
    )


    c1, c2, c3, c4 = (
        st.columns(4)
    )


    c1.metric(
        "Espécies normalizadas",
        f"{records['species_final'].nunique():,}"
    )


    c2.metric(
        "Países conhecidos",
        f"{records.loc[records['country_final'] != 'Unknown', 'country_final'].nunique():,}"
    )


    c3.metric(
        "Registros brasileiros",
        f"{(records['country_final'] == 'Brazil').sum():,}"
    )


    c4.metric(
        "Proteínas hipotéticas",
        f"{(cds['is_hypothetical'] == 'YES').sum():,}"
    )


    st.subheader(
        "Espécies mais representadas"
    )


    top_species = (
        filtered_records[
            "species_final"
        ]
        .value_counts()
        .head(40)
        .sort_values()
        .rename_axis(
            "Espécie"
        )
        .reset_index(
            name="Registros"
        )
    )


    if len(top_species):

        st.plotly_chart(
            px.bar(
                top_species,
                x="Registros",
                y="Espécie",
                orientation="h",
            ),
            width="stretch",
        )


    st.subheader(
        "Espécies / vírus destacados"
    )


    highlighted = (
        filtered_records[
            filtered_records[
                "interest_group"
            ]
            .ne("Other")
        ]
    )


    interest_counts = (
        highlighted[
            "interest_group"
        ]
        .value_counts()
        .sort_values()
        .rename_axis(
            "Grupo"
        )
        .reset_index(
            name="Registros"
        )
    )


    if len(interest_counts):

        st.plotly_chart(
            px.bar(
                interest_counts,
                x="Registros",
                y="Grupo",
                orientation="h",
                text="Registros",
            ),
            width="stretch",
        )


# ============================================================
# ESPÉCIES DE INTERESSE
# ============================================================

with tabs[1]:

    st.header(
        "Espécies e vírus de interesse"
    )

    st.caption(
        "Esta categoria é uma camada de visualização. "
        "A taxonomia normalizada continua preservada em species_final."
    )


    focus = records[
        records[
            "interest_group"
        ]
        .ne("Other")
    ].copy()


    counts = (
        focus[
            "interest_group"
        ]
        .value_counts()
    )


    metric_cols = (
        st.columns(4)
    )


    for index, label in enumerate(
        INTEREST_ORDER
    ):

        metric_cols[
            index % 4
        ].metric(
            label,
            f"{int(counts.get(label, 0)):,}"
        )


    selected_compare = (
        st.multiselect(
            "Comparar",
            INTEREST_ORDER,
            default=[
                label
                for label
                in INTEREST_ORDER
                if counts.get(
                    label,
                    0
                ) > 0
            ],
            key="interest_compare",
        )
    )


    compare = focus[
        focus[
            "interest_group"
        ]
        .isin(
            selected_compare
        )
    ]


    if len(compare):

        st.subheader(
            "País"
        )


        tmp = (
            compare[
                compare[
                    "country_final"
                ]
                .ne("Unknown")
            ]
            .groupby(
                [
                    "country_final",
                    "interest_group",
                ]
            )
            .size()
            .reset_index(
                name="Registros"
            )
        )


        st.plotly_chart(
            px.bar(
                tmp,
                x="country_final",
                y="Registros",
                color="interest_group",
                barmode="group",
                labels={
                    "country_final":
                        "País",
                    "interest_group":
                        "Espécie / vírus",
                },
            ),
            width="stretch",
        )


        st.subheader(
            "Ano"
        )


        tmp = (
            compare[
                compare[
                    "collection_year_num"
                ]
                .notna()
            ]
            .groupby(
                [
                    "collection_year_num",
                    "interest_group",
                ]
            )
            .size()
            .reset_index(
                name="Registros"
            )
        )


        if len(tmp):

            st.plotly_chart(
                px.line(
                    tmp,
                    x="collection_year_num",
                    y="Registros",
                    color="interest_group",
                    markers=True,
                    labels={
                        "collection_year_num":
                            "Ano",
                        "interest_group":
                            "Espécie / vírus",
                    },
                ),
                width="stretch",
            )


        st.subheader(
            "Nomes originais no GenBank"
        )


        original_names = (
            compare
            .groupby(
                [
                    "interest_group",
                    "organism",
                ]
            )
            .size()
            .reset_index(
                name="Registros"
            )
            .sort_values(
                [
                    "interest_group",
                    "Registros",
                ],
                ascending=[
                    True,
                    False,
                ],
            )
        )


        st.dataframe(
            original_names,
            hide_index=True,
            width="stretch",
        )


# ============================================================
# EXPLORADOR TAXONÔMICO
# ============================================================

with tabs[2]:

    st.header(
        "Explorador taxonômico"
    )


    explorer_mode = st.radio(
        "Usar",
        [
            "species_final",
            "organism",
        ],
        horizontal=True,
    )


    explorer_options = sorted(
        records[
            explorer_mode
        ]
        .dropna()
        .astype(str)
        .loc[
            lambda x:
                x.str.strip().ne("")
        ]
        .unique()
        .tolist()
    )


    explorer_selected = (
        st.multiselect(
            "Selecione um ou mais táxons",
            explorer_options,
            key="explorer_taxa",
        )
    )


    if explorer_selected:

        explorer = records[
            records[
                explorer_mode
            ]
            .isin(
                explorer_selected
            )
        ]


        e1, e2, e3, e4 = (
            st.columns(4)
        )


        e1.metric(
            "Registros",
            f"{len(explorer):,}"
        )


        e2.metric(
            "Países",
            f"{explorer.loc[explorer['country_final'] != 'Unknown', 'country_final'].nunique():,}"
        )


        e3.metric(
            "Hospedeiros",
            f"{explorer.loc[~missing_mask(explorer['host']), 'host'].nunique():,}"
        )


        e4.metric(
            "Anos",
            f"{explorer['collection_year_num'].dropna().nunique():,}"
        )


        country = (
            explorer[
                explorer[
                    "country_final"
                ]
                .ne("Unknown")
            ][
                "country_final"
            ]
            .value_counts()
            .head(40)
            .sort_values()
            .rename_axis(
                "País"
            )
            .reset_index(
                name="Registros"
            )
        )


        if len(country):

            st.plotly_chart(
                px.bar(
                    country,
                    x="Registros",
                    y="País",
                    orientation="h",
                ),
                width="stretch",
            )


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


        st.dataframe(
            explorer[
                display
            ],
            hide_index=True,
            width="stretch",
        )


        st.download_button(
            "Baixar registros selecionados",
            download_tsv(
                explorer[
                    display
                ]
            ),
            file_name=(
                "taxon_explorer.tsv"
            ),
            mime=(
                "text/tab-separated-values"
            ),
        )


# ============================================================
# GEOGRAFIA
# ============================================================

with tabs[3]:

    st.header(
        "Distribuição geográfica"
    )


    geo = filtered_records[
        filtered_records[
            "country_final"
        ]
        .ne("Unknown")
    ]


    countries = (
        geo[
            "country_final"
        ]
        .value_counts()
        .head(50)
        .sort_values()
        .rename_axis(
            "País"
        )
        .reset_index(
            name="Registros"
        )
    )


    if len(countries):

        st.plotly_chart(
            px.bar(
                countries,
                x="Registros",
                y="País",
                orientation="h",
            ),
            width="stretch",
        )


    diversity = (
        geo[
            geo[
                "species_final"
            ]
            .ne(
                "Unresolved_species"
            )
        ]
        .groupby(
            "country_final"
        )[
            "species_final"
        ]
        .nunique()
        .sort_values(
            ascending=False
        )
        .head(50)
        .sort_values()
        .rename_axis(
            "País"
        )
        .reset_index(
            name="Espécies distintas"
        )
    )


    st.subheader(
        "Diversidade taxonômica"
    )


    if len(diversity):

        st.plotly_chart(
            px.bar(
                diversity,
                x="Espécies distintas",
                y="País",
                orientation="h",
            ),
            width="stretch",
        )


# ============================================================
# TEMPO
# ============================================================

with tabs[4]:

    st.header(
        "Distribuição temporal"
    )


    dated = filtered_records[
        filtered_records[
            "collection_year_num"
        ]
        .notna()
    ]


    annual = (
        dated
        .groupby(
            "collection_year_num"
        )
        .agg(
            Registros=(
                "accession_version",
                "count",
            ),
            Espécies=(
                "species_final",
                "nunique",
            ),
        )
        .reset_index()
        .sort_values(
            "collection_year_num"
        )
    )


    if len(annual):

        st.plotly_chart(
            px.line(
                annual,
                x="collection_year_num",
                y="Registros",
                markers=True,
                labels={
                    "collection_year_num":
                        "Ano"
                },
            ),
            width="stretch",
        )


        st.plotly_chart(
            px.line(
                annual,
                x="collection_year_num",
                y="Espécies",
                markers=True,
                labels={
                    "collection_year_num":
                        "Ano"
                },
            ),
            width="stretch",
        )


# ============================================================
# BRASIL
# ============================================================

with tabs[5]:

    st.header(
        "Brasil"
    )


    brazil = records[
        records[
            "country_final"
        ]
        .eq("Brazil")
    ]


    b1, b2, b3 = (
        st.columns(3)
    )


    b1.metric(
        "Registros",
        f"{len(brazil):,}"
    )


    b2.metric(
        "Espécies",
        f"{brazil['species_final'].nunique():,}"
    )


    b3.metric(
        "Hospedeiros",
        f"{brazil.loc[~missing_mask(brazil['host']), 'host'].nunique():,}"
    )


    br_species = (
        brazil[
            "species_final"
        ]
        .value_counts()
        .head(50)
        .sort_values()
        .rename_axis(
            "Espécie"
        )
        .reset_index(
            name="Registros"
        )
    )


    if len(br_species):

        st.plotly_chart(
            px.bar(
                br_species,
                x="Registros",
                y="Espécie",
                orientation="h",
            ),
            width="stretch",
        )


    br_interest = (
        brazil[
            brazil[
                "interest_group"
            ]
            .ne("Other")
        ][
            "interest_group"
        ]
        .value_counts()
        .sort_values()
        .rename_axis(
            "Grupo"
        )
        .reset_index(
            name="Registros"
        )
    )


    if len(br_interest):

        st.subheader(
            "Espécies de interesse no Brasil"
        )


        st.plotly_chart(
            px.bar(
                br_interest,
                x="Registros",
                y="Grupo",
                orientation="h",
                text="Registros",
            ),
            width="stretch",
        )


# ============================================================
# HOSPEDEIROS
# ============================================================

with tabs[6]:

    st.header(
        "Hospedeiros"
    )


    host_records = (
        filtered_records[
            ~missing_mask(
                filtered_records[
                    "host"
                ]
            )
        ]
    )


    hosts = (
        host_records[
            "host"
        ]
        .value_counts()
        .head(60)
        .sort_values()
        .rename_axis(
            "Hospedeiro"
        )
        .reset_index(
            name="Registros"
        )
    )


    if len(hosts):

        st.plotly_chart(
            px.bar(
                hosts,
                x="Registros",
                y="Hospedeiro",
                orientation="h",
            ),
            width="stretch",
        )


    selected_host = (
        st.selectbox(
            "Investigar hospedeiro",
            [""]
            +
            sorted(
                host_records[
                    "host"
                ]
                .dropna()
                .unique()
                .tolist()
            ),
        )
    )


    if selected_host:

        host_species = (
            host_records[
                host_records[
                    "host"
                ]
                .eq(
                    selected_host
                )
            ][
                "species_final"
            ]
            .value_counts()
            .head(60)
            .sort_values()
            .rename_axis(
                "Espécie"
            )
            .reset_index(
                name="Registros"
            )
        )


        st.plotly_chart(
            px.bar(
                host_species,
                x="Registros",
                y="Espécie",
                orientation="h",
            ),
            width="stretch",
        )


# ============================================================
# GENES / CDS
# ============================================================

with tabs[7]:

    st.header(
        "Genes e CDS"
    )


    gene_values = (
        filtered_cds.loc[
            ~missing_mask(
                filtered_cds[
                    "gene"
                ]
            ),
            "gene"
        ]
        if "gene"
        in filtered_cds.columns
        else pd.Series(
            dtype=str
        )
    )


    genes = (
        gene_values
        .value_counts()
        .head(60)
        .sort_values()
        .rename_axis(
            "Gene"
        )
        .reset_index(
            name="CDS"
        )
    )


    if len(genes):

        st.subheader(
            "Genes mais representados"
        )


        st.plotly_chart(
            px.bar(
                genes,
                x="CDS",
                y="Gene",
                orientation="h",
            ),
            width="stretch",
        )


    if "product" in filtered_cds.columns:

        products = (
            filtered_cds.loc[
                ~missing_mask(
                    filtered_cds[
                        "product"
                    ]
                ),
                "product"
            ]
            .value_counts()
            .head(60)
            .sort_values()
            .rename_axis(
                "Produto"
            )
            .reset_index(
                name="CDS"
            )
        )


        if len(products):

            st.subheader(
                "Produtos mais representados"
            )


            st.plotly_chart(
                px.bar(
                    products,
                    x="CDS",
                    y="Produto",
                    orientation="h",
                ),
                width="stretch",
            )


    st.subheader(
        "Busca nas anotações"
    )


    query = st.text_input(
        "Pesquisar",
        placeholder=(
            "polymerase, capsid, tegument, ORF, kinase..."
        ),
        key="cds_search",
    )


    cds_search_fields = (
        available_columns(
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
            ],
        )
    )


    search_hits = (
        filtered_cds.copy()
    )


    if query:

        pattern = re.escape(
            query
        )

        mask = pd.Series(
            False,
            index=search_hits.index,
        )


        for column in (
            cds_search_fields
        ):

            mask |= (
                search_hits[
                    column
                ]
                .fillna("")
                .astype(str)
                .str.contains(
                    pattern,
                    case=False,
                    regex=True,
                    na=False,
                )
            )


        search_hits = (
            search_hits[
                mask
            ]
        )


    st.write(
        f"CDS exibidas: "
        f"**{len(search_hits):,}**"
    )


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
            "protein_annotation_class",
        ],
    )


    st.dataframe(
        search_hits[
            display
        ]
        .head(10000),
        hide_index=True,
        width="stretch",
    )


    st.download_button(
        "Baixar resultados",
        download_tsv(
            search_hits[
                display
            ]
        ),
        file_name="CDS_search.tsv",
        mime=(
            "text/tab-separated-values"
        ),
    )


# ============================================================
# PROTEÍNAS HIPOTÉTICAS
# ============================================================

with tabs[8]:

    st.header(
        "Proteínas hipotéticas"
    )


    st.info(
        "Hypothetical protein indica uma CDS prevista "
        "sem atribuição funcional específica suficientemente "
        "estabelecida na anotação. Este painel procura também "
        "pistas adicionais presentes em gene, standard_name, "
        "function, note e inference."
    )


    annotation_classes = (
        filtered_cds[
            "protein_annotation_class"
        ]
        .value_counts()
        .rename_axis(
            "Classe"
        )
        .reset_index(
            name="CDS"
        )
    )


    st.plotly_chart(
        px.bar(
            annotation_classes,
            x="Classe",
            y="CDS",
            text="CDS",
        ),
        width="stretch",
    )


    hyp = filtered_cds[
        filtered_cds[
            "is_hypothetical"
        ]
        .eq("YES")
    ].copy()


    h1, h2, h3, h4 = (
        st.columns(4)
    )


    h1.metric(
        "Hypothetical",
        f"{len(hyp):,}"
    )


    h2.metric(
        "Com pista de anotação",
        f"{(hyp['protein_annotation_class'] == 'Hypothetical with annotation clue').sum():,}"
    )


    h3.metric(
        "Somente identificador",
        f"{(hyp['protein_annotation_class'] == 'Hypothetical with identifier only').sum():,}"
    )


    h4.metric(
        "Sem informação adicional",
        f"{(hyp['protein_annotation_class'] == 'Hypothetical without additional annotation').sum():,}"
    )


    class_filter = (
        st.multiselect(
            "Categoria",
            sorted(
                hyp[
                    "protein_annotation_class"
                ]
                .dropna()
                .unique()
                .tolist()
            ),
        )
    )


    if class_filter:

        hyp = hyp[
            hyp[
                "protein_annotation_class"
            ]
            .isin(
                class_filter
            )
        ]


    hyp_query = (
        st.text_input(
            "Pesquisar nas proteínas hipotéticas",
            placeholder=(
                "ORF, tegument, capsid, locus tag, protein ID..."
            ),
            key="hyp_query",
        )
    )


    if hyp_query:

        pattern = re.escape(
            hyp_query
        )

        mask = pd.Series(
            False,
            index=hyp.index,
        )


        for column in (
            available_columns(
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
            )
        ):

            mask |= (
                hyp[
                    column
                ]
                .fillna("")
                .astype(str)
                .str.contains(
                    pattern,
                    case=False,
                    regex=True,
                    na=False,
                )
            )


        hyp = hyp[
            mask
        ]


    hyp_display = (
        available_columns(
            hyp,
            [
                "accession_version",
                "organism",
                "species_final",
                "country_final",
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
                "protein_annotation_class",
            ],
        )
    )


    st.dataframe(
        hyp[
            hyp_display
        ]
        .head(15000),
        hide_index=True,
        width="stretch",
    )


    st.download_button(
        "Baixar proteínas hipotéticas filtradas",
        download_tsv(
            hyp[
                hyp_display
            ]
        ),
        file_name=(
            "hypothetical_proteins.tsv"
        ),
        mime=(
            "text/tab-separated-values"
        ),
    )


# ============================================================
# COMPARAÇÃO DE PROTEÍNAS
# ============================================================

with tabs[9]:

    st.header(
        "Comparação entre táxons e produtos"
    )


    compare_dimension = (
        st.radio(
            "Comparar por",
            [
                "species_final",
                "organism",
            ],
            horizontal=True,
            key="protein_compare_dimension",
        )
    )


    compare_options = sorted(
        filtered_cds[
            compare_dimension
        ]
        .dropna()
        .astype(str)
        .loc[
            lambda x:
                x.str.strip().ne("")
        ]
        .unique()
        .tolist()
    )


    compare_taxa = (
        st.multiselect(
            "Selecione os táxons",
            compare_options,
            max_selections=12,
            key="protein_compare_taxa",
        )
    )


    matrix_mode = (
        st.radio(
            "Métrica",
            [
                "Contagem de CDS",
                "Presença / ausência",
            ],
            horizontal=True,
        )
    )


    top_products_n = (
        st.slider(
            "Número máximo de produtos",
            10,
            100,
            30,
            5,
        )
    )


    if compare_taxa:

        matrix_source = (
            filtered_cds[
                filtered_cds[
                    compare_dimension
                ]
                .isin(
                    compare_taxa
                )
                &
                ~missing_mask(
                    filtered_cds[
                        "product"
                    ]
                )
            ]
        )


        top_products = (
            matrix_source[
                "product"
            ]
            .value_counts()
            .head(
                top_products_n
            )
            .index
        )


        matrix_source = (
            matrix_source[
                matrix_source[
                    "product"
                ]
                .isin(
                    top_products
                )
            ]
        )


        matrix = pd.crosstab(
            matrix_source[
                "product"
            ],
            matrix_source[
                compare_dimension
            ],
        )


        if (
            matrix_mode
            == "Presença / ausência"
        ):

            matrix = (
                matrix.gt(0)
                .astype(int)
            )


        if not matrix.empty:

            fig = px.imshow(
                matrix,
                aspect="auto",
                labels={
                    "x":
                        compare_dimension,
                    "y":
                        "Produto",
                    "color":
                        matrix_mode,
                },
            )


            st.plotly_chart(
                fig,
                width="stretch",
            )


            matrix_export = (
                matrix
                .reset_index()
            )


            st.dataframe(
                matrix_export,
                hide_index=True,
                width="stretch",
            )


            st.download_button(
                "Baixar matriz",
                download_tsv(
                    matrix_export
                ),
                file_name=(
                    "taxa_protein_matrix.tsv"
                ),
                mime=(
                    "text/tab-separated-values"
                ),
            )


# ============================================================
# QUALIDADE
# ============================================================

with tabs[10]:

    st.header(
        "Completude dos metadados"
    )


    st.plotly_chart(
        px.bar(
            quality,
            x="field",
            y="available_pct",
            text="available_pct",
            labels={
                "field":
                    "Campo",
                "available_pct":
                    "% disponível",
            },
        ),
        width="stretch",
    )


    st.dataframe(
        quality,
        hide_index=True,
        width="stretch",
    )


    if (
        "taxonomy_final_status"
        in records.columns
    ):

        st.subheader(
            "Status taxonômico"
        )


        tax_status = (
            records[
                "taxonomy_final_status"
            ]
            .value_counts()
            .rename_axis(
                "Status"
            )
            .reset_index(
                name="Registros"
            )
        )


        st.dataframe(
            tax_status,
            hide_index=True,
            width="stretch",
        )


# ============================================================
# DPOL
# ============================================================

with tabs[11]:

    st.header(
        "Fluxo DPOL / curadoria"
    )


    flow_clean = (
        dpol_flow
        .dropna(
            subset=[
                "n"
            ]
        )
        .copy()
    )


    flow_clean[
        "n"
    ] = pd.to_numeric(
        flow_clean[
            "n"
        ],
        errors="coerce",
    )


    flow_clean = (
        flow_clean
        .dropna(
            subset=[
                "n"
            ]
        )
    )


    if len(flow_clean):

        fig = go.Figure(
            go.Funnel(
                y=flow_clean[
                    "stage"
                ],
                x=flow_clean[
                    "n"
                ],
                textinfo=(
                    "value"
                ),
            )
        )


        st.plotly_chart(
            fig,
            width="stretch",
        )


    st.dataframe(
        dpol_flow,
        hide_index=True,
        width="stretch",
    )


# ============================================================
# REGISTROS
# ============================================================

with tabs[12]:

    st.header(
        "Tabela geral"
    )


    general_query = (
        st.text_input(
            "Pesquisar nos registros",
            placeholder=(
                "accession, organismo, espécie, país, hospedeiro..."
            ),
            key="general_query",
        )
    )


    table = (
        filtered_records.copy()
    )


    if general_query:

        pattern = re.escape(
            general_query
        )

        mask = pd.Series(
            False,
            index=table.index,
        )


        for column in (
            available_columns(
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
            )
        ):

            mask |= (
                table[
                    column
                ]
                .fillna("")
                .astype(str)
                .str.contains(
                    pattern,
                    case=False,
                    regex=True,
                    na=False,
                )
            )


        table = table[
            mask
        ]


    display = available_columns(
        table,
        [
            "accession_version",
            "taxid",
            "organism",
            "species_final",
            "interest_group",
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


    st.write(
        f"Registros encontrados: "
        f"**{len(table):,}**"
    )


    st.dataframe(
        table[
            display
        ]
        .head(20000),
        hide_index=True,
        width="stretch",
    )


    if len(table) > 20000:

        st.caption(
            "A tela mostra no máximo 20.000 linhas. "
            "O download abaixo contém todas as linhas filtradas."
        )


    st.download_button(
        "Baixar tabela filtrada",
        download_tsv(
            table[
                display
            ]
        ),
        file_name=(
            "Orthoherpes_records_filtered.tsv"
        ),
        mime=(
            "text/tab-separated-values"
        ),
    )


# ============================================================
# RODAPÉ
# ============================================================

st.divider()

st.caption(
    "Representação no dataset corresponde ao número de registros "
    "NCBI/GenBank e não deve ser interpretada como prevalência biológica. "
    "species_final preserva a classificação normalizada; organism preserva "
    "a nomenclatura original do registro."
)
