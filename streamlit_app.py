"""
Orthoherpesviridae Dataset Explorer

Estrutura:
  dashboard_lib.py  -> lógica de dados (pandas puro, testável com pytest)
  streamlit_app.py  -> interface (este arquivo)

Cada seção é uma página do st.navigation: só a seção aberta é executada a cada
interação (antes, as 13 abas rodavam sempre, mesmo as invisíveis).

REGRA: os DataFrames vêm de st.cache_resource e são compartilhados entre
sessões. Nunca altere um DataFrame recebido; sempre derive um novo.
"""

from datetime import datetime
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

import dashboard_lib as L

APP_VERSION = "8.1.0"
from dashboard_lib import (
    INTEREST_ORDER, INTEREST_OTHER, CLUE_CLASS, Filters,
    available_columns, categorical_height, decat, fmt_int, fmt_pct,
    group_count, is_missing, known, value_counts_table,
)

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

# Servidor:  ROOT/DASHBOARD_DATASET_INICIAL/<este arquivo>
# Cloud:     ROOT/streamlit_app.py
if (HERE / "DASHBOARD_DATASET_INICIAL" / "data").exists():
    DATA_DIR = HERE / "DASHBOARD_DATASET_INICIAL" / "data"
else:
    DATA_DIR = HERE / "data"

DPOL_FILE = DATA_DIR / "dpol_flow.tsv"
INTEREST_FILE = DATA_DIR / "interest_summary.tsv"
QUALITY_FILE = DATA_DIR / "metadata_quality.tsv"
BUILD_AUDIT_FILE = DATA_DIR / "BUILD_AUDIT.txt"
MANIFEST_FILE = DATA_DIR / "dataset_manifest_v8_1.json"
GLOBAL_ANNOTATION_FILE = DATA_DIR / "cds_global_annotation_class.tsv"
GLOBAL_GENE_FILE = DATA_DIR / "cds_global_gene_counts.tsv.gz"
GLOBAL_PRODUCT_FILE = DATA_DIR / "cds_global_product_counts.tsv.gz"
GLOBAL_CLUE_FILE = DATA_DIR / "cds_global_clue_counts.tsv.gz"
GLOBAL_CLUE_SOURCE_FILE = DATA_DIR / "cds_global_clue_source_counts.tsv"
GLOBAL_TAXON_FILE = DATA_DIR / "cds_global_taxon_summary.tsv.gz"
LOCAL_VALIDATION_FILE = DATA_DIR / "local_validation_summary.tsv"
NCBI_VALIDATION_FILE = DATA_DIR / "ncbi_validation_summary.tsv"
COUNTRY_EXCEPTIONS_FILE = DATA_DIR / "country_plot_exceptions.tsv"
HOST_CROSSWALK_FILE = DATA_DIR / L.HOST_CROSSWALK_FILENAME
HOST_NORMALIZATION_AUDIT_FILE = DATA_DIR / "HOST_NORMALIZATION_AUDIT.tsv"
COUNTRY_NORMALIZATION_AUDIT_FILE = DATA_DIR / "COUNTRY_NORMALIZATION_AUDIT.tsv"
COUNTRY_CROSSWALK_FILE = DATA_DIR / "country_crosswalk.tsv"
GENE_CONTEXT_FILE = DATA_DIR / "gene_context_legend.tsv.gz"

# Paleta Okabe-Ito (amigável a daltonismo)
COLORS = [
    "#0072B2", "#E69F00", "#009E73", "#CC79A7",
    "#D55E00", "#56B4E9", "#F0E442", "#8A5CF6",
]
INTEREST_COLORS = {
    "CalHV-3": "#0072B2",
    "Cebus albifrons lymphocryptovirus 1": "#E69F00",
    "Cebus apella lymphocryptovirus 1": "#009E73",
    "Callithrix penicillata lymphocryptovirus 1": "#CC79A7",
    "Alouatta macconnelli cytomegalovirus": "#D55E00",
    "Alouatta palliata cytomegalovirus": "#56B4E9",
    "Alouatta seniculus cytomegalovirus": "#F0E442",
    "HSV-1": "#8A5CF6",
    "Other": "#8B95A1",
}
ANNOTATION_CLASS_COLORS = {
    "Annotated product": "#0072B2",
    "Hypothetical with annotation clue": "#009E73",
    "Hypothetical with identifier only": "#E69F00",
    "Hypothetical without additional annotation": "#D55E00",
    "No product annotation": "#999999",
}
PLOTLY_CONFIG = {"displaylogo": False, "responsive": True}
PREVIEW_OPTIONS = [500, 1000, 5000, 10000]
MAX_PLOT_CATEGORIES = 300
MAX_LINE_GROUPS = 8


# ============================================================
# CARREGAMENTO (cache compartilhado, sem cópia por rerun)
# ============================================================

def data_version() -> str:
    """Muda quando dados/crosswalk mudam -> invalida o cache automaticamente."""
    parts = []
    for stem in ("records_dashboard", "cds_core", "cds_details", "cds_dashboard"):
        path = L.find_table(DATA_DIR, stem)
        if path is None:
            parts.append(f"{stem}:ausente")
        else:
            stat = path.stat()
            parts.append(f"{path.name}:{stat.st_mtime_ns}:{stat.st_size}")
    for path in (HOST_CROSSWALK_FILE, HOST_NORMALIZATION_AUDIT_FILE, BUILD_AUDIT_FILE, MANIFEST_FILE, INTEREST_FILE, QUALITY_FILE,
                 LOCAL_VALIDATION_FILE, NCBI_VALIDATION_FILE, COUNTRY_EXCEPTIONS_FILE,
                 COUNTRY_NORMALIZATION_AUDIT_FILE, COUNTRY_CROSSWALK_FILE, GENE_CONTEXT_FILE,
                 GLOBAL_ANNOTATION_FILE, GLOBAL_GENE_FILE, GLOBAL_PRODUCT_FILE, GLOBAL_CLUE_FILE, GLOBAL_CLUE_SOURCE_FILE, GLOBAL_TAXON_FILE):
        if path.exists():
            stat = path.stat()
            parts.append(f"{path.name}:{stat.st_mtime_ns}:{stat.st_size}")
    return "|".join(parts)


@st.cache_resource(show_spinner="Carregando registros...")
def get_records(version: str) -> pd.DataFrame:
    return L.load_records(DATA_DIR)


@st.cache_resource(show_spinner="Carregando CDS...")
def get_cds(version: str) -> pd.DataFrame:
    return L.load_cds(DATA_DIR)


@st.cache_resource(show_spinner="Carregando detalhes de CDS...")
def get_cds_details(version: str) -> pd.DataFrame:
    return L.load_cds_details(DATA_DIR)


@st.cache_resource
def get_manifest(path_str: str, mtime_ns: int) -> dict:
    path = Path(path_str)
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


@st.cache_resource
def get_dpol(mtime: float) -> pd.DataFrame:
    return pd.read_csv(DPOL_FILE, sep="\t")


@st.cache_resource
def get_optional_table(path_str: str, mtime_ns: int) -> pd.DataFrame:
    path = Path(path_str)
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


@st.cache_resource
def get_build_audit(path_str: str, mtime_ns: int) -> dict:
    path = Path(path_str)
    if not path.exists():
        return {}
    meta = {}
    sha = {}
    in_sha = False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line == "OUTPUT SHA256":
            in_sha = True
            continue
        if in_sha and "\t" in raw:
            name, digest = raw.split("\t", 1)
            sha[name.strip()] = digest.strip()
            continue
        if "=" in line and not line.startswith("="):
            key, value = line.split("=", 1)
            meta[key.strip()] = value.strip()
    meta["sha256"] = sha
    return meta


@st.cache_resource
def get_filter_options(version: str, _records: pd.DataFrame) -> dict:
    def uniques(col):
        return sorted(known(_records, col)[col].astype(str).unique().tolist())

    years = _records["collection_year_num"].dropna()
    present = set(_records["interest_group"].dropna().astype(str))
    return {
        "countries": uniques("country_canonical"),
        "species": uniques("species_final"),
        "interest": [g for g in INTEREST_ORDER if g in present],
        "years": (int(years.min()), int(years.max())) if len(years) else None,
    }


VERSION = data_version()
try:
    if not DPOL_FILE.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {DPOL_FILE}")
    records = get_records(VERSION)
    dpol_flow = get_dpol(DPOL_FILE.stat().st_mtime_ns)
    interest_summary_global = (get_optional_table(str(INTEREST_FILE), INTEREST_FILE.stat().st_mtime_ns)
                               if INTEREST_FILE.exists() else pd.DataFrame())
    quality_global = (get_optional_table(str(QUALITY_FILE), QUALITY_FILE.stat().st_mtime_ns)
                      if QUALITY_FILE.exists() else pd.DataFrame())
    build_meta = (get_build_audit(str(BUILD_AUDIT_FILE), BUILD_AUDIT_FILE.stat().st_mtime_ns)
                  if BUILD_AUDIT_FILE.exists() else {})
    local_validation = (get_optional_table(str(LOCAL_VALIDATION_FILE), LOCAL_VALIDATION_FILE.stat().st_mtime_ns)
                        if LOCAL_VALIDATION_FILE.exists() else pd.DataFrame())
    ncbi_validation = (get_optional_table(str(NCBI_VALIDATION_FILE), NCBI_VALIDATION_FILE.stat().st_mtime_ns)
                       if NCBI_VALIDATION_FILE.exists() else pd.DataFrame())
    country_exceptions = (get_optional_table(str(COUNTRY_EXCEPTIONS_FILE), COUNTRY_EXCEPTIONS_FILE.stat().st_mtime_ns)
                          if COUNTRY_EXCEPTIONS_FILE.exists() else pd.DataFrame())
    host_normalization_audit = (get_optional_table(str(HOST_NORMALIZATION_AUDIT_FILE), HOST_NORMALIZATION_AUDIT_FILE.stat().st_mtime_ns)
                                if HOST_NORMALIZATION_AUDIT_FILE.exists() else pd.DataFrame())
    country_normalization_audit = (get_optional_table(str(COUNTRY_NORMALIZATION_AUDIT_FILE), COUNTRY_NORMALIZATION_AUDIT_FILE.stat().st_mtime_ns)
                                   if COUNTRY_NORMALIZATION_AUDIT_FILE.exists() else pd.DataFrame())
    country_crosswalk = (get_optional_table(str(COUNTRY_CROSSWALK_FILE), COUNTRY_CROSSWALK_FILE.stat().st_mtime_ns)
                         if COUNTRY_CROSSWALK_FILE.exists() else pd.DataFrame())
    gene_context = (get_optional_table(str(GENE_CONTEXT_FILE), GENE_CONTEXT_FILE.stat().st_mtime_ns)
                    if GENE_CONTEXT_FILE.exists() else pd.DataFrame())
    dataset_manifest = (get_manifest(str(MANIFEST_FILE), MANIFEST_FILE.stat().st_mtime_ns)
                        if MANIFEST_FILE.exists() else {})
    global_annotation = (get_optional_table(str(GLOBAL_ANNOTATION_FILE), GLOBAL_ANNOTATION_FILE.stat().st_mtime_ns)
                         if GLOBAL_ANNOTATION_FILE.exists() else pd.DataFrame())
    global_gene_counts = (get_optional_table(str(GLOBAL_GENE_FILE), GLOBAL_GENE_FILE.stat().st_mtime_ns)
                          if GLOBAL_GENE_FILE.exists() else pd.DataFrame())
    global_product_counts = (get_optional_table(str(GLOBAL_PRODUCT_FILE), GLOBAL_PRODUCT_FILE.stat().st_mtime_ns)
                             if GLOBAL_PRODUCT_FILE.exists() else pd.DataFrame())
    global_clue_counts = (get_optional_table(str(GLOBAL_CLUE_FILE), GLOBAL_CLUE_FILE.stat().st_mtime_ns)
                          if GLOBAL_CLUE_FILE.exists() else pd.DataFrame())
    global_clue_source_counts = (get_optional_table(str(GLOBAL_CLUE_SOURCE_FILE), GLOBAL_CLUE_SOURCE_FILE.stat().st_mtime_ns)
                                 if GLOBAL_CLUE_SOURCE_FILE.exists() else pd.DataFrame())
    global_taxon_summary = (get_optional_table(str(GLOBAL_TAXON_FILE), GLOBAL_TAXON_FILE.stat().st_mtime_ns)
                            if GLOBAL_TAXON_FILE.exists() else pd.DataFrame())
except (FileNotFoundError, L.SchemaError) as error:
    st.error(f"{error}\n\nExecute primeiro `01_build_dashboard_data.py`.")
    st.stop()


# ============================================================
# BARRA LATERAL: FILTROS GLOBAIS
# ============================================================

options = get_filter_options(VERSION, records)
YEAR_BOUNDS = options["years"]


def init_filters_from_query():
    if st.session_state.get("_filters_initialized"):
        return
    raw = st.query_params.get("filters")
    if raw:
        try:
            state = json.loads(str(raw))
            st.session_state["f_interest"] = [x for x in state.get("interest", []) if x in options["interest"]]
            st.session_state["f_countries"] = [x for x in state.get("countries", []) if x in options["countries"]]
            st.session_state["f_species"] = [x for x in state.get("species", []) if x in options["species"]]
            st.session_state["f_undated"] = bool(state.get("undated", True))
            if YEAR_BOUNDS and isinstance(state.get("year"), list) and len(state["year"]) == 2:
                lo = max(YEAR_BOUNDS[0], int(state["year"][0]))
                hi = min(YEAR_BOUNDS[1], int(state["year"][1]))
                if lo <= hi:
                    st.session_state["f_year"] = (lo, hi)
        except Exception:
            pass
    st.session_state["_filters_initialized"] = True


def reset_filters():
    st.session_state["f_interest"] = []
    st.session_state["f_countries"] = []
    st.session_state["f_species"] = []
    st.session_state["f_undated"] = True
    if YEAR_BOUNDS:
        st.session_state["f_year"] = YEAR_BOUNDS
    st.query_params.clear()


init_filters_from_query()

st.sidebar.header("Filtros globais")
st.sidebar.caption(
    "Seleção vazia = todos. Os filtros valem para todas as seções, exceto onde "
    "houver aviso explícito (ex.: Fluxo DPOL)."
)
st.sidebar.button("Limpar filtros", on_click=reset_filters, width="stretch")

sel_interest = st.sidebar.multiselect(
    "Espécies / vírus de interesse", options["interest"], key="f_interest"
)
sel_countries = st.sidebar.multiselect("País", options["countries"], key="f_countries")
sel_species = st.sidebar.multiselect("species_final", options["species"], key="f_species")

if YEAR_BOUNDS and YEAR_BOUNDS[0] < YEAR_BOUNDS[1]:
    sel_year = st.sidebar.slider(
        "Ano de coleta",
        min_value=YEAR_BOUNDS[0], max_value=YEAR_BOUNDS[1],
        value=YEAR_BOUNDS, key="f_year",
    )
else:
    sel_year = None

sel_undated = st.sidebar.checkbox("Incluir registros sem ano", value=True, key="f_undated")

filters = Filters(
    countries=tuple(sel_countries),
    species=tuple(sel_species),
    interest=tuple(sel_interest),
    year_range=tuple(sel_year) if sel_year else None,
    include_undated=sel_undated,
)
RMASK = L.records_mask(records, filters)
filtered_records = records if RMASK is None else records[RMASK]
SIGNATURE = filters.signature()

# Mantém o estado dos filtros na URL para links reproduzíveis.
_query_state = {
    "countries": list(sel_countries),
    "species": list(sel_species),
    "interest": list(sel_interest),
    "year": list(sel_year) if sel_year else None,
    "undated": bool(sel_undated),
}
_encoded_state = json.dumps(_query_state, ensure_ascii=False, separators=(",", ":"))
if len(_encoded_state) < 6000:
    if st.query_params.get("filters") != _encoded_state:
        st.query_params["filters"] = _encoded_state
elif "filters" in st.query_params:
    del st.query_params["filters"]
st.sidebar.caption("🔗 Os filtros atuais ficam registrados na URL para compartilhamento/reprodutibilidade.")


def current_cds() -> pd.DataFrame:
    """Carrega somente o núcleo das CDS quando a página ativa realmente precisa dele."""
    all_cds = get_cds(VERSION)
    if RMASK is None:
        return all_cds
    accessions = records["accession_version"].to_numpy()[RMASK]
    return all_cds[all_cds["accession_version"].isin(accessions)]


def attach_cds_details(view: pd.DataFrame, columns=None) -> pd.DataFrame:
    """Anexa campos detalhados (note/location/etc.) somente às linhas em uso."""
    if view.empty or "cds_row_id" not in view.columns:
        return view
    details = get_cds_details(VERSION)
    ids = view["cds_row_id"]
    subset = details[details["cds_row_id"].isin(ids)]
    return L.merge_cds_details(view, subset, columns=columns)


def cds_search_with_details(core: pd.DataFrame, query: str, core_fields, detail_fields=("note",)) -> pd.DataFrame:
    """Busca no núcleo e, sob demanda, também nos textos detalhados."""
    query = (query or "").strip()
    if not query:
        return core
    core_mask = L.search_mask(core, query, list(core_fields))
    if detail_fields and "cds_row_id" in core.columns:
        details = get_cds_details(VERSION)
        details = details[details["cds_row_id"].isin(core["cds_row_id"])]
        detail_mask = L.search_mask(details, query, [c for c in detail_fields if c in details.columns])
        detail_ids = details.loc[detail_mask, "cds_row_id"]
        core_mask |= core["cds_row_id"].isin(detail_ids).to_numpy()
    return core[core_mask]


# ============================================================
# HELPERS DE INTERFACE
# ============================================================

def require_records():
    """Aviso claro (em vez de abas em branco) quando os filtros zeram tudo."""
    if len(filtered_records) == 0:
        st.warning(
            "Nenhum registro atende aos filtros atuais. "
            "Amplie a seleção ou use **Limpar filtros** na barra lateral."
        )
        st.stop()


def scope_note(filtered: bool = True):
    if not filtered:
        st.caption("Dados globais desta etapa; filtros laterais não se aplicam.")


def style_fig(fig, height=None, legend_title=None):
    fig.update_layout(
        template="plotly_white",
        separators=",.",  # decimais com vírgula, milhar com ponto (pt-BR)
        font=dict(size=14),
        hoverlabel=dict(font_size=13),
        margin=dict(l=20, r=20, t=60, b=30),
        legend=dict(title=legend_title, orientation="h", yanchor="bottom",
                    y=1.02, xanchor="left", x=0),
    )
    if height is not None:
        fig.update_layout(height=height)
    fig.update_xaxes(showgrid=True, gridcolor="rgba(120,120,120,0.16)",
                     zeroline=False, automargin=True)
    fig.update_yaxes(showgrid=False, zeroline=False, automargin=True)
    return fig


def show_fig(fig, height=None, legend_title=None):
    st.plotly_chart(style_fig(fig, height, legend_title),
                    width="stretch", config=PLOTLY_CONFIG)


def count_limit(label, key, default="50"):
    choices = ["20", "50", "100", "200", "Todas"]
    choice = st.selectbox(label, choices, index=choices.index(default), key=key)
    return None if choice == "Todas" else int(choice)


def top_n(table: pd.DataFrame, n):
    return table.copy() if n is None else table.head(n).copy()


def ranked_bar(table, category, value, *, key, title, label, color=COLORS[0],
               default_n="50", show_text=False):
    """Ranking sem truncamento silencioso: 'Todas' inclui todas as categorias."""
    if table.empty:
        return pd.DataFrame()
    n = count_limit(label, key, default_n)
    if n is None and len(table) > MAX_PLOT_CATEGORIES:
        # Para listas muito longas, um treemap inclui todas as categorias sem criar
        # um gráfico vertical com dezenas de milhares de pixels.
        plot = table.copy()
        st.caption(f"Todas as {fmt_int(len(plot))} categorias estão incluídas no treemap.")
        fig = px.treemap(
            plot, path=[category], values=value, color=value,
            color_continuous_scale="Cividis", title=title,
        )
        show_fig(fig, 720)
        return plot
    plot = top_n(table, n).sort_values(value)
    st.caption(f"Gráfico: {fmt_int(len(plot))} de {fmt_int(len(table))} categorias.")
    fig = px.bar(plot, x=value, y=category, orientation="h", title=title,
                 color_discrete_sequence=[color],
                 text=value if show_text else None)
    show_fig(fig, categorical_height(len(plot)))
    return plot

def show_gene_context_legend(gene_names):
    """Legenda biológica logo abaixo do gráfico, sem assumir função universal pelo símbolo."""
    if gene_context.empty or not gene_names:
        return
    wanted = [str(x) for x in gene_names]
    raw = gene_context[gene_context["gene"].astype(str).isin(wanted)].copy()
    if raw.empty:
        return
    order = {g: i for i, g in enumerate(wanted)}
    raw["_order"] = raw["gene"].astype(str).map(order)
    raw = raw.sort_values("_order")

    st.markdown("#### Legenda biológica dos genes exibidos")
    curated = raw[raw.get("ncbi_note", "").astype(str).str.strip().ne("")]
    if not curated.empty:
        for _, row in curated.iterrows():
            gene = str(row.get("gene", "")).strip()
            note = str(row.get("ncbi_note", "")).strip()
            refs = str(row.get("reference_url", "")).strip()
            ref_links = []
            for i, url in enumerate([u.strip() for u in refs.split("|") if u.strip()], start=1):
                ref_links.append(f"[NCBI {i}]({url})")
            ref_txt = " · " + " · ".join(ref_links) if ref_links else ""
            st.markdown(f"**{gene}** — {note}{ref_txt}")

    cols = [c for c in [
        "gene", "top_products", "top_herpes", "top_associations",
        "ncbi_note", "reference_label", "reference_url", "legend_source"
    ] if c in raw.columns]
    legend = raw[cols].rename(columns={
        "gene": "Gene",
        "top_products": "Proteínas/produtos mais observados",
        "top_herpes": "Herpesvírus/táxons mais associados",
        "top_associations": "Associações gene × vírus × proteína",
        "ncbi_note": "Interpretação validada / nota",
        "reference_label": "Fonte",
        "reference_url": "Referência NCBI",
        "legend_source": "Base da legenda",
    })
    with st.expander("Ver contexto completo dos genes no dataset", expanded=False):
        st.caption(
            "Símbolos UL/US podem ter funções diferentes entre herpesvírus. As associações "
            "da tabela são as observadas no próprio dataset; as notas externas são identificadas como NCBI."
        )
        st.dataframe(legend, hide_index=True, width="stretch", height=min(680, 42 * len(legend) + 90))


def safe_slider(label, lo, hi, value, key):
    """st.slider levanta exceção se max <= min; aqui vira valor fixo."""
    if hi <= lo:
        st.caption(f"{label}: {hi} (só há {hi} disponível)")
        return hi
    return st.slider(label, min_value=lo, max_value=hi,
                     value=min(max(value, lo), hi), key=key)


def small_download(label, df, stem, key):
    """Para tabelas pequenas (resumos): serializar é barato."""
    payload, ext, mime = L.to_download(df)
    st.download_button(label, payload, file_name=stem + ext, mime=mime, key=key)


def deferred_download(label, make_df, stem, key, extra=""):
    """
    Para tabelas grandes: só serializa quando o usuário pede, e guarda o
    arquivo em session_state enquanto filtros/consulta não mudarem.
    """
    sig = f"{SIGNATURE}|{extra}"
    state_key = f"dl::{key}"
    ready = st.session_state.get(state_key)
    if ready and ready["sig"] == sig:
        st.download_button(
            f"⬇️ {label}", ready["data"], file_name=stem + ready["ext"],
            mime=ready["mime"], key=f"dlb::{key}",
        )
        return
    if st.button(f"Preparar download: {label}", key=f"prep::{key}"):
        with st.spinner("Gerando arquivo..."):
            payload, ext, mime = L.to_download(make_df())
        st.session_state[state_key] = dict(sig=sig, data=payload, ext=ext, mime=mime)
        st.rerun()


def show_table(df, columns, key, options=PREVIEW_OPTIONS, index=1, caption_noun="linhas"):
    cols = available_columns(df, columns)
    n = st.selectbox("Linhas na prévia", options, index=min(index, len(options) - 1), key=key)
    st.dataframe(df[cols].head(n), hide_index=True, width="stretch")
    if len(df) > n:
        st.caption(
            f"A prévia mostra {fmt_int(n)} {caption_noun}; o download contém "
            f"todas as {fmt_int(len(df))}."
        )
    return cols


def search_form(key, label, placeholder):
    """Busca em st.form: só executa ao enviar, não a cada tecla."""
    with st.form(f"form_{key}"):
        query = st.text_input(label, placeholder=placeholder, key=f"q_{key}")
        st.form_submit_button("Buscar")
    return query.strip()


# ============================================================
# CABEÇALHO E MÉTRICAS GLOBAIS
# ============================================================

st.title("🧬 Orthoherpesviridae Dataset Explorer")
st.caption("NCBI/GenBank • metadados • hospedeiros • proteínas • curadoria DPOL")
st.markdown(
    """
<div style="padding:0.75rem 1rem;border:1px solid rgba(86,180,233,0.32);border-radius:12px;
            background:rgba(0,114,178,0.06);margin:0.25rem 0 0.85rem 0;line-height:1.45;">
<strong>Idealizado por Sara Cândida Ferreira dos Santos</strong><br>
Mestranda em Bioinformática - UFMG<br>
Bióloga - Universidade Federal de Minas Gerais<br>
Técnica em Biotecnologia - COLTEC - UFMG<br>
Currículo Lattes: <a href="http://lattes.cnpq.br/5384262726323885" target="_blank">http://lattes.cnpq.br/5384262726323885</a><br>
LinkedIn: <a href="https://www.linkedin.com/in/sara-cândida-santos/" target="_blank">https://www.linkedin.com/in/sara-cândida-santos/</a><br>
Email: <a href="mailto:saracandida@ufmg.br">saracandida@ufmg.br</a><br>
<small>Dashboard desenvolvido com dataset obtido do NCBI no dia 21/09/2026.</small>
</div>
""",
    unsafe_allow_html=True,
)

PAGE_LABELS = [
    "🏠 Visão geral", "⭐ Espécies de interesse", "🔎 Explorador taxonômico",
    "🌍 Geografia", "📅 Tempo", "🇧🇷 Brasil", "🐒 Hospedeiros",
    "🧫 Genes / CDS", "❓ Proteínas hipotéticas", "🧬 Comparação de proteínas",
    "✅ Qualidade", "🔐 Proveniência / NCBI", "🧪 DPOL / Curadoria", "📋 Registros",
]
m1, m2, m3, m4 = st.columns(4)
m1.metric("Registros no dataset", fmt_int(len(records)))
m2.metric("Registros após filtros", fmt_int(len(filtered_records)))
_total_cds = int(dataset_manifest.get("counts", {}).get("cds", build_meta.get("cds", 0)) or 0)
m3.metric("CDS no dataset", fmt_int(_total_cds) if _total_cds else "—")
m4.metric("Espécies resolvidas no recorte", fmt_int(known(filtered_records, "species_final")["species_final"].nunique()))

# Barra superior semelhante ao modelo anterior, agora com execução preguiçosa.
TOP_TABS = st.tabs(PAGE_LABELS, key="top_tabs", on_change="rerun")


# ============================================================
# SEÇÕES
# ============================================================

def page_overview():
    st.header("Visão geral")
    require_records()
    scope_note()
    recs = filtered_records

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("species_final distintos (resolvidos)",
              fmt_int(known(recs, "species_final")["species_final"].nunique()))
    c2.metric("Países conhecidos",
              fmt_int(known(recs, "country_canonical")["country_canonical"].nunique()))
    c3.metric("Registros brasileiros", fmt_int((recs["country_canonical"] == "Brazil").sum()))
    c4.metric("Hospedeiros conhecidos", fmt_int(known(recs, "host_final")["host_final"].nunique()))

    st.subheader("Espécies mais representadas")
    species_all = value_counts_table(recs["species_final"], "Espécie", "Registros")
    n_unresolved = int(is_missing(recs["species_final"]).sum())
    if n_unresolved:
        st.caption(f"{fmt_int(n_unresolved)} registros sem espécie resolvida "
                   "(Unresolved_species) não entram neste ranking.")
    ranked_bar(species_all, "Espécie", "Registros", key="overview_species_n",
               title="Representação por species_final",
               label="Quantidade de espécies no gráfico")
    with st.expander("Ver tabela completa de espécies"):
        st.dataframe(species_all, hide_index=True, width="stretch")
        small_download("Baixar contagem completa de espécies", species_all,
                       "species_counts", "dl_species_counts")

    st.subheader("Espécies / vírus destacados")
    highlighted = recs[recs["interest_group"] != INTEREST_OTHER]
    interest_counts = value_counts_table(
        highlighted["interest_group"], "Grupo", "Registros"
    ).sort_values("Registros")
    if len(interest_counts):
        fig = px.bar(interest_counts, x="Registros", y="Grupo", orientation="h",
                     color="Grupo", color_discrete_map=INTEREST_COLORS,
                     text="Registros", title="Grupos de interesse")
        show_fig(fig, categorical_height(len(interest_counts)))


def page_interest():
    st.header("Espécies e vírus de interesse")
    require_records()
    scope_note()
    st.caption("Categoria apenas de visualização; a taxonomia normalizada permanece em species_final.")

    focus = filtered_records[filtered_records["interest_group"] != INTEREST_OTHER]
    if focus.empty:
        st.info("Nenhum registro dos grupos de interesse com os filtros atuais.")
        return

    filtered_counts = value_counts_table(focus["interest_group"], "Grupo", "Registros_filtrados")
    if not interest_summary_global.empty and {"interest_group", "records"}.issubset(interest_summary_global.columns):
        global_counts = interest_summary_global[["interest_group", "records"]].rename(
            columns={"interest_group": "Grupo", "records": "Registros_globais"})
        interest_table = filtered_counts.merge(global_counts, on="Grupo", how="outer").fillna(0)
        interest_table["% do global"] = np.where(
            interest_table["Registros_globais"] > 0,
            100 * interest_table["Registros_filtrados"] / interest_table["Registros_globais"], 0.0)
    else:
        interest_table = filtered_counts.assign(Registros_globais=filtered_counts["Registros_filtrados"], **{"% do global": 100.0})

    count_map = dict(zip(interest_table["Grupo"], interest_table["Registros_filtrados"]))
    metric_cols = st.columns(4)
    for i, label in enumerate(INTEREST_ORDER):
        metric_cols[i % 4].metric(label, fmt_int(count_map.get(label, 0)))

    st.subheader("Representação no recorte")
    scale_mode = st.radio(
        "Escala do eixo de registros",
        ["Logarítmica", "Linear"], horizontal=True, key="interest_scale",
        help="A escala logarítmica torna visíveis grupos raros quando HSV-1 domina a contagem.",
    )
    plot_interest = interest_table[interest_table["Registros_filtrados"] > 0].sort_values("Registros_filtrados")
    fig = px.bar(
        plot_interest, x="Registros_filtrados", y="Grupo", orientation="h",
        color="Grupo", color_discrete_map=INTEREST_COLORS, text="Registros_filtrados",
        hover_data={"Registros_globais": True, "% do global": ":.2f"},
        labels={"Registros_filtrados": "Registros filtrados"}, title="Grupos de interesse",
        log_x=(scale_mode == "Logarítmica"),
    )
    show_fig(fig, categorical_height(len(plot_interest)), "Grupo")
    with st.expander("Ver resumo filtrado × global"):
        st.dataframe(interest_table.sort_values("Registros_filtrados", ascending=False), hide_index=True, width="stretch")

    selected = st.multiselect(
        "Comparar", INTEREST_ORDER,
        default=[g for g in INTEREST_ORDER if count_map.get(g, 0) > 0],
        max_selections=8, key="interest_compare",
    )
    compare = focus[focus["interest_group"].isin(selected)]
    if compare.empty:
        return

    st.subheader("Distribuição geográfica")
    by_country = group_count(known(compare, "country_canonical"), ["country_canonical", "interest_group"])
    if len(by_country):
        totals = by_country.groupby("country_canonical")["Registros"].sum().sort_values(ascending=False)
        n = count_limit("Países no gráfico", "interest_country_n")
        cap = len(totals) if n is None else n
        keep = totals.head(cap).index
        fig = px.bar(
            by_country[by_country["country_canonical"].isin(keep)],
            x="country_canonical", y="Registros", color="interest_group", barmode="group",
            color_discrete_map=INTEREST_COLORS,
            labels={"country_canonical": "País", "interest_group": "Espécie / vírus"},
            title="País × grupo de interesse",
        )
        fig.update_xaxes(tickangle=-45)
        show_fig(fig, 560, "Espécie / vírus")

    st.subheader("Distribuição temporal")
    by_year = group_count(compare[compare["collection_year_num"].notna()],
                          ["collection_year_num", "interest_group"])
    if len(by_year):
        fig = px.line(
            by_year, x="collection_year_num", y="Registros", color="interest_group",
            line_dash="interest_group", symbol="interest_group", markers=True,
            color_discrete_map=INTEREST_COLORS,
            labels={"collection_year_num": "Ano", "interest_group": "Espécie / vírus"},
            title="Registros por ano",
        )
        fig.update_xaxes(tickformat="d")
        show_fig(fig, 520, "Espécie / vírus")

    st.subheader("Hospedeiros dos grupos selecionados")
    host_col = "host_final"
    by_host = group_count(known(compare, host_col), [host_col, "interest_group"])
    if len(by_host):
        totals = by_host.groupby(host_col)["Registros"].sum().sort_values(ascending=False)
        n = count_limit("Hospedeiros no gráfico", "interest_hosts_n")
        if n is None and len(totals) > MAX_PLOT_CATEGORIES:
            by_host_plot = by_host.copy()
            fig = px.sunburst(
                by_host_plot, path=[host_col, "interest_group"], values="Registros",
                color="Registros", color_continuous_scale="Cividis",
                title="Todos os hospedeiros × grupo de interesse",
            )
            show_fig(fig, 720)
        else:
            cap = len(totals) if n is None else n
            keep = totals.head(cap).index
            by_host_plot = by_host[by_host[host_col].isin(keep)]
            fig = px.bar(
                by_host_plot, x="Registros", y=host_col, color="interest_group",
                orientation="h", barmode="stack", color_discrete_map=INTEREST_COLORS,
                labels={host_col: "Hospedeiro", "interest_group": "Espécie / vírus"},
                title="Hospedeiro × grupo de interesse",
            )
            show_fig(fig, categorical_height(by_host_plot[host_col].nunique()), "Espécie / vírus")

    st.subheader("Nomes originais no GenBank")
    original = (group_count(compare, ["interest_group", "organism"])
                .sort_values(["interest_group", "Registros"], ascending=[True, False]))
    st.dataframe(original, hide_index=True, width="stretch")

def page_taxon_explorer():
    st.header("Explorador taxonômico")
    require_records()

    mode = st.pills("Campo taxonômico", ["species_final", "organism"], default="species_final", key="explorer_mode")
    base = known(filtered_records, mode)
    taxa_options = sorted(base[mode].astype(str).unique().tolist())
    use_all = st.checkbox("Usar todos os táxons do recorte", value=False, key="explorer_all")
    if use_all:
        explorer = base
        selected = taxa_options
    else:
        selected = st.multiselect(
            "Selecione um ou mais táxons", taxa_options, key="explorer_taxa",
            help="A opção 'Selecionar tudo' do multiselect pode ser usada sem limite artificial de 20 táxons.",
        )
        if not selected:
            st.info("Selecione táxons ou marque 'Usar todos os táxons do recorte'.")
            return
        explorer = base[base[mode].astype(str).isin(selected)]

    if explorer.empty:
        st.info("Nenhum registro para a seleção atual.")
        return

    e1, e2, e3, e4 = st.columns(4)
    e1.metric("Registros", fmt_int(len(explorer)))
    e2.metric("Táxons", fmt_int(explorer[mode].astype(str).nunique()))
    e3.metric("Países", fmt_int(known(explorer, "country_canonical")["country_canonical"].nunique()))
    e4.metric("Hospedeiros", fmt_int(known(explorer, "host_final")["host_final"].nunique()))

    taxon_summary = value_counts_table(explorer[mode], "Táxon", "Registros")
    if len(taxon_summary):
        st.subheader("Representação dos táxons")
        if use_all or len(taxon_summary) > MAX_PLOT_CATEGORIES:
            fig = px.treemap(
                taxon_summary, path=["Táxon"], values="Registros", color="Registros",
                color_continuous_scale="Cividis", title="Todos os táxons do recorte",
            )
            show_fig(fig, 760)
            st.caption(f"O gráfico inclui todos os {fmt_int(len(taxon_summary))} táxons.")
        else:
            fig = px.bar(
                taxon_summary.sort_values("Registros"), x="Registros", y="Táxon",
                orientation="h", color_discrete_sequence=[COLORS[0]],
                title="Táxons selecionados",
            )
            show_fig(fig, categorical_height(len(taxon_summary)))
        with st.expander("Tabela completa de táxons", expanded=use_all):
            st.dataframe(taxon_summary, hide_index=True, width="stretch", height=560)

    ranked_bar(
        value_counts_table(explorer["country_canonical"], "País", "Registros"),
        "País", "Registros", key="explorer_country_n", label="Países no gráfico",
        title="Distribuição geográfica",
    )

    by_year = group_count(explorer[explorer["collection_year_num"].notna()], ["collection_year_num", mode])
    if len(by_year):
        groups = by_year.groupby(mode)["Registros"].sum().sort_values(ascending=False).index.astype(str).tolist()
        if len(groups) > MAX_LINE_GROUPS:
            line_groups = st.multiselect(
                f"Táxons na série temporal (até {MAX_LINE_GROUPS})", groups,
                default=groups[:MAX_LINE_GROUPS], max_selections=MAX_LINE_GROUPS,
                key="explorer_line_taxa",
            )
            by_year_plot = by_year[by_year[mode].astype(str).isin(line_groups)]
        else:
            by_year_plot = by_year
        if len(by_year_plot):
            fig = px.line(
                by_year_plot, x="collection_year_num", y="Registros", color=mode,
                line_dash=mode, symbol=mode, markers=True, color_discrete_sequence=COLORS,
                labels={"collection_year_num": "Ano"}, title="Distribuição temporal",
            )
            fig.update_xaxes(tickformat="d")
            show_fig(fig, 520, mode)

    ranked_bar(
        value_counts_table(explorer["host_final"], "Hospedeiro", "Registros"),
        "Hospedeiro", "Registros", key="explorer_hosts_n", label="Hospedeiros no gráfico",
        title="Hospedeiros", color=COLORS[2],
    )

    cols = show_table(
        explorer,
        ["accession_version", "taxid", "organism", "species_final", "taxonomy_final_status",
         "country_final", "country_canonical", "collection_date", "collection_year_num",
         "host", "host_final", "host_group", "isolate", "strain", "definition"],
        key="explorer_preview_n",
    )
    deferred_download(
        "Baixar registros selecionados", lambda: explorer[cols],
        "taxon_explorer", "dl_explorer", extra=f"{mode}|all={use_all}|{len(selected)}",
    )

def page_geo():
    st.header("Distribuição geográfica")
    require_records()
    scope_note()
    geo = known(filtered_records, "country_canonical")
    if geo.empty:
        st.info("Nenhum registro com país conhecido nos filtros atuais.")
        return

    countries_all = value_counts_table(geo["country_canonical"], "País", "Registros")
    ranked_bar(
        countries_all, "País", "Registros", key="geo_country_n",
        label="Países no gráfico de registros", title="Registros por país (nomes atuais normalizados)",
    )

    mapped, unmapped = L.add_iso3(countries_all, "País", "Registros")
    if len(mapped):
        neon_blue_scale = [
            [0.00, "rgba(10, 92, 155, 0.12)"],
            [0.18, "rgba(0, 153, 219, 0.28)"],
            [0.45, "rgba(0, 210, 255, 0.48)"],
            [0.72, "rgba(0, 145, 255, 0.72)"],
            [1.00, "rgba(0, 82, 255, 0.96)"],
        ]
        fig = go.Figure(go.Choropleth(
            locations=mapped["iso3"], z=mapped["Registros"], text=mapped["País"],
            customdata=mapped[["País", "Registros"]],
            hovertemplate="<b>%{customdata[0]}</b><br>Registros: %{customdata[1]:,.0f}<extra></extra>",
            colorscale=neon_blue_scale,
            colorbar=dict(title="Registros", thickness=14, len=0.72),
            marker=dict(line=dict(color="rgba(78, 226, 255, 0.88)", width=0.9)),
        ))
        fig.update_geos(
            projection_type="natural earth", showframe=False,
            showcoastlines=True, coastlinecolor="rgba(80, 222, 255, 0.85)", coastlinewidth=0.8,
            showcountries=True, countrycolor="rgba(76, 214, 255, 0.58)", countrywidth=0.6,
            showland=True, landcolor="rgba(20, 55, 82, 0.10)",
            showocean=True, oceancolor="rgba(5, 24, 42, 0.03)",
            showlakes=True, lakecolor="rgba(5, 24, 42, 0.04)",
            bgcolor="rgba(0,0,0,0)",
        )
        fig.update_layout(
            title="Mapa mundial de registros", paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)", margin=dict(l=0, r=0, t=55, b=0),
        )
        show_fig(fig, 650)
    st.caption("country_final é preservado; o dashboard usa country_current/country_canonical com nomenclatura atual quando a equivalência é inequívoca.")
    if unmapped:
        lost = countries_all.loc[countries_all["País"].isin(unmapped), "Registros"].sum()
        st.caption(
            f"Fora do mapa por não terem código ISO atual: {', '.join(unmapped)} "
            f"({fmt_int(lost)} registros)."
        )

    st.subheader("Diversidade taxonômica")
    resolved = known(geo, "species_final")
    diversity_all = (
        resolved.groupby("country_canonical", observed=True)["species_final"].nunique()
        .sort_values(ascending=False).rename_axis("País").reset_index(name="Espécies distintas")
    )
    diversity_all = decat(diversity_all)
    ranked_bar(
        diversity_all, "País", "Espécies distintas", key="geo_diversity_n",
        label="Países no gráfico de diversidade", color=COLORS[4],
        title="Número de species_final distintos por país",
    )

    with st.expander("Auditar country_final original × país canônico"):
        audit = group_count(geo, ["country_final", "country_canonical"])
        st.dataframe(audit.sort_values("Registros", ascending=False), hide_index=True, width="stretch")

def page_time():
    st.header("Distribuição temporal")
    require_records()
    scope_note()
    dated = filtered_records[filtered_records["collection_year_num"].notna()]
    if dated.empty:
        st.info("Nenhum registro com ano de coleta nos filtros atuais.")
        return

    # Espécie/país ausentes viram NaN para não inflar a contagem de distintos.
    tmp = pd.DataFrame({
        "Ano": dated["collection_year_num"].to_numpy().astype(int),
        "accession": dated["accession_version"].to_numpy(),
        "species": dated["species_final"].where(~is_missing(dated["species_final"]).to_numpy()),
        "country": dated["country_canonical"].where(~is_missing(dated["country_canonical"]).to_numpy()),
    })
    annual = (
        tmp.groupby("Ano")
        .agg(Registros=("accession", "count"), Espécies=("species", "nunique"),
             Países=("country", "nunique"))
        .reset_index().sort_values("Ano")
    )

    for column, color, title in [
        ("Registros", COLORS[0], "Registros por ano de coleta"),
        ("Espécies", COLORS[4], "Diversidade de species_final por ano"),
    ]:
        fig = px.line(annual, x="Ano", y=column, markers=True,
                      color_discrete_sequence=[color], title=title)
        fig.update_xaxes(tickformat="d")
        show_fig(fig, 500)

    with st.expander("Ver série temporal completa"):
        st.dataframe(annual, hide_index=True, width="stretch")


def page_brazil():
    st.header("Brasil")
    require_records()
    scope_note()
    brazil = filtered_records[filtered_records["country_canonical"] == "Brazil"]
    if brazil.empty:
        st.info("Nenhum registro do Brasil com os filtros atuais (o filtro de país pode estar excluindo o Brasil).")
        return

    b1, b2, b3, b4 = st.columns(4)
    b1.metric("Registros", fmt_int(len(brazil)))
    b2.metric("species_final distintos", fmt_int(known(brazil, "species_final")["species_final"].nunique()))
    b3.metric("Hospedeiros distintos", fmt_int(known(brazil, "host_final")["host_final"].nunique()))
    b4.metric("% dos registros filtrados", fmt_pct(100 * len(brazil) / len(filtered_records), 2))

    ranked_bar(
        value_counts_table(brazil["species_final"], "Espécie", "Registros"),
        "Espécie", "Registros", key="br_species_n", label="Espécies no gráfico",
        title="Espécies mais representadas no Brasil",
    )

    br_interest = value_counts_table(
        brazil[brazil["interest_group"] != INTEREST_OTHER]["interest_group"], "Grupo", "Registros"
    ).sort_values("Registros")
    if len(br_interest):
        st.subheader("Espécies de interesse no Brasil")
        scale_mode = st.radio("Escala", ["Logarítmica", "Linear"], horizontal=True, key="br_interest_scale")
        fig = px.bar(
            br_interest, x="Registros", y="Grupo", orientation="h", color="Grupo",
            color_discrete_map=INTEREST_COLORS, text="Registros", log_x=(scale_mode == "Logarítmica"),
        )
        show_fig(fig, categorical_height(len(br_interest)), "Grupo")

    br_hosts = value_counts_table(brazil["host_final"], "Hospedeiro", "Registros")
    if len(br_hosts):
        st.subheader("Hospedeiros no Brasil")
        ranked_bar(
            br_hosts, "Hospedeiro", "Registros", key="br_hosts_n",
            label="Hospedeiros no gráfico", color=COLORS[2], title="",
        )

def page_hosts():
    st.header("Hospedeiros")
    require_records()
    scope_note()
    crosswalk_loaded = HOST_CROSSWALK_FILE.exists()
    if crosswalk_loaded:
        st.caption("host_final usa a padronização aprovada; host original permanece disponível para auditoria.")
    else:
        st.caption(
            "Nenhum host_crosswalk.tsv foi fornecido. host_final replica o campo host original; nenhuma fusão taxonômica automática foi feita."
        )

    host_col = "host_final"
    host_records = known(filtered_records, host_col)
    h1, h2, h3, h4 = st.columns(4)
    h1.metric("Registros com hospedeiro", fmt_int(len(host_records)))
    h2.metric("Hospedeiros distintos", fmt_int(host_records[host_col].nunique()))
    h3.metric("species_final ligados a hospedeiro",
              fmt_int(known(host_records, "species_final")["species_final"].nunique()))
    h4.metric("Completude", fmt_pct(100 * len(host_records) / len(filtered_records)))
    if host_records.empty:
        st.info("Nenhum registro com hospedeiro nos filtros atuais.")
        return

    host_counts_all = value_counts_table(host_records[host_col], "Hospedeiro", "Registros")

    st.subheader("Ranking de hospedeiros")
    ranked_bar(
        host_counts_all, "Hospedeiro", "Registros", key="hosts_ranking_n",
        label="Quantidade de hospedeiros no ranking", title="Número de registros por hospedeiro",
    )

    st.subheader("Visão proporcional")
    n = count_limit("Hospedeiros no treemap", "hosts_treemap_n")
    n_effective = len(host_counts_all) if n is None else n
    treemap_data = top_n(host_counts_all, n_effective)
    st.caption(f"Treemap: {fmt_int(len(treemap_data))} de {fmt_int(len(host_counts_all))} hospedeiros.")
    fig = px.treemap(
        treemap_data, path=["Hospedeiro"], values="Registros", color="Registros",
        color_continuous_scale="Cividis", title="Proporção relativa dos hospedeiros",
    )
    show_fig(fig, 650)

    groups = known(host_records, "host_group") if "host_group" in host_records.columns else host_records.iloc[0:0]
    if not groups.empty:
        st.subheader("Grupos taxonômicos de hospedeiro (crosswalk curado)")
        group_counts = value_counts_table(groups["host_group"], "Grupo", "Registros")
        ranked_bar(
            group_counts, "Grupo", "Registros", key="host_groups_n",
            label="Grupos no gráfico", title="Registros por grupo de hospedeiro", color=COLORS[2],
        )

    st.subheader("Hospedeiro × táxon")
    dimension = st.radio(
        "Agrupar vírus por", ["species_final", "interest_group"], horizontal=True,
        key="host_matrix_dimension",
    )
    mode = st.radio(
        "Hospedeiros da matriz", ["Mais frequentes", "Escolher manualmente"], horizontal=True,
        key="host_selection_mode",
    )
    hosts_list = host_counts_all["Hospedeiro"].tolist()
    if mode == "Mais frequentes":
        if not hosts_list:
            st.info("Sem hospedeiros para montar a matriz.")
            return
        n_hosts = safe_slider("Número de hospedeiros", 1, min(30, len(hosts_list)), 10, "host_matrix_n")
        chosen = hosts_list[:n_hosts]
    else:
        chosen = st.multiselect(
            "Selecione os hospedeiros", hosts_list, max_selections=30, key="host_matrix_manual"
        )
    if not chosen:
        return

    hm = host_records[host_records[host_col].astype(str).isin(chosen)]
    if dimension == "interest_group":
        hm = hm[hm["interest_group"] != INTEREST_OTHER]
    hm = known(hm, dimension)
    taxa_totals = value_counts_table(hm[dimension], dimension, "n")
    if taxa_totals.empty:
        st.info("Sem táxons para os hospedeiros escolhidos.")
        return
    n_taxa = safe_slider("Táxons na matriz", 1, min(25, len(taxa_totals)), 12, "host_matrix_taxa_n")
    keep = taxa_totals[dimension].head(n_taxa).astype(str)
    hm = hm[hm[dimension].astype(str).isin(keep)]
    matrix = pd.crosstab(hm[host_col].astype(str), hm[dimension].astype(str))

    fig = px.imshow(
        matrix, aspect="auto", color_continuous_scale="Cividis",
        labels={"x": dimension, "y": "Hospedeiro", "color": "Registros"},
        title="Matriz de associação hospedeiro × táxon",
    )
    show_fig(fig, max(500, min(1600, 42 * len(matrix) + 220)))
    st.dataframe(matrix.reset_index(), hide_index=True, width="stretch")

    with st.expander("Auditar nomes originais e normalizados"):
        audit = group_count(host_records, ["host", "host_final", "host_group"])
        st.dataframe(audit.sort_values("Registros", ascending=False), hide_index=True, width="stretch")
        small_download("Baixar auditoria de hospedeiros", audit, "host_normalization_audit", "dl_host_audit")

    with st.expander("Tabela completa de hospedeiros"):
        st.dataframe(host_counts_all, hide_index=True, width="stretch")
        small_download("Baixar contagem completa de hospedeiros", host_counts_all,
                       "host_counts", "dl_host_counts")

def page_genes():
    st.header("Genes e CDS")
    require_records()
    scope_note()

    # Sem filtros, rankings simples vêm de agregados do build e NÃO carregam 447 mil CDS.
    use_global = RMASK is None and not global_gene_counts.empty and not global_product_counts.empty
    cds_f = None
    if use_global:
        counts = dataset_manifest.get("counts", {}) if isinstance(dataset_manifest, dict) else {}
        g1, g2, g3, g4 = st.columns(4)
        g1.metric("CDS no recorte", fmt_int(counts.get("cds", _total_cds)))
        g2.metric("Registros/accessions com CDS", fmt_int(counts.get("cds_unique_accessions", 0)))
        g3.metric("Genes anotados distintos", fmt_int(len(global_gene_counts)))
        g4.metric("Produtos distintos", fmt_int(len(global_product_counts)))
        annotation = global_annotation.copy()
    else:
        cds_f = current_cds()
        if cds_f.empty:
            st.info("Nenhuma CDS com os filtros atuais.")
            return
        g1, g2, g3, g4 = st.columns(4)
        g1.metric("CDS no recorte", fmt_int(len(cds_f)))
        g2.metric("Registros/accessions com CDS", fmt_int(cds_f["accession_version"].nunique()))
        g3.metric("Genes anotados distintos", fmt_int(known(cds_f, "gene")["gene"].nunique()))
        g4.metric("Produtos distintos", fmt_int(known(cds_f, "product")["product"].nunique()))
        annotation = value_counts_table(
            cds_f["protein_annotation_class"], "Classe", "CDS", drop_missing=False
        )

    if len(annotation):
        annotation = annotation.rename(columns={"protein_annotation_class": "Classe"})
        fig = px.bar(
            annotation, x="Classe", y="CDS", color="Classe",
            color_discrete_map=ANNOTATION_CLASS_COLORS, text="CDS",
            title="Classes de anotação proteica",
        )
        fig.update_xaxes(tickangle=-20)
        show_fig(fig, 500, "Classe")

    ranking_metric = st.radio(
        "Base dos rankings de genes/produtos",
        ["Accessions únicos", "CDS"], horizontal=True, key="genes_ranking_metric",
        help="Accessions únicos reduzem o efeito de múltiplas CDS/duplicações no mesmo registro.",
    )
    by_accessions = ranking_metric == "Accessions únicos"
    value_label = "Accessions" if by_accessions else "CDS"

    if use_global:
        metric_col = "accession_count" if by_accessions else "cds_count"
        genes = global_gene_counts[["gene", metric_col]].rename(columns={"gene": "Gene", metric_col: value_label})
        products = global_product_counts[["product", metric_col]].rename(columns={"product": "Produto", metric_col: value_label})
    else:
        genes = L.feature_frequency(cds_f, "gene", by_genomes=by_accessions).rename(columns={"n": value_label, "gene": "Gene"})
        products = L.feature_frequency(cds_f, "product", by_genomes=by_accessions).rename(columns={"n": value_label, "product": "Produto"})

    if len(genes):
        st.subheader("Genes mais representados")
        plotted_genes = ranked_bar(
            genes, "Gene", value_label, key="genes_n", label="Genes no gráfico",
            title=f"Genes por {value_label.lower()}",
        )
        if plotted_genes is not None and not plotted_genes.empty:
            show_gene_context_legend(plotted_genes["Gene"].astype(str).tolist())
    if len(products):
        st.subheader("Produtos mais representados")
        ranked_bar(
            products, "Produto", value_label, key="products_n", label="Produtos no gráfico",
            color=COLORS[2], title=f"Produtos por {value_label.lower()}",
        )

    st.subheader("Busca e tabela de CDS")
    st.caption(
        "Os campos detalhados (por exemplo, note/location) ficam em um arquivo separado e só são carregados quando você executa uma busca ou solicita a prévia."
    )
    with st.form("form_cds"):
        query = st.text_input("Pesquisar", placeholder="polymerase, capsid, tegument, ORF, kinase...", key="q_cds")
        include_details = st.checkbox("Incluir busca em note e exibir campos detalhados", value=True, key="cds_include_details")
        submitted = st.form_submit_button("Buscar / carregar prévia")

    if not submitted:
        st.info("Os gráficos acima já estão disponíveis. Use o formulário para carregar a tabela de CDS sob demanda.")
        return

    if cds_f is None:
        cds_f = current_cds()
    core_fields = [
        "gene", "gene_synonym", "product", "standard_name", "function",
        "protein_id", "locus_tag", "organism", "species_final", "annotation_clue",
    ]
    hits = cds_search_with_details(cds_f, query, core_fields, ("note",) if include_details else ()) if query else cds_f
    if include_details:
        hits = attach_cds_details(hits, columns=["note", "location", "collection_date", "strand"])

    st.write(
        f"CDS encontradas: **{fmt_int(len(hits))}** em **{fmt_int(hits['accession_version'].nunique())}** registros/accessions."
    )
    display_cols = [
        "accession_version", "organism", "species_final", "country_final", "host", "gene",
        "gene_synonym", "product", "standard_name", "function", "note", "protein_id",
        "locus_tag", "annotation_clue", "annotation_clue_source", "protein_annotation_class",
        "location", "collection_date", "strand",
    ]
    cols = show_table(hits, display_cols, key="cds_preview_n", caption_noun="linhas")
    deferred_download(
        "Baixar resultados completos", lambda: hits[cols],
        "CDS_search", "dl_cds", extra=f"{query}|details={include_details}",
    )

def page_hypothetical():
    st.header("Proteínas hipotéticas")
    require_records()
    scope_note()
    st.info(
        "'Hypothetical protein' indica uma CDS prevista sem atribuição funcional suficientemente "
        "estabelecida. O painel separa pistas de anotação funcional de identificadores como "
        "protein_id/locus_tag e mostra também suporte por accession/registro. Nenhuma proteína é "
        "renomeada automaticamente a partir dessas pistas."
    )

    use_global = RMASK is None and not global_taxon_summary.empty and not global_annotation.empty
    cds_f = None
    hyp = None
    if use_global:
        counts = dataset_manifest.get("counts", {}) if isinstance(dataset_manifest, dict) else {}
        total_cds = int(counts.get("cds", _total_cds) or 0)
        hyp_n = int(counts.get("hypothetical_cds", 0) or 0)
        hyp_acc = int(counts.get("hypothetical_accessions", 0) or 0)
        clue_n = int(counts.get("hypothetical_with_annotation_clue", 0) or 0)
        classes = global_annotation.rename(columns={"protein_annotation_class": "Classe"}).copy()
    else:
        cds_f = current_cds()
        if cds_f.empty:
            st.info("Nenhuma CDS com os filtros atuais.")
            return
        classes = value_counts_table(
            cds_f["protein_annotation_class"], "Classe", "CDS", drop_missing=False
        )
        hyp = cds_f[cds_f["is_hypothetical"]]
        total_cds = len(cds_f)
        hyp_n = len(hyp)
        hyp_acc = hyp["accession_version"].nunique()
        clue_n = int((hyp["protein_annotation_class"] == CLUE_CLASS).sum())

    if len(classes):
        fig = px.bar(
            classes, x="Classe", y="CDS", color="Classe",
            color_discrete_map=ANNOTATION_CLASS_COLORS, text="CDS",
            title="Estado da anotação proteica",
        )
        fig.update_xaxes(tickangle=-20)
        show_fig(fig, 520, "Classe")

    h1, h2, h3, h4 = st.columns(4)
    h1.metric("CDS hypothetical", fmt_int(hyp_n))
    h2.metric("Accessions com ≥1 hypothetical", fmt_int(hyp_acc))
    h3.metric("Com pista de anotação", fmt_int(clue_n))
    h4.metric("% das CDS filtradas", fmt_pct(100 * hyp_n / total_cds if total_cds else 0, 2))
    if hyp_n == 0:
        st.info("Nenhuma CDS hipotética com os filtros atuais.")
        return

    st.subheader("Onde estão as proteínas hipotéticas?")
    dimension = st.radio(
        "Agrupar por", ["species_final", "organism", "interest_group"],
        horizontal=True, key="hyp_dimension",
    )
    if use_global:
        summary = global_taxon_summary[global_taxon_summary["dimension"].astype(str).eq(dimension)].copy()
        if "taxon" in summary.columns:
            summary = summary.rename(columns={"taxon": dimension})
    else:
        summary = L.taxon_protein_summary(cds_f, dimension)
    if dimension == "interest_group" and not summary.empty:
        summary = summary[summary[dimension].astype(str) != INTEREST_OTHER]

    default_min = 20
    min_cds = st.number_input(
        "Mínimo de CDS no grupo para calcular/ranquear taxas",
        min_value=1, value=default_min, step=1, key="hyp_min_cds",
        help="O padrão 20 reduz taxas extremas causadas por grupos com pouquíssimas CDS. "
             "Reduza para 1 quando quiser inspecionar táxons raros."
    )
    summary = summary[summary["CDS_total"] >= min_cds] if not summary.empty else summary
    if summary.empty:
        st.warning("Nenhum grupo atinge o mínimo de CDS selecionado. Reduza o limite para explorar grupos raros.")
        return

    rank_options = {
        "Número de CDS hypothetical": ("CDS_hypothetical", "CDS hypothetical"),
        "% das CDS hypothetical": ("Taxa_hypothetical_pct", "% das CDS"),
        "Accessions com ≥1 hypothetical": ("Accessions_com_hypothetical", "Accessions com hypothetical"),
        "% dos accessions com ≥1 hypothetical": ("Taxa_accessions_hypothetical_pct", "% dos accessions"),
    }
    rank_mode = st.radio(
        "Ordenar por", list(rank_options), horizontal=True, key="hyp_rank_mode",
    )
    x_col, x_label = rank_options[rank_mode]
    summary = summary.sort_values(x_col, ascending=False)
    n = count_limit("Táxons no gráfico", "hyp_taxa_n")
    plot = top_n(summary, n).sort_values(x_col)
    if len(plot):
        if n is None and len(plot) > MAX_PLOT_CATEGORIES:
            fig = px.treemap(
                plot, path=[dimension], values=x_col, color="Hypothetical_com_pista",
                color_continuous_scale="Cividis",
                hover_data={"CDS_total": True, "CDS_hypothetical": True, "Accessions_total": True},
                title="Todos os grupos elegíveis — proteínas hipotéticas",
            )
            show_fig(fig, 740)
        else:
            fig = px.bar(
                plot, x=x_col, y=dimension, orientation="h",
                color="Hypothetical_com_pista", color_continuous_scale="Cividis",
                hover_data={
                    "CDS_total": True, "CDS_hypothetical": True, "Hypothetical_com_pista": True,
                    "Accessions_total": True, "Accessions_com_hypothetical": True,
                    "Taxa_hypothetical_pct": ":.2f", "Taxa_accessions_hypothetical_pct": ":.2f",
                },
                labels={x_col: x_label, dimension: "Táxon/grupo"},
                title="Distribuição das proteínas hipotéticas por táxon/grupo",
            )
            fig.update_layout(coloraxis_colorbar=dict(title="Com pista"))
            show_fig(fig, categorical_height(len(plot)))
        st.caption(f"Gráfico: {fmt_int(len(plot))} de {fmt_int(len(summary))} grupos elegíveis.")

    sparse = summary[summary["Accessions_total"] < 3]
    if len(sparse):
        st.warning(
            f"{fmt_int(len(sparse))} grupos exibidos têm menos de 3 accessions. "
            "Taxas nesses grupos devem ser interpretadas com cautela."
        )

    st.subheader("Tipos de pista encontrados")
    if use_global:
        sources = global_clue_source_counts.rename(columns={"annotation_clue_source": "Campo de origem", "cds_count": "CDS"})
        clue_counts = global_clue_counts.rename(columns={"annotation_clue": "Pista de anotação", "cds_count": "CDS"})
    else:
        clue_hyp = hyp[hyp["protein_annotation_class"] == CLUE_CLASS]
        sources = value_counts_table(clue_hyp["annotation_clue_source"], "Campo de origem", "CDS")
        clue_counts = value_counts_table(clue_hyp["annotation_clue"], "Pista de anotação", "CDS")
    if len(sources):
        fig = px.bar(
            sources.sort_values("CDS"), x="CDS", y="Campo de origem", orientation="h",
            color="Campo de origem", color_discrete_sequence=COLORS,
            title="Campo que forneceu a primeira pista adicional",
        )
        show_fig(fig, categorical_height(len(sources), minimum=420), "Campo")
    ranked_bar(
        clue_counts, "Pista de anotação", "CDS", key="hyp_clues_n", label="Pistas no gráfico",
        color=COLORS[2], title="Pistas de anotação mais recorrentes entre hypothetical proteins",
    )

    st.subheader("Explorar CDS hipotéticas")
    st.caption("A tabela detalhada só é carregada quando você solicita esta exploração.")
    class_options = [
        "Hypothetical with annotation clue",
        "Hypothetical with identifier only",
        "Hypothetical without additional annotation",
    ]
    with st.form("form_hyp"):
        class_filter = st.multiselect("Categoria", class_options, key="hyp_class_filter")
        query = st.text_input(
            "Pesquisar nas proteínas hipotéticas",
            placeholder="ORF, tegument, capsid, locus tag, protein ID, espécie...",
            key="q_hyp",
        )
        include_details = st.checkbox("Incluir busca em note e exibir campos detalhados", value=True, key="hyp_include_details")
        submitted = st.form_submit_button("Buscar / carregar prévia")
    if not submitted:
        return

    if cds_f is None:
        cds_f = current_cds()
    view = cds_f[cds_f["is_hypothetical"]]
    if class_filter:
        view = view[view["protein_annotation_class"].astype(str).isin(class_filter)]
    if query:
        view = cds_search_with_details(
            view, query,
            ["organism", "species_final", "gene", "gene_synonym", "annotation_clue",
             "function", "protein_id", "locus_tag"],
            ("note",) if include_details else (),
        )
    if include_details:
        view = attach_cds_details(view, columns=["note", "location", "collection_date", "strand"])

    st.write(
        f"CDS hipotéticas exibidas: **{fmt_int(len(view))}** em "
        f"**{fmt_int(view['accession_version'].nunique())}** registros/accessions."
    )
    display_cols = [
        "accession_version", "organism", "species_final", "interest_group", "country_final",
        "host", "gene", "gene_synonym", "product", "standard_name", "function", "note",
        "protein_id", "locus_tag", "annotation_clue", "annotation_clue_source",
        "protein_annotation_class", "location", "collection_date", "strand",
    ]
    cols = show_table(view, display_cols, key="hyp_preview_n")
    deferred_download(
        "Baixar proteínas hipotéticas filtradas", lambda: view[cols],
        "hypothetical_proteins", "dl_hyp", extra=f"{class_filter}|{query}|details={include_details}",
    )

def page_compare():
    st.header("Comparação de proteínas")
    require_records()
    scope_note()
    st.caption(
        "Compare produtos, genes ou pistas entre táxons. Para interpretação biológica, "
        "as métricas por accession evitam que táxons com mais CDS dominem automaticamente a matriz. Um accession não é assumido como genoma completo."
    )

    cds_f = current_cds()
    if cds_f.empty:
        st.info("Nenhuma CDS com os filtros atuais.")
        return

    dimension = st.radio(
        "Comparar táxons por", ["species_final", "organism", "interest_group"],
        horizontal=True, key="protein_compare_dimension",
    )
    base = known(cds_f, dimension)
    if dimension == "interest_group":
        base = base[base["interest_group"] != INTEREST_OTHER]
    taxa_options = sorted(base[dimension].astype(str).unique().tolist())
    default_taxa = [g for g in INTEREST_ORDER if g in taxa_options] if dimension == "interest_group" else []
    taxa = st.multiselect(
        "Selecione os táxons", taxa_options, default=default_taxa,
        max_selections=12, key="protein_compare_taxa",
    )
    if not taxa:
        st.info("Selecione pelo menos um táxon para iniciar a comparação.")
        return

    selected = base[base[dimension].astype(str).isin(taxa)]
    support = (
        selected.groupby(dimension, observed=True)
        .agg(Accessions=("accession_version", "nunique"), CDS=("accession_version", "size"))
        .reset_index()
    )
    support = decat(support)
    support_display = support.copy()
    support_display["Rótulo"] = support_display.apply(
        lambda r: f"{r[dimension]} (n={fmt_int(r['Accessions'])})", axis=1
    )
    st.dataframe(support_display.sort_values("Accessions", ascending=False), hide_index=True, width="stretch")
    sparse = support[support["Accessions"] < 3]
    if len(sparse):
        st.warning(
            "Táxons com n<3 accessions: "
            + ", ".join(f"{r[dimension]} (n={int(r['Accessions'])})" for _, r in sparse.iterrows())
            + ". Presença/ausência e percentuais devem ser interpretados com cautela."
        )

    feature = st.radio(
        "Comparar características por", ["product", "gene", "annotation_clue"],
        horizontal=True, key="feature_dimension",
    )
    exclude_hyp = st.checkbox(
        "Excluir CDS classificadas como hypothetical da matriz",
        value=(feature == "product"), key="matrix_exclude_hyp",
        help="Útil para impedir que 'hypothetical protein' domine a escala. O perfil de anotação abaixo continua completo.",
    )
    metric_labels = {
        "% dos accessions do táxon": "accession_pct",
        "Accessions com a característica": "accessions",
        "Contagem de CDS": "cds",
        "Presença / ausência": "presence",
    }
    metric_label = st.radio(
        "Métrica da matriz", list(metric_labels), horizontal=True,
        key="matrix_metric", index=0,
    )
    metric = metric_labels[metric_label]

    rank_basis = st.radio(
        "Como escolher as características mais frequentes",
        ["Accessions únicos", "CDS"], horizontal=True, key="feature_rank_basis",
    )
    feature_source = selected.loc[~selected["is_hypothetical"]] if exclude_hyp else selected
    feature_counts = L.feature_frequency(feature_source, feature, by_genomes=(rank_basis == "Accessions únicos"))
    all_features = feature_counts[feature].astype(str).tolist() if len(feature_counts) else []
    if not all_features:
        st.info("Nenhuma característica disponível para a configuração selecionada.")
        return

    mode = st.radio(
        "Escolha das características", ["Mais frequentes", "Selecionar manualmente"],
        horizontal=True, key="feature_select_mode",
    )
    if mode == "Mais frequentes":
        n_feat = safe_slider(
            "Número de características", 1, min(200, len(all_features)), 30,
            "protein_compare_feature_n",
        )
        keep = all_features[:n_feat]
    else:
        keep = st.multiselect(
            "Selecione produtos/genes/pistas", all_features,
            max_selections=100, key="protein_compare_features_manual",
        )
    if not keep:
        st.info("Selecione pelo menos uma característica.")
        return

    matrix, support_matrix = L.taxon_feature_matrix(
        selected, dimension, feature, taxa, features=keep,
        metric=metric, exclude_hypothetical=exclude_hyp,
    )
    if not matrix.empty:
        support_map = support_matrix.set_index(dimension)["Accessions"].to_dict()
        matrix = matrix.rename(columns={c: f"{c} (n={fmt_int(support_map.get(c, 0))})" for c in matrix.columns})
        if metric == "presence":
            scale = [[0.0, "#F2F2F2"], [0.49, "#F2F2F2"], [0.5, "#0072B2"], [1.0, "#0072B2"]]
            label = "Presença"
            range_color = (0, 1)
        elif metric == "accession_pct":
            scale, label, range_color = "Cividis", "% dos accessions", (0, 100)
        elif metric == "accessions":
            scale, label, range_color = "Cividis", "Accessions", None
        else:
            scale, label, range_color = "Cividis", "CDS", None

        fig = px.imshow(
            matrix, aspect="auto", color_continuous_scale=scale,
            labels={"x": f"{dimension} (n=accessions)", "y": feature, "color": label},
            title=f"{feature} × {dimension}", range_color=range_color,
        )
        fig.update_xaxes(side="top")
        show_fig(fig, max(560, min(5000, 32 * len(matrix) + 260)))

        export = matrix.reset_index()
        st.dataframe(export, hide_index=True, width="stretch")
        small_download("Baixar matriz", export, "taxa_protein_matrix", "dl_matrix")

    st.subheader("Perfil de anotação por táxon")
    by_class = group_count(selected, [dimension, "protein_annotation_class"], "CDS")
    if len(by_class):
        support_map = support.set_index(dimension)["Accessions"].to_dict()
        by_class["taxon_label"] = by_class[dimension].astype(str).map(
            lambda x: f"{x} (n={fmt_int(support_map.get(x, 0))})"
        )
        fig = px.bar(
            by_class, x="taxon_label", y="CDS", color="protein_annotation_class",
            barmode="stack", color_discrete_map=ANNOTATION_CLASS_COLORS,
            labels={"protein_annotation_class": "Classe de anotação", "taxon_label": dimension},
            title="Composição das classes de anotação",
        )
        fig.update_xaxes(tickangle=-30)
        show_fig(fig, 560, "Classe de anotação")

        pct = by_class.assign(
            Percentual=100 * by_class["CDS"] / by_class.groupby("taxon_label")["CDS"].transform("sum")
        )
        fig = px.bar(
            pct, x="taxon_label", y="Percentual", color="protein_annotation_class",
            barmode="stack", color_discrete_map=ANNOTATION_CLASS_COLORS,
            labels={"protein_annotation_class": "Classe de anotação", "Percentual": "% das CDS",
                    "taxon_label": dimension},
            title="Composição percentual das classes de anotação",
        )
        fig.update_yaxes(range=[0, 100])
        fig.update_xaxes(tickangle=-30)
        show_fig(fig, 560, "Classe de anotação")

def page_quality():
    st.header("Completude dos metadados")
    require_records()
    scope_note()
    fields = ["country_final", "species_final", "collection_date", "host_final", "isolate", "strain"]
    quality = L.metadata_completeness(filtered_records, fields)
    quality_global_live = L.metadata_completeness(records, fields).rename(columns={
        "available": "available_global", "missing": "missing_global", "available_pct": "available_pct_global"
    })
    comparison = quality.merge(quality_global_live, on="field", how="left")

    fig = px.bar(
        quality, x="field", y="available_pct", text="available_pct",
        color="available_pct", color_continuous_scale="Cividis",
        labels={"field": "Campo", "available_pct": "% disponível"},
        title="Completude dos principais campos no recorte filtrado", range_color=(0, 100),
    )
    fig.update_layout(coloraxis_colorbar=dict(title="% disponível"))
    fig.update_yaxes(range=[0, 100])
    show_fig(fig, 520)
    st.dataframe(comparison, hide_index=True, width="stretch")
    st.caption("Completude recalculada diretamente dos registros com regra única para valores ausentes.")

    st.subheader("Status taxonômico")
    status = (filtered_records["taxonomy_final_status"].astype(str).str.strip()
              .replace({"": "Missing", "nan": "Missing"}))
    tax_status = value_counts_table(status, "Status", "Registros", drop_missing=False)
    fig = px.bar(
        tax_status, x="Registros", y="Status", orientation="h", color="Status",
        color_discrete_sequence=COLORS,
    )
    show_fig(fig, categorical_height(len(tax_status), minimum=420), "Status")
    st.dataframe(tax_status, hide_index=True, width="stretch")

def page_provenance():
    st.header("Proveniência e validação")
    scope_note(filtered=False)
    st.caption("Fontes originais preservadas; transformações derivadas permanecem auditáveis.")

    st.subheader("Versão e política do dataset")
    if dataset_manifest:
        counts = dataset_manifest.get("counts", {})
        policy = dataset_manifest.get("policy", {})
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Versão", dataset_manifest.get("dashboard_data_version", "—"))
        c2.metric("Registros", fmt_int(counts.get("records", 0)))
        c3.metric("CDS", fmt_int(counts.get("cds", 0)))
        c4.metric("Espécies resolvidas", fmt_int(counts.get("species_final_distinct_resolved", 0)))
        st.caption("Arquivos-fonte preservados; normalizações ficam em campos/crosswalks derivados.")
        with st.expander("Ver manifest completo"):
            st.json(dataset_manifest)
    else:
        st.caption("dataset_manifest_v8.json não encontrado; usando apenas BUILD_AUDIT.")

    st.subheader("BUILD_AUDIT")
    if build_meta:
        rows = []
        for key, value in build_meta.items():
            if key == "sha256":
                continue
            rows.append({"Campo": key, "Valor": value})
        if rows:
            st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        sha = build_meta.get("sha256", {})
        if sha:
            st.dataframe(
                pd.DataFrame([{"Arquivo": k, "SHA256": v} for k, v in sha.items()]),
                hide_index=True, width="stretch"
            )
    else:
        st.warning("BUILD_AUDIT.txt não encontrado.")

    st.subheader("Validação local do build")
    if local_validation.empty:
        st.caption(
            "local_validation_summary.tsv ainda não foi gerado. Rode validate_local_data.py "
            "antes do próximo deploy para registrar a auditoria completa."
        )
    else:
        st.dataframe(local_validation, hide_index=True, width="stretch")
        if "status" in local_validation.columns:
            bad = local_validation[~local_validation["status"].astype(str).str.upper().isin(["OK", "INFO"])]
            if len(bad):
                st.warning(f"{fmt_int(len(bad))} verificações locais requerem atenção.")
            else:
                st.success("As verificações locais registradas estão OK/INFO.")


    st.subheader("Padronização de hospedeiros")
    if not host_normalization_audit.empty:
        st.dataframe(host_normalization_audit, hide_index=True, width="stretch")
    if HOST_CROSSWALK_FILE.exists():
        st.caption("Regras aplicadas: correção Homo spaiens→Homo sapiens; corte no primeiro ';'; padronização de capitalização dos dois primeiros termos quando o par é observado no próprio dataset.")

    st.subheader("Normalização atual dos nomes de países")
    if country_normalization_audit.empty:
        st.caption("Auditoria de normalização de países não encontrada.")
    else:
        st.dataframe(country_normalization_audit, hide_index=True, width="stretch")
    if not country_crosswalk.empty:
        changed_countries = country_crosswalk[country_crosswalk.get("changed", "NO").astype(str).eq("YES")] if "changed" in country_crosswalk.columns else country_crosswalk
        if len(changed_countries):
            with st.expander("Ver nomes de países alterados para a nomenclatura atual"):
                st.dataframe(changed_countries, hide_index=True, width="stretch")

    st.subheader("Exceções conservadoras de país para plotagem")
    if country_exceptions.empty:
        st.caption("country_plot_exceptions.tsv não encontrado; a política continua definida no código.")
    else:
        st.dataframe(country_exceptions, hide_index=True, width="stretch")

    st.subheader("Snapshot de validação NCBI")
    if ncbi_validation.empty:
        st.caption(
            "ncbi_validation_summary.tsv ainda não foi gerado. Rode validate_ncbi_current.py. "
            "A ausência desse arquivo NÃO altera os dados locais; apenas significa que não há um "
            "snapshot externo recente anexado ao dashboard."
        )
    else:
        st.dataframe(ncbi_validation, hide_index=True, width="stretch")
        if "status" in ncbi_validation.columns:
            changed = ncbi_validation[~ncbi_validation["status"].astype(str).str.upper().isin(["MATCH", "MATCH_FORMATTING", "PRESENT", "OK", "INFO"])]
            if len(changed):
                st.warning(
                    f"{fmt_int(len(changed))} itens diferem ou requerem revisão em relação ao snapshot NCBI. "
                    "Essas diferenças não são aplicadas automaticamente ao dataset arquivado."
                )
            else:
                st.success("O snapshot NCBI anexado não registra divergências nos itens verificados.")


def page_dpol():
    st.header("Fluxo DPOL / curadoria")
    flow = dpol_flow.copy()
    flow["n"] = pd.to_numeric(flow["n"], errors="coerce")
    flow = flow.dropna(subset=["n"])
    if len(flow):
        funnel_colors = [COLORS[i % len(COLORS)] for i in range(len(flow))]
        fig = go.Figure(go.Funnel(
            y=flow["stage"], x=flow["n"], textinfo="value",
            marker=dict(color=funnel_colors),
            connector=dict(line=dict(color="rgba(0,114,178,0.28)", width=1.2)),
            hovertemplate="<b>%{y}</b><br>n = %{x:,.0f}<extra></extra>",
        ))
        show_fig(fig, 600)
    st.dataframe(dpol_flow, hide_index=True, width="stretch")

def page_records():
    st.header("Tabela geral")
    require_records()
    scope_note()
    query = search_form(
        "records", "Pesquisar nos registros",
        "accession, organismo, espécie, país, hospedeiro...",
    )
    table = filtered_records
    if query:
        table = table[L.search_mask(table, query, [
            "accession_version", "definition", "organism", "species_final", "interest_group",
            "country_final", "country_current", "country_canonical", "host", "host_final", "host_group", "isolate", "strain",
        ])]

    st.write(f"Registros encontrados: **{fmt_int(len(table))}**")
    cols = show_table(
        table,
        ["accession_version", "taxid", "definition", "length_nt", "sequence_status", "organism",
         "species_final", "interest_group", "taxonomy_final_status", "country_final", "country_current", "country_canonical",
         "country_source", "collection_date", "collection_year_num", "geo_loc_name", "host",
         "host_final", "host_group", "lab_host", "isolate", "strain"],
        key="record_preview_n", options=[500, 1000, 5000, 10000], index=1,
    )
    deferred_download(
        "Baixar tabela filtrada completa", lambda: table[cols],
        "Orthoherpes_records_filtered", "dl_records", extra=query,
    )


# ============================================================
# NAVEGAÇÃO SUPERIOR — somente a página escolhida é executada
# ============================================================

PAGE_RENDERERS = {
    "🏠 Visão geral": page_overview,
    "⭐ Espécies de interesse": page_interest,
    "🔎 Explorador taxonômico": page_taxon_explorer,
    "🌍 Geografia": page_geo,
    "📅 Tempo": page_time,
    "🇧🇷 Brasil": page_brazil,
    "🐒 Hospedeiros": page_hosts,
    "🧫 Genes / CDS": page_genes,
    "❓ Proteínas hipotéticas": page_hypothetical,
    "🧬 Comparação de proteínas": page_compare,
    "✅ Qualidade": page_quality,
    "🔐 Proveniência / NCBI": page_provenance,
    "🧪 DPOL / Curadoria": page_dpol,
    "📋 Registros": page_records,
}
for _label, _tab in zip(PAGE_LABELS, TOP_TABS):
    if _tab.open:
        with _tab:
            PAGE_RENDERERS[_label]()
        break


# ============================================================
# RODAPÉ
# ============================================================

st.divider()
files_info = []
for _stem in ("records_dashboard", "cds_core"):
    _path = L.find_table(DATA_DIR, _stem)
    if _path:
        files_info.append(f"{_path.name} ({datetime.fromtimestamp(_path.stat().st_mtime):%d/%m/%Y})")
_sha = build_meta.get("sha256", {}) if isinstance(build_meta, dict) else {}
_sha_short = ", ".join(f"{name}:{digest[:12]}…" for name, digest in _sha.items() if digest)
st.caption(
    f"v{APP_VERSION} • registros NCBI/GenBank ≠ prevalência biológica • "
    f"Dados: {'; '.join(files_info)}" + (f" • SHA256: {_sha_short}" if _sha_short else "")
)
