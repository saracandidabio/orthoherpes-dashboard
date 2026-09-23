"""
Lógica de dados do Orthoherpesviridae Dataset Explorer.

Este módulo NÃO importa Streamlit: tudo aqui é pandas/numpy puro, o que permite
testar com pytest sem subir o app. A interface fica em streamlit_app.py.
"""

from __future__ import annotations

import gzip
import io
import json
import hashlib
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd

# ============================================================
# CONSTANTES
# ============================================================

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
INTEREST_OTHER = "Other"
CLUE_CLASS = "Hypothetical with annotation clue"

# Normalizações conservadoras e auditáveis. Os valores originais continuam
# disponíveis em country_final/host; estas camadas servem apenas para agregação.
COUNTRY_CANONICAL_ALIASES = {
    # Nomes atuais / formas normalizadas para visualização. country_final original
    # permanece intacto no dataset. Casos ambíguos (ex.: Korea, USSR) NÃO são
    # convertidos automaticamente.
    "Turkey": "Türkiye",
    "Turkiye": "Türkiye",
    "Türkiye": "Türkiye",
    "Zaire": "Democratic Republic of the Congo",
    "Democratic Republic of Congo": "Democratic Republic of the Congo",
    "Macedonia": "North Macedonia",
    "Cape Verde": "Cabo Verde",
    "Republic of Korea": "South Korea",
    "Korea, South": "South Korea",
    "Viet Nam": "Vietnam",
    "United States": "USA",
    "United States of America": "USA",
}
HOST_CROSSWALK_FILENAME = "host_crosswalk.tsv"

# Termos que significam "sem informação" (comparados em minúsculas, sem espaços
# nas pontas). "unknown" cobre o "Unknown" de country_final.
MISSING_TERMS = frozenset({
    "", "unknown", "unknown_country", "unknown year", "unknown_year",
    "na", "n/a", "none", "null", "missing", "not provided", "not_provided",
    "not available", "not_available", ".", "-",
})

# Sentinelas específicas de cada coluna, além de MISSING_TERMS.
COLUMN_MISSING_EXTRA = {
    "species_final": frozenset({"unresolved_species"}),
}

# ---- Esquema: coluna -> tipo. Só estas colunas são carregadas. -------------
# "cat"     -> category (colunas repetitivas: maior economia de memória)
# "str"     -> texto livre / alta cardinalidade
# "int"     -> inteiro (downcast)      "float" -> float32
# "bool_yes"-> lido como category e convertido para bool (YES -> True)
#
# Colunas do TSV que o app não usa e que ficam de fora de propósito:
#   CDS: location (~30 MB), country_base, collection_date, collection_year,
#        interest_rule, strand, taxid,
#        identifier_clue / identifier_clue_source (idênticas a protein_id).
RECORDS_SPEC = {
    "accession_version": "str", "taxid": "int", "definition": "str",
    "length_nt": "int", "sequence_status": "cat", "organism": "cat",
    "species_final": "cat", "taxonomy_final_status": "cat",
    "country_final": "cat", "country_source": "cat", "collection_date": "str",
    "collection_year_num": "float", "geo_loc_name": "cat", "host": "cat",
    "lab_host": "cat", "isolate": "str", "strain": "str",
    "interest_group": "cat",
}
CDS_SPEC = {
    "cds_row_id": "int", "accession_version": "cat", "organism": "cat", "host": "cat",
    "gene": "cat", "gene_synonym": "cat", "product": "cat",
    "protein_id": "str", "locus_tag": "cat", "standard_name": "cat",
    "function": "cat", "species_final": "cat",
    "taxonomy_final_status": "cat", "country_final": "cat",
    "collection_year_num": "float", "interest_group": "cat",
    "annotation_clue": "cat", "annotation_clue_source": "cat",
    "protein_annotation_class": "cat", "is_hypothetical": "bool_yes",
    "gene_normalized": "cat", "product_normalized": "cat",
}


# Campos de alta cardinalidade/uso eventual ficam em cds_details.* quando o build v7
# está disponível. Isso reduz RAM nas páginas que só precisam do núcleo analítico.
CDS_DETAILS_SPEC = {
    "cds_row_id": "int",
    "collection_date": "str",
    "country_base": "cat",
    "location": "str",
    "strand": "cat",
    "note": "str",
    "collection_year": "cat",
    "interest_rule": "cat",
    "identifier_clue": "str",
    "identifier_clue_source": "cat",
}

# Colunas sem as quais o app não funciona (o resto é opcional).
REQUIRED_RECORDS = [
    "accession_version", "organism", "species_final", "taxonomy_final_status",
    "country_final", "collection_year_num", "host", "interest_group",
]
REQUIRED_CDS = [
    "accession_version", "organism", "species_final", "country_final",
    "host", "gene", "product", "interest_group", "annotation_clue",
    "annotation_clue_source", "protein_annotation_class", "is_hypothetical",
    "protein_id",
]


class SchemaError(RuntimeError):
    """Arquivo de dados sem as colunas exigidas pelo app."""


# ============================================================
# CARREGAMENTO
# ============================================================

def _build_audit_sha(data_dir: Path) -> dict[str, str]:
    path = Path(data_dir) / "BUILD_AUDIT.txt"
    if not path.exists():
        return {}
    out = {}
    in_sha = False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if raw.strip() == "OUTPUT SHA256":
            in_sha = True
            continue
        if in_sha and "\t" in raw:
            name, digest = raw.split("\t", 1)
            out[name.strip()] = digest.strip()
    return out


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def find_table(data_dir: Path, stem: str) -> Path | None:
    """
    Prefere Parquet. Se existir parquet_manifest.json, valida o SHA da fonte
    contra BUILD_AUDIT.txt para não usar um Parquet gerado de um build antigo.
    Sem manifest, usa mtime como fallback conservador.
    """
    data_dir = Path(data_dir)
    parquet = data_dir / f"{stem}.parquet"
    sources = [p for p in (data_dir / f"{stem}.tsv.gz", data_dir / f"{stem}.tsv") if p.exists()]
    newest_source = max(sources, key=lambda p: p.stat().st_mtime_ns) if sources else None

    if parquet.exists():
        manifest_path = data_dir / "parquet_manifest.json"
        if manifest_path.exists():
            try:
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                info = manifest.get(stem, {})
                expected = info.get("source_sha256", "")
                source_name = info.get("source_name", "")
                audit_sha = _build_audit_sha(data_dir).get(source_name, "")
                if expected and audit_sha:
                    if expected == audit_sha:
                        return parquet
                    return newest_source if newest_source is not None else parquet
            except Exception:
                pass
        if newest_source is None or parquet.stat().st_mtime_ns >= newest_source.stat().st_mtime_ns:
            return parquet
    if newest_source is not None:
        return newest_source
    return parquet if parquet.exists() else None


def _finalize(df: pd.DataFrame, spec: dict[str, str]) -> pd.DataFrame:
    """Aplica os tipos do esquema. Idempotente (serve para TSV e Parquet)."""
    for col, kind in spec.items():
        if col not in df.columns:
            continue
        if kind == "cat":
            if not isinstance(df[col].dtype, pd.CategoricalDtype):
                df[col] = df[col].astype("category")
        elif kind == "int":
            df[col] = pd.to_numeric(df[col], errors="coerce", downcast="integer")
        elif kind == "float":
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("float32")
        elif kind == "bool_yes":
            if df[col].dtype != bool:
                df[col] = df[col].astype("object").eq("YES")
    return df


def read_table(path: Path, spec: dict[str, str]) -> pd.DataFrame:
    wanted = set(spec)
    path = Path(path)
    if path.suffix == ".parquet":
        import pyarrow.parquet as pq  # import tardio: só se houver Parquet

        available = pq.ParquetFile(path).schema.names
        df = pd.read_parquet(path, columns=[c for c in spec if c in available])
    else:
        dtypes = {c: "category" for c, k in spec.items() if k in ("cat", "bool_yes")}
        dtypes.update({c: str for c, k in spec.items() if k == "str"})
        df = pd.read_csv(
            path,
            sep="\t",
            compression="gzip" if path.suffix == ".gz" else None,
            usecols=lambda c: c in wanted,
            dtype=dtypes,
            low_memory=False,
        )
    return _finalize(df, spec)


def missing_columns(df: pd.DataFrame, required: list[str]) -> list[str]:
    return [c for c in required if c not in df.columns]


def validate_schema(df: pd.DataFrame, required: list[str], name: str) -> None:
    missing = missing_columns(df, required)
    if missing:
        raise SchemaError(
            f"A tabela '{name}' não tem as colunas obrigatórias: "
            f"{', '.join(missing)}. Rode novamente 01_build_dashboard_data.py."
        )


def canonicalize_country(series: pd.Series) -> pd.Series:
    """Camada de agregação; nunca apaga o valor original de country_final."""
    obj = series.astype("object").where(series.notna(), "")
    return obj.replace(COUNTRY_CANONICAL_ALIASES)


def load_host_crosswalk(data_dir: Path) -> pd.DataFrame:
    """
    Crosswalk opcional, curado pelo usuário. Colunas aceitas:
    host_original (obrigatória), host_final e host_group (opcionais).
    Linhas sem host_final não alteram o nome original.
    """
    path = Path(data_dir) / HOST_CROSSWALK_FILENAME
    if not path.exists():
        return pd.DataFrame(columns=["host_original", "host_final", "host_group"])
    cw = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "host_original" not in cw.columns:
        raise SchemaError(f"{HOST_CROSSWALK_FILENAME} precisa da coluna host_original.")
    for col in ("host_final", "host_group"):
        if col not in cw.columns:
            cw[col] = ""
    cw = cw[["host_original", "host_final", "host_group"]].drop_duplicates("host_original", keep="last")
    return cw


def apply_host_crosswalk(df: pd.DataFrame, crosswalk: pd.DataFrame) -> pd.DataFrame:
    """Adiciona host_final/host_group sem modificar host."""
    raw = df["host"].astype("object").where(df["host"].notna(), "")
    if crosswalk.empty:
        final = raw
        group = pd.Series("", index=df.index, dtype="object")
    else:
        final_map = dict(zip(crosswalk["host_original"], crosswalk["host_final"]))
        group_map = dict(zip(crosswalk["host_original"], crosswalk["host_group"]))
        mapped = raw.map(final_map).fillna("")
        final = mapped.where(mapped.astype(str).str.strip().ne(""), raw)
        group = raw.map(group_map).fillna("")
    df["host_final"] = pd.Series(final, index=df.index).astype("category")
    df["host_group"] = pd.Series(group, index=df.index).astype("category")
    return df


def load_records(data_dir: Path) -> pd.DataFrame:
    path = find_table(data_dir, "records_dashboard")
    if path is None:
        raise FileNotFoundError(f"records_dashboard.* não encontrado em {data_dir}")
    df = read_table(path, RECORDS_SPEC)
    validate_schema(df, REQUIRED_RECORDS, "records")
    current_country = canonicalize_country(df["country_final"]).astype("category")
    df["country_current"] = current_country
    # Compatibilidade com versões anteriores do dashboard.
    df["country_canonical"] = current_country
    df = apply_host_crosswalk(df, load_host_crosswalk(data_dir))
    return df


def load_cds(data_dir: Path) -> pd.DataFrame:
    """Carrega o núcleo analítico das CDS.

    No build v7, prefere cds_core.* (sem textos/colunas de uso eventual).
    Em builds antigos, faz fallback para cds_dashboard.* preservando compatibilidade.
    """
    path = find_table(data_dir, "cds_core")
    if path is None:
        path = find_table(data_dir, "cds_dashboard")
    if path is None:
        raise FileNotFoundError(f"cds_core.* ou cds_dashboard.* não encontrado em {data_dir}")
    df = read_table(path, CDS_SPEC)
    if "cds_row_id" not in df.columns:
        df.insert(0, "cds_row_id", np.arange(len(df), dtype=np.int64))
    validate_schema(df, REQUIRED_CDS, "CDS")
    df = apply_host_crosswalk(df, load_host_crosswalk(data_dir))
    return df


def load_cds_details(data_dir: Path) -> pd.DataFrame:
    """Carrega campos detalhados somente quando uma página realmente os solicita."""
    path = find_table(data_dir, "cds_details")
    if path is not None:
        df = read_table(path, CDS_DETAILS_SPEC)
        if "cds_row_id" not in df.columns:
            df.insert(0, "cds_row_id", np.arange(len(df), dtype=np.int64))
        return df

    # Compatibilidade com build antigo: extrai apenas as colunas detalhadas do arquivo completo.
    path = find_table(data_dir, "cds_dashboard")
    if path is None:
        return pd.DataFrame(columns=list(CDS_DETAILS_SPEC))
    df = read_table(path, CDS_DETAILS_SPEC)
    if "cds_row_id" not in df.columns:
        df.insert(0, "cds_row_id", np.arange(len(df), dtype=np.int64))
    return df


def merge_cds_details(core: pd.DataFrame, details: pd.DataFrame, columns: list[str] | None = None) -> pd.DataFrame:
    """Anexa detalhes por cds_row_id sem alterar o DataFrame cacheado de origem."""
    if core.empty or details.empty or "cds_row_id" not in core.columns or "cds_row_id" not in details.columns:
        return core.copy()
    if columns is None:
        columns = [c for c in details.columns if c != "cds_row_id"]
    keep = ["cds_row_id"] + [c for c in columns if c in details.columns]
    if len(keep) == 1:
        return core.copy()
    return core.merge(details[keep], on="cds_row_id", how="left", validate="many_to_one")


# ============================================================
# VALORES AUSENTES (fonte única de verdade)
# ============================================================

def _normalize(values: pd.Series) -> pd.Series:
    return values.fillna("").astype(str).str.strip().str.casefold()


def is_missing(series: pd.Series, column: str | None = None) -> pd.Series:
    """
    True onde o valor é nulo ou equivale a "sem informação".

    Vale para QUALQUER coluna (país, espécie, hospedeiro, gene...). Sentinelas
    próprias de cada coluna (ex.: Unresolved_species) vêm de COLUMN_MISSING_EXTRA.
    Em colunas category a normalização roda só nas categorias, não nas linhas.
    """
    column = column or series.name
    terms = MISSING_TERMS | COLUMN_MISSING_EXTRA.get(column, frozenset())

    if isinstance(series.dtype, pd.CategoricalDtype):
        cats = pd.Series(series.cat.categories.astype(str))
        lookup = np.append(_normalize(cats).isin(terms).to_numpy(), True)  # código -1 (NaN) -> ausente
        return pd.Series(lookup[series.cat.codes.to_numpy()], index=series.index)

    return _normalize(series).isin(terms) | series.isna()


def known(df: pd.DataFrame, column: str) -> pd.DataFrame:
    """Linhas em que `column` tem informação de verdade."""
    return df[~is_missing(df[column], column).to_numpy()]


# ============================================================
# BUSCA DE TEXTO
# ============================================================

def contains_mask(series: pd.Series, query: str) -> pd.Series:
    """Substring sem distinção de caixa (regex=False, então sem re.escape)."""
    if isinstance(series.dtype, pd.CategoricalDtype):
        cats = pd.Series(series.cat.categories.astype(str))
        hit = cats.str.contains(query, case=False, regex=False, na=False).to_numpy()
        lookup = np.append(hit, False)  # código -1 (NaN) -> não casa
        return pd.Series(lookup[series.cat.codes.to_numpy()], index=series.index)
    return series.astype(str).str.contains(query, case=False, regex=False, na=False) & series.notna()


def search_mask(df: pd.DataFrame, query: str, columns: list[str]) -> np.ndarray:
    query = (query or "").strip()
    if not query:
        return np.ones(len(df), dtype=bool)
    mask = np.zeros(len(df), dtype=bool)
    for col in columns:
        if col in df.columns:
            mask |= contains_mask(df[col], query).to_numpy()
    return mask


# ============================================================
# FILTROS
# ============================================================

@dataclass(frozen=True)
class Filters:
    countries: tuple = ()
    species: tuple = ()
    interest: tuple = ()
    year_range: tuple | None = None
    include_undated: bool = True

    def signature(self) -> str:
        return repr(self)


def records_mask(records: pd.DataFrame, f: Filters) -> np.ndarray | None:
    """Máscara booleana dos registros, ou None se nenhum filtro está ativo."""
    mask = None

    def combine(current, new):
        new = np.asarray(new, dtype=bool)
        return new if current is None else current & new

    if f.countries:
        country_col = "country_canonical" if "country_canonical" in records.columns else "country_final"
        mask = combine(mask, records[country_col].isin(f.countries))
    if f.species:
        mask = combine(mask, records["species_final"].isin(f.species))
    if f.interest:
        mask = combine(mask, records["interest_group"].isin(f.interest))

    years = records["collection_year_num"]
    full_range = (years.min(), years.max())
    if f.year_range is not None and (
        tuple(f.year_range) != full_range or not f.include_undated
    ):
        in_range = years.between(f.year_range[0], f.year_range[1], inclusive="both")
        if f.include_undated:
            in_range = in_range | years.isna()
        mask = combine(mask, in_range)

    return mask


class FilteredData:
    """
    Aplica os filtros globais UMA vez por execução.

    - Sem filtro ativo, devolve os DataFrames originais (zero cópia): quem
      consome NÃO pode alterá-los (eles vivem no cache compartilhado).
    - O CDS é filtrado pelos accession_version dos registros filtrados (o CDS
      herda os metadados dos registros; ver test_cds_filter_matches_columns) e
      só é materializado quando alguma seção o usa.
    """

    def __init__(self, records: pd.DataFrame, cds: pd.DataFrame, filters: Filters):
        self.all_records = records
        self.all_cds = cds
        self.filters = filters
        self._rmask = records_mask(records, filters)
        self.records = records if self._rmask is None else records[self._rmask]

    @cached_property
    def _cmask(self) -> np.ndarray | None:
        if self._rmask is None:
            return None
        accessions = self.all_records["accession_version"].to_numpy()[self._rmask]
        return self.all_cds["accession_version"].isin(accessions).to_numpy()

    @cached_property
    def cds(self) -> pd.DataFrame:
        return self.all_cds if self._cmask is None else self.all_cds[self._cmask]

    @property
    def n_records(self) -> int:
        return len(self.records)

    @property
    def n_cds(self) -> int:
        return len(self.all_cds) if self._cmask is None else int(self._cmask.sum())


# ============================================================
# AGREGAÇÕES
# ============================================================

def decat(df: pd.DataFrame) -> pd.DataFrame:
    """Converte colunas category em texto, para o Plotly/crosstab não mostrarem
    categorias que não existem mais depois de filtrar."""
    cat_cols = [c for c in df.columns if isinstance(df[c].dtype, pd.CategoricalDtype)]
    return df.astype({c: "object" for c in cat_cols}) if cat_cols else df


def value_counts_table(series: pd.Series, name: str, value_name: str,
                       drop_missing: bool = True) -> pd.DataFrame:
    s = series[~is_missing(series).to_numpy()] if drop_missing else series.dropna()
    counts = s.value_counts()
    counts = counts[counts > 0]  # category devolve categorias vazias com 0
    out = counts.rename_axis(name).reset_index(name=value_name)
    out[name] = out[name].astype(str)
    return out


def group_count(df: pd.DataFrame, keys: list[str], name: str = "Registros") -> pd.DataFrame:
    return decat(df.groupby(keys, observed=True).size().reset_index(name=name))


def taxon_protein_summary(cds: pd.DataFrame, dimension: str) -> pd.DataFrame:
    """Resumo por táxon com contagens de CDS e de accessions únicos."""
    cols = [dimension, "accession_version", "is_hypothetical", "protein_annotation_class"]
    base = cds.loc[~is_missing(cds[dimension]).to_numpy(), cols]
    if base.empty:
        return pd.DataFrame(columns=[
            dimension, "CDS_total", "CDS_hypothetical", "Hypothetical_com_pista",
            "Taxa_hypothetical_pct", "Accessions_total", "Accessions_com_hypothetical",
            "Taxa_accessions_hypothetical_pct",
        ])
    tmp = pd.DataFrame({
        dimension: base[dimension],
        "accession_version": base["accession_version"],
        "CDS_hypothetical": base["is_hypothetical"].astype(int),
        "Hypothetical_com_pista": (base["protein_annotation_class"] == CLUE_CLASS).astype(int),
    })
    cds_summary = (
        tmp.groupby(dimension, observed=True)
        .agg(CDS_total=("CDS_hypothetical", "size"),
             CDS_hypothetical=("CDS_hypothetical", "sum"),
             Hypothetical_com_pista=("Hypothetical_com_pista", "sum"))
        .reset_index()
    )
    genomes_total = (base.groupby(dimension, observed=True)["accession_version"]
                     .nunique().rename("Accessions_total").reset_index())
    genomes_hyp = (base.loc[base["is_hypothetical"]]
                   .groupby(dimension, observed=True)["accession_version"]
                   .nunique().rename("Accessions_com_hypothetical").reset_index())
    summary = cds_summary.merge(genomes_total, on=dimension, how="left").merge(
        genomes_hyp, on=dimension, how="left")
    summary["Accessions_com_hypothetical"] = summary["Accessions_com_hypothetical"].fillna(0).astype(int)
    summary = decat(summary)
    summary["Taxa_hypothetical_pct"] = np.where(
        summary["CDS_total"] > 0, 100 * summary["CDS_hypothetical"] / summary["CDS_total"], 0.0)
    summary["Taxa_accessions_hypothetical_pct"] = np.where(
        summary["Accessions_total"] > 0,
        100 * summary["Accessions_com_hypothetical"] / summary["Accessions_total"], 0.0)
    return summary


def feature_frequency(cds: pd.DataFrame, feature: str, by_genomes: bool = True) -> pd.DataFrame:
    """Ranking de características por accessions únicos (padrão) ou CDS."""
    base = known(cds, feature)
    if base.empty:
        return pd.DataFrame(columns=[feature, "n"])
    if by_genomes:
        out = (base.groupby(feature, observed=True)["accession_version"].nunique()
               .sort_values(ascending=False).rename("n").reset_index())
    else:
        out = value_counts_table(base[feature], feature, "n")
    return decat(out)


def taxon_feature_matrix(
    cds: pd.DataFrame, dimension: str, feature: str, taxa: list[str] | tuple[str, ...],
    features: list[str] | tuple[str, ...] | None = None, metric: str = "accession_pct",
    exclude_hypothetical: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Matriz táxon × característica. metric: cds, accessions, accession_pct, presence.
    O denominador de accession_pct é o nº de accession_version únicos do táxon no recorte.
    Isso NÃO implica que cada accession seja um genoma completo.
    Retorna também a tabela de suporte (n de accessions e CDS por táxon).
    """
    base = cds.loc[~is_missing(cds[dimension]).to_numpy()]
    base = base[base[dimension].astype(str).isin(list(taxa))]
    if base.empty:
        return pd.DataFrame(), pd.DataFrame(columns=[dimension, "Accessions", "CDS"])
    support = (base.groupby(dimension, observed=True)
               .agg(Accessions=("accession_version", "nunique"), CDS=("accession_version", "size"))
               .reset_index())
    support = decat(support)
    source = base
    if exclude_hypothetical and "is_hypothetical" in source.columns:
        source = source.loc[~source["is_hypothetical"]]
    source = source.loc[~is_missing(source[feature]).to_numpy()]
    if features is not None:
        source = source[source[feature].astype(str).isin(list(features))]
    if source.empty:
        return pd.DataFrame(), support

    if metric == "cds":
        grouped = source.groupby([feature, dimension], observed=True).size().rename("value").reset_index()
    else:
        uniq = source[[feature, dimension, "accession_version"]].drop_duplicates()
        grouped = (uniq.groupby([feature, dimension], observed=True).size()
                   .rename("value").reset_index())

    grouped = decat(grouped)
    matrix = grouped.pivot(index=feature, columns=dimension, values="value").fillna(0)
    if features is not None:
        matrix = matrix.reindex([f for f in features if f in matrix.index])
    if metric == "accession_pct":
        den = support.set_index(dimension)["Accessions"]
        for col in matrix.columns:
            matrix[col] = 100 * matrix[col].astype(float) / max(float(den.get(col, 0)), 1.0)
    elif metric == "presence":
        matrix = matrix.gt(0).astype(int)
    return matrix, support


def metadata_completeness(records: pd.DataFrame, fields: list[str]) -> pd.DataFrame:
    rows = []
    total = len(records)
    for field in fields:
        if field not in records.columns:
            continue
        available = int((~is_missing(records[field], field)).sum())
        rows.append({
            "field": field,
            "available": available,
            "missing": total - available,
            "available_pct": round(100 * available / total, 2) if total else 0.0,
        })
    return pd.DataFrame(rows)


# ============================================================
# PAÍSES -> ISO-3 (para o mapa)
# ============================================================
# Nomes exatamente como aparecem em country_final (INSDC). Nomes duplicados
# Aliases geográficos inequívocos podem ser somados no mapa (ex.: Turkey/Turkiye,
# Zaire/RD Congo). Termos historicamente/semanticamente ambíguos, como "Korea",
# NÃO são forçados para um país moderno. Sem entrada de propósito: "USSR" e "Korea".

COUNTRY_ISO3 = {
    "Afghanistan": "AFG", "Algeria": "DZA", "Angola": "AGO", "Argentina": "ARG",
    "Australia": "AUS", "Austria": "AUT", "Bahamas": "BHS", "Bangladesh": "BGD",
    "Barbados": "BRB", "Belarus": "BLR", "Belgium": "BEL", "Bolivia": "BOL",
    "Bosnia and Herzegovina": "BIH", "Botswana": "BWA", "Brazil": "BRA",
    "Bulgaria": "BGR", "Cambodia": "KHM", "Cameroon": "CMR", "Canada": "CAN",
    "Cape Verde": "CPV", "Cabo Verde": "CPV", "Cayman Islands": "CYM", "Central African Republic": "CAF",
    "Chile": "CHL", "China": "CHN", "Colombia": "COL", "Costa Rica": "CRI",
    "Croatia": "HRV", "Cuba": "CUB", "Czechia": "CZE", "Côte d'Ivoire": "CIV",
    "Democratic Republic of the Congo": "COD", "Denmark": "DNK", "Ecuador": "ECU",
    "Egypt": "EGY", "Ethiopia": "ETH", "Finland": "FIN", "France": "FRA",
    "French Guiana": "GUF", "Gabon": "GAB", "Gambia": "GMB", "Georgia": "GEO",
    "Germany": "DEU", "Ghana": "GHA", "Greece": "GRC", "Greenland": "GRL",
    "Grenada": "GRD", "Guatemala": "GTM", "Guinea": "GIN", "Guinea-Bissau": "GNB",
    "Haiti": "HTI", "Hong Kong": "HKG", "Hungary": "HUN", "Iceland": "ISL",
    "India": "IND", "Indonesia": "IDN", "Iran": "IRN", "Iraq": "IRQ",
    "Ireland": "IRL", "Israel": "ISR", "Italy": "ITA", "Jamaica": "JAM",
    "Japan": "JPN", "Kazakhstan": "KAZ", "Kenya": "KEN",
    "Kuwait": "KWT", "Kyrgyzstan": "KGZ", "Laos": "LAO", "Libya": "LBY",
    "Luxembourg": "LUX", "Macedonia": "MKD", "North Macedonia": "MKD", "Madagascar": "MDG", "Malawi": "MWI",
    "Malaysia": "MYS", "Mali": "MLI", "Martinique": "MTQ", "Mauritius": "MUS",
    "Mexico": "MEX", "Mongolia": "MNG", "Morocco": "MAR", "Myanmar": "MMR",
    "Namibia": "NAM", "Netherlands": "NLD", "New Caledonia": "NCL",
    "New Zealand": "NZL", "Nicaragua": "NIC", "Nigeria": "NGA", "Norway": "NOR",
    "Pakistan": "PAK", "Panama": "PAN", "Papua New Guinea": "PNG", "Peru": "PER",
    "Philippines": "PHL", "Poland": "POL", "Portugal": "PRT", "Puerto Rico": "PRI",
    "Republic of the Congo": "COG", "Russia": "RUS", "Rwanda": "RWA",
    "Saint Kitts and Nevis": "KNA", "Sao Tome and Principe": "STP",
    "Saudi Arabia": "SAU", "Senegal": "SEN", "Serbia": "SRB", "Singapore": "SGP",
    "Slovakia": "SVK", "Slovenia": "SVN", "Solomon Islands": "SLB",
    "South Africa": "ZAF", "South Korea": "KOR", "Spain": "ESP", "Sri Lanka": "LKA",
    "Sudan": "SDN", "Sweden": "SWE", "Switzerland": "CHE", "Taiwan": "TWN",
    "Tanzania": "TZA", "Thailand": "THA", "Tunisia": "TUN", "Turkey": "TUR", "Türkiye": "TUR",
    "Turkiye": "TUR", "USA": "USA", "Uganda": "UGA", "Ukraine": "UKR",
    "United Kingdom": "GBR", "Uruguay": "URY", "Vanuatu": "VUT", "Vietnam": "VNM",
    "Wallis and Futuna": "WLF", "Zaire": "COD", "Zambia": "ZMB", "Zimbabwe": "ZWE",
}


def canonical_country_counts(counts: pd.DataFrame, country_col: str, value_col: str) -> pd.DataFrame:
    """Soma aliases conservadores também nos rankings, preservando os dados brutos."""
    work = counts.copy()
    work[country_col] = work[country_col].replace(COUNTRY_CANONICAL_ALIASES)
    return (work.groupby(country_col, as_index=False)[value_col].sum()
            .sort_values(value_col, ascending=False))


def add_iso3(counts: pd.DataFrame, country_col: str, value_col: str
             ) -> tuple[pd.DataFrame, list[str]]:
    """
    Converte uma tabela país -> valor em ISO-3, somando nomes que apontam para
    o mesmo país. Devolve (tabela mapeada, países sem código ISO).
    """
    iso = counts[country_col].map(COUNTRY_ISO3)
    unmapped = sorted(counts.loc[iso.isna(), country_col].astype(str).tolist())
    mapped = counts.assign(iso3=iso).dropna(subset=["iso3"])
    out = (
        mapped.groupby("iso3", as_index=False)
        .agg(**{value_col: (value_col, "sum"),
                country_col: (country_col, lambda s: " / ".join(sorted(s)))})
    )
    return out, unmapped


# ============================================================
# UTILITÁRIOS
# ============================================================

def available_columns(df: pd.DataFrame, columns: list[str]) -> list[str]:
    return [c for c in columns if c in df.columns]


def fmt_int(n) -> str:
    """77099 -> '77.099' (separador de milhar pt-BR)."""
    return f"{int(n):,}".replace(",", ".")


def fmt_pct(x: float, digits: int = 1) -> str:
    return f"{x:.{digits}f}".replace(".", ",") + "%"


def categorical_height(n_items: int, minimum: int = 440, per_item: int = 27,
                       maximum: int = 5000) -> int:
    return max(minimum, min(maximum, int(n_items) * per_item + 160))


def to_download(df: pd.DataFrame, max_plain_rows: int = 100_000) -> tuple[bytes, str, str]:
    """
    Serializa sob demanda. Para tabelas grandes grava direto em gzip/BytesIO,
    evitando manter simultaneamente uma cópia TSV bruta e outra comprimida.
    """
    if len(df) > max_plain_rows:
        buffer = io.BytesIO()
        with gzip.GzipFile(fileobj=buffer, mode="wb", compresslevel=5) as gz:
            with io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:
                df.to_csv(text, sep="\t", index=False)
        return buffer.getvalue(), ".tsv.gz", "application/gzip"
    text = io.StringIO()
    df.to_csv(text, sep="\t", index=False)
    return text.getvalue().encode("utf-8"), ".tsv", "text/tab-separated-values"


def memory_usage_table(df: pd.DataFrame) -> pd.DataFrame:
    """Tabela auditável de memória por coluna, em MB."""
    m = (df.memory_usage(deep=True) / (1024 ** 2)).sort_values(ascending=False)
    return m.rename("MB").rename_axis("coluna").reset_index()
