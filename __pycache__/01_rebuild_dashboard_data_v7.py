from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd

VERSION = "7.0.0"

MISSING_TERMS = {
    "", "unknown", "unknown_country", "unknown year", "unknown_year",
    "na", "n/a", "none", "null", "missing", "not provided", "not_provided",
    "not available", "not_available", ".", "-",
}

COUNTRY_PLOT_POLICY = [
    {
        "country_final": "Turkiye",
        "action": "CANONICALIZE_FOR_AGGREGATION_AND_MAP_ONLY",
        "plot_value": "Turkey",
        "reason": "Alias inequívoco; country_final original é preservado.",
    },
    {
        "country_final": "Zaire",
        "action": "MAP_HISTORICAL_GEOGRAPHY_ONLY",
        "plot_value": "Democratic Republic of the Congo",
        "reason": "Nome histórico; usado apenas para posicionamento geográfico, preservando country_final.",
    },
    {
        "country_final": "Korea",
        "action": "DO_NOT_FORCE_MODERN_COUNTRY",
        "plot_value": "",
        "reason": "Ambíguo; não converter automaticamente para South Korea.",
    },
    {
        "country_final": "USSR",
        "action": "DO_NOT_FORCE_MODERN_COUNTRY",
        "plot_value": "",
        "reason": "Entidade histórica; não atribuir a país atual sem evidência adicional.",
    },
    {
        "country_final": "missing",
        "action": "TREAT_AS_MISSING",
        "plot_value": "",
        "reason": "Valor sentinela; não é país disponível.",
    },
]

CORE_CDS_COLUMNS = [
    "accession_version", "taxid", "organism", "host", "gene", "gene_synonym",
    "product", "protein_id", "locus_tag", "standard_name", "function",
    "species_final", "taxonomy_final_status", "country_final", "collection_year_num",
    "interest_group", "annotation_clue", "annotation_clue_source",
    "protein_annotation_class", "is_hypothetical", "gene_normalized",
    "product_normalized",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def human_bytes(n: int) -> str:
    value = float(n)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024:
            return f"{value:.2f} {unit}"
        value /= 1024
    return f"{value:.2f} PB"


def norm(series: pd.Series) -> pd.Series:
    return (
        series.fillna("").astype(str).str.strip().str.casefold()
    )


def missing_mask(series: pd.Series, column: str | None = None) -> pd.Series:
    terms = set(MISSING_TERMS)
    if column == "species_final":
        terms.add("unresolved_species")
    return norm(series).isin(terms) | series.isna()


def available_count(series: pd.Series, column: str | None = None) -> int:
    return int((~missing_mask(series, column)).sum())


def parse_source_audit(path: Path) -> dict:
    out = {"meta": {}, "sha256": {}}
    if not path.exists():
        return out
    in_sha = False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line == "OUTPUT SHA256":
            in_sha = True
            continue
        if in_sha and "\t" in raw:
            name, digest = raw.split("\t", 1)
            out["sha256"][name.strip()] = digest.strip()
            continue
        if "=" in line and not line.startswith("="):
            key, value = line.split("=", 1)
            out["meta"][key.strip()] = value.strip()
    return out


def metadata_quality(records: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for field in ["country_final", "species_final", "collection_date", "host", "isolate", "strain", "geo_loc_name", "definition"]:
        if field not in records.columns:
            continue
        available = available_count(records[field], field)
        rows.append({
            "field": field,
            "available": available,
            "missing": int(len(records) - available),
            "available_pct": round(100 * available / len(records), 4),
        })
    return pd.DataFrame(rows)


def interest_summary(records: pd.DataFrame) -> pd.DataFrame:
    focus = records.loc[records["interest_group"].astype(str).ne("Other")].copy()
    if focus.empty:
        return pd.DataFrame(columns=["interest_group", "records", "species_final", "countries", "hosts"])
    rows = []
    for group, part in focus.groupby("interest_group", dropna=False):
        rows.append({
            "interest_group": group,
            "records": int(len(part)),
            "species_final": int(part.loc[~missing_mask(part["species_final"], "species_final"), "species_final"].nunique()),
            "countries": int(part.loc[~missing_mask(part["country_final"], "country_final"), "country_final"].nunique()),
            "hosts": int(part.loc[~missing_mask(part["host"], "host"), "host"].nunique()),
        })
    return pd.DataFrame(rows).sort_values("records", ascending=False).reset_index(drop=True)


def protein_taxon_summary(cds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    clue = "Hypothetical with annotation clue"
    for dimension in ["species_final", "organism", "interest_group"]:
        if dimension not in cds.columns:
            continue
        valid = cds.loc[~missing_mask(cds[dimension], dimension)].copy()
        if valid.empty:
            continue
        valid["_hyp"] = norm(valid["is_hypothetical"]).eq("yes")
        valid["_clue"] = valid["protein_annotation_class"].astype(str).eq(clue)
        for taxon, part in valid.groupby(dimension, dropna=False):
            acc_total = int(part["accession_version"].nunique())
            cds_total = int(len(part))
            hyp = part.loc[part["_hyp"]]
            rows.append({
                "dimension": dimension,
                "taxon": str(taxon),
                "CDS_total": cds_total,
                "CDS_hypothetical": int(part["_hyp"].sum()),
                "Hypothetical_com_pista": int(part["_clue"].sum()),
                "Accessions_total": acc_total,
                "Accessions_com_hypothetical": int(hyp["accession_version"].nunique()),
                "Taxa_hypothetical_pct": round(100 * part["_hyp"].sum() / cds_total, 6) if cds_total else 0.0,
                "Taxa_accessions_hypothetical_pct": round(100 * hyp["accession_version"].nunique() / acc_total, 6) if acc_total else 0.0,
            })
    return pd.DataFrame(rows)


def feature_counts(cds: pd.DataFrame, column: str, out_name: str) -> pd.DataFrame:
    valid = cds.loc[~missing_mask(cds[column], column), [column, "accession_version"]]
    if valid.empty:
        return pd.DataFrame(columns=[out_name, "cds_count", "accession_count"])
    counts = (
        valid.groupby(column, dropna=False)
        .agg(cds_count=("accession_version", "size"), accession_count=("accession_version", "nunique"))
        .reset_index()
        .rename(columns={column: out_name})
        .sort_values(["accession_count", "cds_count"], ascending=False)
        .reset_index(drop=True)
    )
    return counts


def memory_profile(df: pd.DataFrame, table: str) -> pd.DataFrame:
    mem = df.memory_usage(deep=True)
    rows = []
    for col in df.columns:
        rows.append({
            "table": table,
            "column": col,
            "dtype_source": str(df[col].dtype),
            "nunique": int(df[col].nunique(dropna=False)),
            "memory_mb_object": round(float(mem[col]) / 1024**2, 4),
        })
    return pd.DataFrame(rows).sort_values("memory_mb_object", ascending=False)


def local_validation(records: pd.DataFrame, cds: pd.DataFrame, source_audit: dict, data_dir: Path) -> pd.DataFrame:
    rows = []
    def add(check, status, observed, expected="", note=""):
        rows.append({"check": check, "status": status, "observed": str(observed), "expected": str(expected), "note": note})

    source_names = ["records_dashboard.tsv.gz", "cds_dashboard.tsv.gz"]
    for name in source_names:
        expected = source_audit.get("sha256", {}).get(name, "")
        observed = sha256_file(data_dir / name)
        add(f"sha256:{name}", "OK" if (not expected or observed == expected) else "FAIL", observed, expected,
            "Hash atual versus hash registrado no BUILD_AUDIT de origem.")

    add("records:rows", "OK" if len(records) == 77099 else "REVIEW", len(records), 77099)
    add("records:unique_accession_version", "OK" if records["accession_version"].nunique() == len(records) else "REVIEW",
        records["accession_version"].nunique(), len(records))
    add("cds:rows", "OK" if len(cds) == 446993 else "REVIEW", len(cds), 446993)
    missing_cds_acc = set(cds["accession_version"]) - set(records["accession_version"])
    add("cds:accessions_subset_of_records", "OK" if not missing_cds_acc else "REVIEW", len(missing_cds_acc), 0,
        "Número de accession_version presentes em CDS mas ausentes da tabela de registros.")
    add("cds:unique_accessions", "INFO", cds["accession_version"].nunique(), "", "Accessions com pelo menos uma CDS no inventário.")
    all_species = records["species_final"].nunique(dropna=False)
    resolved_species = records.loc[~missing_mask(records["species_final"], "species_final"), "species_final"].nunique()
    unresolved = int(missing_mask(records["species_final"], "species_final").sum())
    add("taxonomy:species_final_distinct_all", "INFO", all_species, source_audit.get("meta", {}).get("species_final_distinct", ""),
        "Inclui Unresolved_species quando presente.")
    add("taxonomy:species_final_distinct_resolved", "INFO", resolved_species, "", "Exclui valores sem espécie resolvida.")
    add("taxonomy:unresolved_records", "INFO", unresolved, "", "Registros sem species_final resolvido.")
    available_country = available_count(records["country_final"], "country_final")
    old_quality = pd.read_csv(data_dir / "archive_source" / "metadata_quality_source.tsv", sep="\t") if (data_dir / "archive_source" / "metadata_quality_source.tsv").exists() else pd.DataFrame()
    expected_country = ""
    if not old_quality.empty and "field" in old_quality.columns:
        m = old_quality.loc[old_quality["field"].astype(str).eq("country_final"), "available"]
        if len(m): expected_country = int(m.iloc[0])
    add("metadata_quality:country_available_recomputed", "INFO" if expected_country and available_country != expected_country else "OK",
        available_country, expected_country, "Recalcula missing de forma consistente; o literal 'missing' é ausente.")
    hyp = norm(cds["is_hypothetical"]).eq("yes")
    add("protein_annotation:hypothetical_count", "OK" if int(hyp.sum()) == 25086 else "REVIEW", int(hyp.sum()), 25086)
    add("protein_annotation:classes_sum_to_cds", "OK" if cds["protein_annotation_class"].value_counts(dropna=False).sum() == len(cds) else "REVIEW",
        cds["protein_annotation_class"].value_counts(dropna=False).sum(), len(cds))
    for value, note in [
        ("Korea", "Ambíguo para mapa moderno; NÃO mapear automaticamente para South Korea."),
        ("USSR", "Entidade histórica; não atribuir a país atual sem evidência adicional."),
        ("Zaire", "Nome histórico; mapear apenas na camada cartográfica, preservando o original."),
        ("Turkiye", "Alias inequívoco; unificar apenas na camada de plotagem."),
        ("missing", "Valor sentinela; tratar como ausente."),
    ]:
        add(f"country_value:{value}", "INFO", int((records["country_final"] == value).sum()), "", note)
    return pd.DataFrame(rows)


def ncbi_snapshot() -> pd.DataFrame:
    checked_at = "2026-09-23T19:11:26.601564+00:00"
    return pd.DataFrame([
        {
            "validation_type": "taxonomy",
            "key": "focus_taxonomy_names",
            "local_value": "checked=14",
            "ncbi_value": "MATCH=14",
            "status": "MATCH",
            "checked_at_utc": checked_at,
            "note": "TaxIDs/names listed in the focus validation matched NCBI current taxonomy; local archived values were not changed.",
        },
        {
            "validation_type": "nuccore",
            "key": "accession_version:focus",
            "local_value": "checked=176",
            "ncbi_value": "MATCH=176",
            "status": "MATCH",
            "checked_at_utc": checked_at,
            "note": "Current NCBI record queried by accession root; local archived values were not changed.",
        },
        {
            "validation_type": "nuccore",
            "key": "taxid:focus",
            "local_value": "checked=176",
            "ncbi_value": "MATCH=176",
            "status": "MATCH",
            "checked_at_utc": checked_at,
            "note": "Current NCBI TaxID matched for all focus records checked.",
        },
        {
            "validation_type": "nuccore",
            "key": "length_nt:focus",
            "local_value": "checked=176",
            "ncbi_value": "MATCH=176",
            "status": "MATCH",
            "checked_at_utc": checked_at,
            "note": "Current NCBI length matched for all focus records checked.",
        },
        {
            "validation_type": "nuccore",
            "key": "definition:focus",
            "local_value": "v6 flagged 176/176",
            "ncbi_value": "validator formatting artifact identified",
            "status": "INFO",
            "checked_at_utc": checked_at,
            "note": "The v6 comparison was sensitive to terminal punctuation. v6.1 validator distinguishes MATCH_FORMATTING. No local definition was overwritten; rerun v6.1 for a fresh definition status when desired.",
        },
    ])


def main(data_dir: Path) -> None:
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    archive = data_dir / "archive_source"
    archive.mkdir(exist_ok=True)

    source_records = data_dir / "records_dashboard.tsv.gz"
    source_cds = data_dir / "cds_dashboard.tsv.gz"
    source_dpol = data_dir / "dpol_flow.tsv"
    source_audit_path = data_dir / "BUILD_AUDIT_SOURCE.txt"
    for path in (source_records, source_cds, source_dpol):
        if not path.exists():
            raise FileNotFoundError(path)

    # Preserve source summaries if they are present next to this build script.
    script_dir = Path(__file__).resolve().parent
    original_root = script_dir.parent if (script_dir.parent / "metadata_quality.tsv").exists() else script_dir
    for src_name, dst_name in [
        ("metadata_quality.tsv", "metadata_quality_source.tsv"),
        ("interest_summary.tsv", "interest_summary_source.tsv"),
    ]:
        candidate = original_root / src_name
        if candidate.exists() and not (archive / dst_name).exists():
            shutil.copy2(candidate, archive / dst_name)

    source_audit = parse_source_audit(source_audit_path)
    print("Carregando records...")
    records = pd.read_csv(source_records, sep="\t", compression="gzip", dtype=str, low_memory=False).fillna("")
    print("Carregando CDS...")
    cds = pd.read_csv(source_cds, sep="\t", compression="gzip", dtype=str, low_memory=False).fillna("")

    # Integrity first: do not rebuild from an altered source silently.
    for name, path in [(source_records.name, source_records), (source_cds.name, source_cds)]:
        expected = source_audit.get("sha256", {}).get(name)
        observed = sha256_file(path)
        if expected and observed != expected:
            raise RuntimeError(f"SHA256 divergente para {name}: {observed} != {expected}")

    # Split CDS into a light core and details. Every source column is preserved in one of them.
    cds = cds.reset_index(drop=True)
    cds.insert(0, "cds_row_id", np.arange(len(cds), dtype=np.int64))
    core_cols = ["cds_row_id"] + [c for c in CORE_CDS_COLUMNS if c in cds.columns]
    detail_cols = ["cds_row_id"] + [c for c in cds.columns if c not in core_cols and c != "cds_row_id"]
    core = cds[core_cols]
    details = cds[detail_cols]
    core.to_csv(data_dir / "cds_core.tsv.gz", sep="\t", index=False, compression="gzip")
    details.to_csv(data_dir / "cds_details.tsv.gz", sep="\t", index=False, compression="gzip")

    # Correct derived summaries without modifying archived source records/CDS.
    metadata_quality(records).to_csv(data_dir / "metadata_quality.tsv", sep="\t", index=False)
    interest_summary(records).to_csv(data_dir / "interest_summary.tsv", sep="\t", index=False)

    # Conservative country plot policy, with observed counts.
    country_rows = []
    for row in COUNTRY_PLOT_POLICY:
        row = dict(row)
        row["observed_records"] = int((records["country_final"] == row["country_final"]).sum())
        country_rows.append(row)
    pd.DataFrame(country_rows).to_csv(data_dir / "country_plot_exceptions.tsv", sep="\t", index=False)

    # Host candidates for future manual curation only. No host is renamed here.
    host_known = records.loc[~missing_mask(records["host"], "host")].copy()
    host_candidates = (
        host_known.groupby("host", dropna=False)
        .agg(records=("accession_version", "size"),
             species_final_distinct=("species_final", "nunique"),
             interest_group_distinct=("interest_group", "nunique"))
        .reset_index().rename(columns={"host": "host_original"})
        .sort_values("records", ascending=False)
    )
    host_candidates["host_final"] = ""
    host_candidates["host_group"] = ""
    host_candidates["curation_note"] = ""
    host_candidates.to_csv(data_dir / "host_crosswalk_candidates.tsv", sep="\t", index=False)

    # Global pre-aggregates: used when there are no filters, avoiding a 447k-row load for simple charts.
    feature_counts(cds, "gene", "gene").to_csv(data_dir / "cds_global_gene_counts.tsv.gz", sep="\t", index=False, compression="gzip")
    feature_counts(cds, "product", "product").to_csv(data_dir / "cds_global_product_counts.tsv.gz", sep="\t", index=False, compression="gzip")
    feature_counts(cds, "annotation_clue", "annotation_clue").to_csv(data_dir / "cds_global_clue_counts.tsv.gz", sep="\t", index=False, compression="gzip")
    feature_counts(cds, "annotation_clue_source", "annotation_clue_source").to_csv(data_dir / "cds_global_clue_source_counts.tsv", sep="\t", index=False)
    protein_taxon_summary(cds).to_csv(data_dir / "cds_global_taxon_summary.tsv.gz", sep="\t", index=False, compression="gzip")
    cds["protein_annotation_class"].value_counts(dropna=False).rename_axis("protein_annotation_class").reset_index(name="CDS").to_csv(
        data_dir / "cds_global_annotation_class.tsv", sep="\t", index=False)

    # Memory/source profile documents why the online app loads selected columns/categories.
    pd.concat([memory_profile(records, "records_source"), memory_profile(cds.drop(columns=["cds_row_id"]), "cds_source")], ignore_index=True).to_csv(
        data_dir / "DATA_PROFILE.tsv", sep="\t", index=False)

    local = local_validation(records, cds, source_audit, data_dir)
    local.to_csv(data_dir / "local_validation_summary.tsv", sep="\t", index=False)

    snapshot = ncbi_snapshot()
    snapshot.to_csv(data_dir / "ncbi_validation_summary.tsv", sep="\t", index=False)

    # Short human-readable NCBI snapshot; it is a frozen audit, not an automatic correction layer.
    snapshot_md = """# NCBI focus validation snapshot

Validation timestamp from the provided run: `2026-09-23T19:11:26.601564+00:00`.

- 176/176 focus records: accession_version matched.
- 176/176 focus records: TaxID matched.
- 176/176 focus records: length_nt matched.
- Focus taxonomy names shown in the run matched NCBI current taxonomy.
- The v6 definition comparator flagged all 176 because terminal punctuation/formatting was not normalized. The corrected v6.1 validator is bundled; definitions were not overwritten.

This snapshot documents validation only. It does not mutate the archived dataset.
"""
    (data_dir / "NCBI_FOCUS_VALIDATION_SNAPSHOT.md").write_text(snapshot_md, encoding="utf-8")

    generated_files = [
        "records_dashboard.tsv.gz", "cds_dashboard.tsv.gz", "cds_core.tsv.gz", "cds_details.tsv.gz",
        "metadata_quality.tsv", "interest_summary.tsv", "dpol_flow.tsv", "country_plot_exceptions.tsv",
        "host_crosswalk_candidates.tsv", "cds_global_gene_counts.tsv.gz", "cds_global_product_counts.tsv.gz",
        "cds_global_clue_counts.tsv.gz", "cds_global_clue_source_counts.tsv", "cds_global_taxon_summary.tsv.gz", "cds_global_annotation_class.tsv",
        "DATA_PROFILE.tsv", "local_validation_summary.tsv", "ncbi_validation_summary.tsv",
        "NCBI_FOCUS_VALIDATION_SNAPSHOT.md",
    ]
    hashes = {name: sha256_file(data_dir / name) for name in generated_files if (data_dir / name).exists()}
    source_hashes = {
        "records_dashboard.tsv.gz": sha256_file(source_records),
        "cds_dashboard.tsv.gz": sha256_file(source_cds),
    }
    manifest = {
        "dashboard_data_version": VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "policy": {
            "source_records_mutated": False,
            "source_cds_mutated": False,
            "country_canonical_is_plot_layer_only": True,
            "host_normalization_applied": False,
            "ncbi_snapshot_mutates_local_data": False,
        },
        "counts": {
            "records": int(len(records)),
            "unique_accession_version": int(records["accession_version"].nunique()),
            "cds": int(len(cds)),
            "cds_unique_accessions": int(cds["accession_version"].nunique()),
            "species_final_distinct_all": int(records["species_final"].nunique()),
            "species_final_distinct_resolved": int(records.loc[~missing_mask(records["species_final"], "species_final"), "species_final"].nunique()),
            "unresolved_records": int(missing_mask(records["species_final"], "species_final").sum()),
            "hypothetical_cds": int(norm(cds["is_hypothetical"]).eq("yes").sum()),
            "hypothetical_accessions": int(cds.loc[norm(cds["is_hypothetical"]).eq("yes"), "accession_version"].nunique()),
            "hypothetical_with_annotation_clue": int((cds["protein_annotation_class"] == "Hypothetical with annotation clue").sum()),
        },
        "source_sha256": source_hashes,
        "derived_sha256": hashes,
        "cds_split": {"core_columns": core_cols, "detail_columns": detail_cols},
    }
    (data_dir / "dataset_manifest_v7.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    # Audit v7, designed to be human-readable in GitHub/Streamlit.
    lines = [
        f"DASHBOARD BUILD AUDIT V7", "=" * 80,
        f"dashboard_data_version={VERSION}",
        f"records={len(records)}",
        f"unique_accession_version={records['accession_version'].nunique()}",
        f"cds={len(cds)}",
        f"cds_unique_accessions={cds['accession_version'].nunique()}",
        f"species_final_distinct_all={records['species_final'].nunique()}",
        f"species_final_distinct_resolved={records.loc[~missing_mask(records['species_final'], 'species_final'), 'species_final'].nunique()}",
        f"unresolved_records={missing_mask(records['species_final'], 'species_final').sum()}",
        f"country_available_recomputed={available_count(records['country_final'], 'country_final')}",
        f"hypothetical_cds={norm(cds['is_hypothetical']).eq('yes').sum()}",
        "", "SOURCE INTEGRITY", "-" * 80,
    ]
    for name, digest in source_hashes.items():
        expected = source_audit.get("sha256", {}).get(name, "")
        lines.append(f"{name}\t{digest}\t{'MATCH_SOURCE_AUDIT' if expected and digest == expected else 'NO_SOURCE_HASH_OR_REVIEW'}")
    lines += ["", "OUTPUT SHA256", "-" * 80]
    for name, digest in sorted(hashes.items()):
        lines.append(f"{name}\t{digest}")
    lines += ["", "POLICY", "-" * 80,
              "Original records/CDS are preserved byte-for-byte.",
              "Corrections for missing values/country mapping are derived visualization layers only.",
              "Korea and USSR are not forced into modern countries.",
              "No host normalization is applied without a manually curated host_crosswalk.tsv.",
              "NCBI validation is provenance only and never overwrites archived fields."]
    audit_text = "\n".join(lines) + "\n"
    (data_dir / "BUILD_AUDIT_V7.txt").write_text(audit_text, encoding="utf-8")
    (data_dir / "BUILD_AUDIT.txt").write_text(audit_text, encoding="utf-8")

    print("=" * 80)
    print("DASHBOARD DATA V7 REBUILT")
    print("=" * 80)
    print(f"records={len(records):,}; CDS={len(cds):,}")
    print(f"resolved species={manifest['counts']['species_final_distinct_resolved']:,}; unresolved records={manifest['counts']['unresolved_records']:,}")
    print(f"country available (consistent missing rule)={manifest['counts']['records'] - metadata_quality(records).set_index('field').loc['country_final','missing']:,}")
    print(f"hypothetical CDS={manifest['counts']['hypothetical_cds']:,}")
    print("\nCore/details split:")
    print(f"  cds_core: {len(core_cols)} columns -> {human_bytes((data_dir / 'cds_core.tsv.gz').stat().st_size)}")
    print(f"  cds_details: {len(detail_cols)} columns -> {human_bytes((data_dir / 'cds_details.tsv.gz').stat().st_size)}")
    print("\nSource SHA preserved:")
    for name, digest in source_hashes.items():
        print(f"  {name}: {digest}")
    print("\nGenerated:", data_dir)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", nargs="?", default="data")
    args = parser.parse_args()
    main(Path(args.data_dir))
