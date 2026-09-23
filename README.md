# Orthoherpesviridae Dataset Explorer

Interactive Streamlit dashboard for exploration of the
Orthoherpesviridae NCBI/GenBank dataset and curated DPOL workflow.

## Entrypoint

`streamlit_app.py`

## Main datasets

- 77,099 normalized NCBI/GenBank records
- 446,993 CDS
- Taxonomic normalization
- Geographic and temporal metadata
- Hosts
- Genes and protein products
- Hypothetical proteins
- Species of interest
- DPOL curation flow

## Important interpretation

Record representation in the dataset must not be interpreted
as biological prevalence.

`species_final` preserves the normalized taxonomy.

`organism` preserves the original organism designation associated
with the GenBank record.

## Deployment

Deploy with Streamlit Community Cloud using:

- Branch: `main`
- Main file: `streamlit_app.py`
- Python: 3.11

