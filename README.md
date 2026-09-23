# Orthoherpesviridae Dataset Explorer v7

Dashboard Streamlit para exploração do universo inicial de Orthoherpesviridae do NCBI/GenBank,
metadados normalizados, hospedeiros, genes/CDS, proteínas hipotéticas e fluxo de curadoria DPOL.

## Princípio central

O dashboard é uma camada de exploração. Ele **não reescreve silenciosamente os dados arquivados**.
Correções para plotagem, valores ausentes, países e hospedeiros são documentadas como camadas derivadas.

## Arquivos principais

- `streamlit_app.py`: interface;
- `dashboard_lib.py`: lógica testável;
- `DASHBOARD_DATASET_INICIAL/data/records_dashboard.tsv.gz`: registros normalizados;
- `DASHBOARD_DATASET_INICIAL/data/cds_core.tsv.gz`: núcleo analítico das CDS;
- `DASHBOARD_DATASET_INICIAL/data/cds_details.tsv.gz`: detalhes sob demanda;
- `DASHBOARD_DATASET_INICIAL/data/BUILD_AUDIT.txt`: contagens/hashes;
- `DASHBOARD_DATASET_INICIAL/data/dataset_manifest_v7.json`: versão/política;
- `DASHBOARD_DATASET_INICIAL/data/local_validation_summary.tsv`: auditoria local;
- `DASHBOARD_DATASET_INICIAL/data/ncbi_validation_summary.tsv`: snapshot externo de validação.

## Dependências

```text
streamlit==1.59.1
pandas==2.3.1
plotly==7.1.0
numpy==2.4.6
pyarrow==24.0.0
```

## Community Cloud

Use Python 3.11 e `streamlit_app.py` como entrypoint.

## Parquet opcional/recomendado

No servidor com `pyarrow` instalado:

```bash
python convert_to_parquet.py DASHBOARD_DATASET_INICIAL/data
```

O app prefere Parquet quando o manifest comprova que ele corresponde ao build atual; caso contrário faz fallback seguro para TSV.gz.

## Validação local

No pacote completo/auditável:

```bash
python validate_local_data.py DASHBOARD_DATASET_INICIAL/data
```

## NCBI

A validação NCBI é separada do build e nunca corrige automaticamente os arquivos locais:

```bash
python validate_ncbi_current.py DASHBOARD_DATASET_INICIAL/data \
  --email SEU_EMAIL \
  --mode focus
```

## Interpretação

Número de registros no dataset não equivale a prevalência biológica. `accession_version` também não é automaticamente chamado de genoma completo.
