# Orthoherpesviridae Dataset Explorer — v8

Dashboard Streamlit do universo inicial de Orthoherpesviridae/DPOL.

## Publicação
- Python 3.11
- arquivo principal: `streamlit_app.py`
- dependências fixadas em `requirements.txt`

## Estrutura
- `streamlit_app.py`: interface
- `dashboard_lib.py`: lógica de dados
- `DASHBOARD_DATASET_INICIAL/data/`: dados de deploy
- `BUILD_AUDIT.txt`: contagens e SHA256
- `host_crosswalk.tsv`: padronização aprovada de hospedeiros
- `HOST_NORMALIZATION_AUDIT.tsv`: resumo quantitativo das mudanças

## v8
- navegação superior por abas, com execução preguiçosa;
- explorador taxonômico sem limite de 20 táxons;
- "Todas" não trunca silenciosamente;
- mapa azul/transparente;
- funil DPOL restaurado;
- `host_final` padronizado mantendo `host` original.

Para reproduzir os arquivos derivados, use `01_rebuild_dashboard_data_v8.py`. Para gerar o pacote do Streamlit Cloud, use `03_prepare_streamlit_cloud_v8.py`.
