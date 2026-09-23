# Changelog v8.0.0 — 23/09/2026

## Interface
- Restaurada a navegação por barra superior com `st.tabs`.
- As abas usam `on_change="rerun"` e `.open`; apenas a aba ativa executa.
- Reduzidas notas repetitivas e textos explicativos na interface.
- O Fluxo DPOL voltou ao funil horizontal/vertical do modelo anterior, com a mesma lógica de estágios e contagens.

## Explorador taxonômico
- Removido o limite artificial de 20 táxons no multiselect.
- Adicionada a opção **Usar todos os táxons do recorte**.
- Quando todos os táxons são usados, o treemap inclui todo o conjunto elegível e a tabela completa é exibida sem truncamento taxonômico.
- Em controles Top N, **Todas** passa a significar todas as categorias. Listas muito longas mudam para treemap/sunburst em vez de truncar silenciosamente.

## Geografia
- Mapa redesenhado com fundo transparente, gradação em tons de azul e linhas de fronteira/costa em azul-ciano de alto contraste.
- Mantida a política conservadora de país: `country_final` original não é sobrescrito.

## Hospedeiros — padronização aprovada
O campo `host` original é preservado. O dashboard usa `host_final` derivado por `host_crosswalk.tsv`.

Regras aplicadas:
1. `Homo spaiens` → `Homo sapiens`.
2. Quando um rótulo contém `;`, `host_final` considera somente o texto antes do primeiro `;`.
3. Os dois primeiros termos são padronizados para `Gênero espécie` somente quando o mesmo par de dois termos existe como rótulo independente no próprio dataset. Isso evita inventar uma identificação taxonômica nova.

Resultado desta versão:
- 1.249 rótulos originais não vazios.
- 972 rótulos distintos em `host_final`.
- 372 rótulos alterados por pelo menos uma regra aprovada.
- 2.527 registros associados a rótulos alterados.
- 285 rótulos tinham metadados após `;` (850 registros).
- 409 registros `Homo spaiens` passam a `host_final=Homo sapiens`.

## Integridade
- `records_dashboard.tsv.gz` continua com SHA256 `fee5fcd76e59ce995c5e6954fd562ee7baeb8a6fff618ca3610b1b06891568b0`.
- `cds_dashboard.tsv.gz` continua com SHA256 `2c47da5eb2753becf77d476b31d14c18135c8e53b7003b47610f6050c8a013fa`.
- As correções de hospedeiro são uma camada derivada e auditável; os arquivos-fonte não foram sobrescritos.
