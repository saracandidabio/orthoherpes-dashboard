# Orthoherpesviridae Dataset Explorer — v7.0.0

Data: 2026-09-23

## Objetivo

Consolidar em uma única versão as correções de desempenho, visualização, casos de borda,
coerência de filtros, interpretação biológica, proveniência e auditoria discutidas durante a
construção do dashboard, sem reescrever silenciosamente os dados arquivados.

## Integridade do dado-fonte

Os dois arquivos principais usados para reconstrução foram preservados byte a byte:

- `records_dashboard.tsv.gz`: SHA256 `fee5fcd76e59ce995c5e6954fd562ee7baeb8a6fff618ca3610b1b06891568b0`
- `cds_dashboard.tsv.gz`: SHA256 `2c47da5eb2753becf77d476b31d14c18135c8e53b7003b47610f6050c8a013fa`

Contagens validadas:

- 77.099 registros;
- 77.099 `accession_version` únicos;
- 446.993 CDS;
- 48.994 accessions com pelo menos uma CDS;
- 835 `species_final` resolvidos;
- 441 registros `Unresolved_species`;
- 25.086 CDS classificadas como hypothetical.

## 1. Desempenho

- Mantido `st.navigation`: somente a página selecionada executa.
- `st.cache_resource` para tabelas de leitura compartilhadas.
- Colunas repetitivas são carregadas como `category`.
- CDS foram separadas em:
  - `cds_core.tsv.gz`: campos usados em gráficos/comparações;
  - `cds_details.tsv.gz`: textos e campos de uso eventual (`note`, `location`, datas etc.).
- `note` e outros detalhes só são carregados quando o usuário executa busca ou solicita a prévia detalhada.
- Tabelas agregadas globais foram pré-calculadas para genes, produtos, pistas, classes de anotação e resumo de proteínas hipotéticas por táxon.
- Quando não há filtros, páginas de Genes/CDS e Proteínas hipotéticas usam esses agregados antes de carregar 446.993 CDS.
- Downloads grandes continuam sob demanda.
- Busca usa `regex=False` e formulários, evitando rerun a cada tecla.
- Conversor Parquet atualizado para `records_dashboard`, `cds_core` e `cds_details`.

## 2. Bugs/casos de borda

- Sliders seguros quando existe somente um item.
- Filtro vazio gera aviso explícito.
- Esquema obrigatório validado antes da análise.
- `Unresolved_species` é tratado consistentemente como espécie não resolvida.
- `missing` é tratado consistentemente como ausência, inclusive em `country_final`.
- Limite de categorias em gráficos impede figuras gigantes; tabelas/downloads preservam o universo completo.

## 3. Filtros

Filtros globais valem para as páginas analíticas. Seções pré-calculadas/auditáveis que ignoram filtros informam isso explicitamente.

## 4. Visualização e interpretação

- Paleta Okabe–Ito para categorias e Cividis para escalas quantitativas.
- Números formatados em padrão pt-BR.
- HSV-1 versus grupos raros pode ser visto em escala logarítmica.
- Comparação de proteínas usa `accession`/registro, nunca assume que todo accession seja genoma completo.
- Rótulos de comparação incluem suporte `n` de accessions.
- Opção para retirar hypothetical protein da matriz comparativa.
- Taxa de hypothetical usa mínimo padrão de 20 CDS por grupo, reduzível pelo usuário.
- Seleção de linhas/grupos é limitada quando necessário para evitar repetição de cores; linhas usam também traço/marcador quando aplicável.

## 5. Hospedeiros

Nenhum hospedeiro foi renomeado automaticamente.

- `host` continua sendo o valor original.
- Sem `host_crosswalk.tsv`, `host_final == host`.
- `host_crosswalk_candidates.tsv` é fornecido apenas como planilha de candidatos para curadoria manual futura.
- `host_group` só recebe conteúdo se existir crosswalk manualmente curado.

## 6. Países/mapa

Correções são somente camadas derivadas de plotagem:

- `Turkiye` pode ser agregado como `Turkey`;
- `Zaire` pode ser localizado no mapa na atual República Democrática do Congo, mantendo o original;
- `Korea` NÃO é convertido automaticamente em South Korea;
- `USSR` NÃO é atribuído automaticamente a um país moderno;
- `missing` é ausência.

## 7. Qualidade

`metadata_quality.tsv` foi recalculado com a mesma regra de missing do app.

País disponível:

- resumo anterior: 44.525;
- regra consistente atual: 44.523;
- diferença: 2 registros cujo `country_final` literal é `missing`.

Nenhum valor-fonte foi alterado para obter essa correção.

## 8. NCBI

Snapshot de validação fornecido em 2026-09-23:

- 176/176 accessions de foco: versão/accession = MATCH;
- 176/176: TaxID = MATCH;
- 176/176: `length_nt` = MATCH;
- nomes taxonômicos de foco exibidos no run = MATCH.

A comparação `definition` da v6 foi sensível a pontuação terminal. O validador v6.1 corrigido está incluído e distingue `MATCH_FORMATTING` de mudança substantiva. Nenhum `definition` local foi sobrescrito.

## 9. Proveniência

O dashboard mostra:

- versão dos dados;
- contagens do build;
- SHA256;
- relatório de validação local;
- exceções de país;
- snapshot NCBI;
- política de preservação.

## 10. Testes

- `py_compile`: OK;
- 20 testes originais de lógica: OK;
- 4 testes adicionais de integridade v7: OK;
- total: 24 testes aprovados.
