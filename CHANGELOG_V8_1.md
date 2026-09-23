# Orthoherpesviridae Dataset Explorer — v8.1

Data: 2026-09-23

## Alterações solicitadas

### Países

A coluna-fonte `country_final` continua intacta. O dashboard usa a camada derivada
`country_current`/`country_canonical` para filtros, rankings e mapa.

Normalizações atuais aplicadas apenas quando inequívocas:

- `Turkey` → `Türkiye`
- `Turkiye` → `Türkiye`
- `Zaire` → `Democratic Republic of the Congo`
- `Macedonia` → `North Macedonia`
- `Cape Verde` → `Cabo Verde`

`Korea` e `USSR` são preservados porque não há correspondência moderna 1:1 segura a partir
do rótulo isolado.

Arquivos de auditoria:

- `country_crosswalk.tsv`
- `COUNTRY_NORMALIZATION_AUDIT.tsv`
- `country_plot_exceptions.tsv`

Resultado: 5 rótulos alterados, afetando 645 registros. Os hashes dos arquivos-fonte
`records_dashboard.tsv.gz` e `cds_dashboard.tsv.gz` permanecem os mesmos.

### Genes / CDS

Foi criado `gene_context_legend.tsv.gz`, com uma linha para cada gene observado no dataset.
A legenda mostra:

- produtos/proteínas mais frequentes no próprio dataset;
- herpesvírus/táxons mais associados;
- combinações gene × vírus × produto mais frequentes;
- nota NCBI curada para símbolos frequentes/ambíguos (UL30, UL55, UL44, UL54, UL23,
  UL40, UL27, UL42, UL5, UL52, US3, US6 e US8);
- referência NCBI usada na nota curada.

As notas NCBI curadas aparecem imediatamente abaixo do gráfico; o contexto completo observado no dataset fica em um expander logo abaixo. Importante:
símbolos UL/US não têm significado universal entre todos os herpesvírus; por exemplo,
UL55 é gB em HCMV e uma proteína nuclear em HSV-1.

### Cores

O preto foi removido da paleta categórica. `HSV-1` recebeu cor violeta explícita
`#8A5CF6`, visível em temas claro e escuro. Os grupos de interesse agora usam um mapa de
cores fixo para manter consistência entre gráficos.

### Cabeçalho

Foi adicionado o bloco de autoria solicitado, com formação, Currículo Lattes, LinkedIn,
e-mail e a data de obtenção do dataset no NCBI (`21/09/2026`).

## Integridade

- `records_dashboard.tsv.gz`: SHA256 preservado.
- `cds_dashboard.tsv.gz`: SHA256 preservado.
- 77.099 registros e 446.993 CDS continuam inalterados.
- 835 `species_final` resolvidos + `Unresolved_species` permanecem inalterados.

