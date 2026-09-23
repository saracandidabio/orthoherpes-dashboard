# Política de dados e plotagem — v8

1. `records_dashboard.tsv.gz` e `cds_dashboard.tsv.gz` são fontes arquivadas e não são sobrescritas.
2. `species_final` não é alterado pela interface.
3. `country_canonical` existe apenas para agregação/mapa; `country_final` permanece disponível.
4. `host_final` é uma camada derivada de `host`, com regras aprovadas e registradas em `host_crosswalk.tsv`.
5. `Homo spaiens` é corrigido para `Homo sapiens`; NCBI Taxonomy mantém `Homo sapiens` como nome atual, TaxID 9606.
6. Metadados após o primeiro `;` no campo host não entram em `host_final`.
7. Capitalização é padronizada sem criar nova identidade taxonômica: o par é usado apenas quando já ocorre como rótulo independente no dataset.
8. "Todas" significa todas as categorias. Para grande cardinalidade o dashboard troca barras por treemap/sunburst, evitando truncamento silencioso.
9. Registros NCBI/GenBank no dashboard representam registros no banco, não prevalência biológica.
