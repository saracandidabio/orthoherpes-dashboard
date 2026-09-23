# Política de plotagem e preservação — v7

1. **Fonte arquivada primeiro.** `records_dashboard.tsv.gz` e `cds_dashboard.tsv.gz` têm SHA256 preservado.
2. **Camadas derivadas não substituem a fonte.** Ex.: `country_canonical`, `host_final`, agregados, classes de missing.
3. **Sem inferência de país ambíguo.** `Korea` e `USSR` permanecem sem conversão moderna automática.
4. **Sem normalização automática de hospedeiro.** Crosswalk futuro exige revisão manual.
5. **Sem renomear hypothetical protein por similaridade.** `annotation_clue` é pista, não identificação confirmada.
6. **Sem tratar todo accession como genoma.** Comparações usam `accession`/registro.
7. **Gráficos podem limitar renderização, nunca o dado disponível.** Tabelas/downloads mantêm todas as categorias/linhas do recorte.
8. **Taxas exigem denominador explícito.** A página de hypothetical usa mínimo padrão de 20 CDS e mostra suporte de accessions.
9. **Validação NCBI é snapshot externo.** Divergência nunca é aplicada automaticamente ao arquivo local.
