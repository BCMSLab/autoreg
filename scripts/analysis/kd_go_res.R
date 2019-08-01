# loading requried librar
library(tidyverse)
library(goseq)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)

# loading data
cebpb_kd_res <- read_rds('data/cebpb_kd_res.rds')
pparg_kd_res <- read_rds('data/pparg_kd_res.rds')

# extract differentially expressed genes
deg <- list()

deg$cebpb <- cebpb_kd_res %>%
  na.omit() %>%
  mutate(sig = padj < .2) %>%
  with(split(., time)) %>%
  map(function(x) {
    vec <- as.integer(x$sig)
    names(vec) <- x$row
    vec
  })

deg$pparg <- pparg_kd_res %>%
  na.omit() %>%
  mutate(sig = adj.P.Val < .2) %>%
  with(split(., time)) %>%
  map(function(x) {
    vec <- as.integer(x$sig)
    x$gene_id[vec == 1]
  })

# prepare gene to gene onology terms
gene_go <- AnnotationDbi::select(org.Mm.eg.db,
                                 unique(c(cebpb_kd_res$row, pparg_kd_res$gene_id)),
                                 'GO',
                                 'SYMBOL') %>%
  dplyr::select(gene = SYMBOL, cat = GO) %>%
  unique()

# perform gene set enrichment analysis
go_test <- list()

go_test$cebpb <- map(deg$cebpb, function(y) {
    nullp(y, 'mm10', 'geneSymbol', plot.fit = FALSE) %>%
      goseq(gene2cat = gene_go)
  }) %>%
    bind_rows(.id = 'time') %>%
    as_tibble() %>%
    mutate(ratio = numDEInCat/numInCat) %>%
    dplyr::select(time, go_id = category, ratio, pvalue = over_represented_pvalue)

go_test$pparg <- map(deg$pparg, function(x) {
  enricher(x,
           TERM2GENE = dplyr::select(gene_go, cat, gene),
           minGSSize = 1, 
           universe = unique(gene_go$gene),
           pAdjustMethod = 'none')@result
}) %>%
  bind_rows(.id = 'time') %>%
  as_tibble() %>%
  separate(GeneRatio, c('de', 'bg')) %>%
  mutate(ratio = as.numeric(de)/as.numeric(bg)) %>%
  dplyr::select(time, go_id = ID, ratio, pvalue)

# tidy data
kd_go_res <- go_test %>%
  bind_rows(.id = 'factor') %>%
  mutate(term = AnnotationDbi::select(GO.db,
                                      go_id, 
                                      'TERM', 
                                      'GOID')$TERM)

write_rds(kd_go_res, 'data/kd_go_res.rds')

