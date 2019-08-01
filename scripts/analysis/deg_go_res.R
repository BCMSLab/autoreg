# loading required libraries
library(tidyverse)
library(org.Mm.eg.db)
library(goseq)

# loading data
deg_res <- read_rds('data/deg_res.rds')

# extract differentially espressed genes
deg <- deg_res %>%
  na.omit() %>%
  mutate(sig = (padj < .2) & (abs(log2FoldChange) > .5)) %>%
  with(split(., contrast)) %>%
  map(function(x) {
    vec <- as.integer(x$sig)
    names(vec) <- x$row
    vec
  })

# prepare gene to gene ontology terms
gene_go <- AnnotationDbi::select(org.Mm.eg.db,
                                 unique(deg_res$row),
                                 'GO',
                                 'SYMBOL') %>%
  dplyr::select(gene = SYMBOL, cat = GO) %>%
  unique()

# perform gene erichment analysis
deg_go_res <- map(deg, function(y) {
  nullp(y, 'mm10', 'geneSymbol', plot.fit = FALSE) %>%
    goseq(gene2cat = gene_go)
}) %>%
  bind_rows(.id = 'contrast') %>%
  as_tibble() 

write_rds(deg_go_res, 'data/deg_go_res.rds')

