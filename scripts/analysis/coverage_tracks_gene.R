library(rtracklayer)
library(EnrichedHeatmap)
library(tidyverse)
library(reshape2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

peak_counts <- read_rds('data/peak_counts.rds')
pd <- tibble(factor = peak_counts$factor,
             sample = peak_counts$id,
             time = peak_counts$time,
             group = peak_counts$group)

go_annotation <- read_rds('data/go_annotation.rds')
tf <- c('Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')

factor_track_signal <- read_rds('data/factor_track_signal.rds')
hm_track_signal <- read_rds('data/hm_track_signal.rds')

goi <- AnnotationDbi::select(org.Mm.eg.db,
                             c(unique(go_annotation$SYMBOL), tf),
                             'ENTREZID', 'SYMBOL')

tss <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,  
                 filter = list(gene_id = goi$ENTREZID),
                 upstream = 0,
                 downstream = 1,
                 columns = c('tx_name', 'gene_id'))

gene_trans <- as_tibble(tss) %>%
  dplyr::select(Var1 = tx_name, geneId = gene_id) %>%
  mutate_all(as.character) %>%
  unique()

coverage_tracks_gene <- map(list(tf = factor_track_signal,
                                 hm = hm_track_signal), function(t) {
                                   map(t, function(x) {
                                     normalizeToMatrix(x,
                                                       value_column = 'score',
                                                       target = tss, 
                                                       extend = 3000,
                                                       mean_mode = 'w0',
                                                       w = 50) %>%
                                       melt()
                                   }) 
                                 })

coverage_tracks_gene_df <- map(coverage_tracks_gene, function(x) {
  x %>%
    bind_rows(.id = 'sample') %>%
    inner_join(gene_trans) %>%
    inner_join(pd) %>%
    mutate(var2 = factor(Var2, levels = paste0(c('u', 'd'), c(60:1, 1:60)))) %>%
    na.omit() %>%
    group_by(Var1, Var2, geneId, factor, group) %>%
    summarise(value = mean(value)) %>%
    ungroup()
})

write_rds(coverage_tracks_gene_df, 'data/coverage_tracks_gene.rds')
