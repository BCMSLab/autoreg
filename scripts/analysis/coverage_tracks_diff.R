library(rtracklayer)
library(EnrichedHeatmap)
library(tidyverse)
library(reshape2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

peak_counts <- read_rds('data/peak_counts.rds')
pd <- tibble(factor = peak_counts$factor,
             sample = peak_counts$id,
             time = peak_counts$time)

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

coverage_tracks_diff <- map(list(tf = factor_track_signal,
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

deg_res <- read_rds('data/deg_res.rds') %>%
  filter(contrast != 'late_vs_early') %>%
  dplyr::select(group = contrast, SYMBOL = row, fc = log2FoldChange) %>%
  mutate(group = ifelse(group == 'early_vs_non', 'early', 'late'),
         dir = case_when(fc > 1 ~ 'UP',
                         fc < -1 ~ 'Down',
                         TRUE ~ 'None'))

gene_diff <- as_tibble(mcols(tss)) %>%
  mutate_all(as.character) %>%
  left_join(goi, by = c('gene_id'='ENTREZID')) %>%
  dplyr::select(-gene_id) %>%
  left_join(deg_res)

coverage_tracks_diff <- map(coverage_tracks_diff, function(x) {
  x %>%
    bind_rows(.id = 'sample') %>%
    inner_join(gene_diff, by = c('Var1'='tx_name')) %>%
    inner_join(pd) %>%
    group_by(factor, group, dir, Var2) %>%
    summarise(value = mean(value, na.rm = TRUE)) %>%
    mutate(Var2 = str_sub(Var2, start = 2) %>% as.numeric(),
           Var2 = ifelse(duplicated(Var2), Var2, Var2-60)) %>%
    na.omit() %>%
    ungroup()
})

write_rds(coverage_tracks_diff, 'data/coverage_tracks_diff.rds')
