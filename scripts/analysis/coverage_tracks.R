library(rtracklayer)
library(EnrichedHeatmap)
library(tidyverse)
library(reshape2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

peak_counts <- read_rds('data/peak_counts.rds')
pd <- tibble(factor = peak_counts$factor,
             sample = peak_counts$id,
             group = peak_counts$group,
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
                 downstream = 1)

coverage_tracks <- map(list(tf = factor_track_signal,
                            hm = hm_track_signal), function(t) {
                              map(t, function(x) {
                                normalizeToMatrix(x,
                                                  value_column = 'score',
                                                  target = tss, 
                                                  extend = 3000,
                                                  mean_mode = 'w0',
                                                  w = 50) %>%
                                  melt()
                              }) %>%
                                bind_rows(.id = 'sample') %>%
                                inner_join(pd) %>%
                                group_by(factor, group, Var2) %>%
                                summarise(value = mean(value)) %>%
                                mutate(Var2 = str_sub(Var2, start = 2) %>% as.numeric(),
                                       Var2 = ifelse(duplicated(Var2), Var2, Var2-60))
                            })

write_rds(coverage_tracks, 'data/coverage_tracks.rds')
