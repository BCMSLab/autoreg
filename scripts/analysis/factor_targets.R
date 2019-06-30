library(tidyverse)
library(SummarizedExperiment)

load('data/peak_counts.rda')

peak_counts$group <- cut(peak_counts$time,
                         breaks = c(-50, 0, 48, 240),
                         labels = c('non', 'early', 'late'))

peak_counts$group <- ifelse(is.na(peak_counts$group) & peak_counts$stage == 0,
                            'non',
                            as.character(peak_counts$group))

peak_counts$group <- ifelse(is.na(peak_counts$group) & peak_counts$stage == 3,
                            'late',
                            as.character(peak_counts$group))

tf <- c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')

tf_gsm <- data.frame(peak = colnames(peak_counts),
                     factor = peak_counts$factor,
                     time = peak_counts$time,
                     group = peak_counts$group) %>%
  filter(factor %in% toupper(tf))

factor_targets <- mcols(peak_counts) %>%
  as_tibble() %>%
  dplyr::select(peak, annotation, distanceToTSS, geneId) %>%
  unnest(peak) %>%
  unique() %>%
  left_join(tf_gsm) %>%
  na.omit()

write_rds(factor_targets, 'data/factor_targets.rds')
