library(tidyverse)
library(SummarizedExperiment)

peak_counts <- read_rds('data/peak_counts.rds')

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
  mutate(annotation = str_split(annotation, ' \\(', simplify = TRUE)[, 1],
         annotation = ifelse(annotation %in% c("3' UTR", "5' UTR", 'Promoter'), annotation, 'Other')) %>%
  na.omit() %>%
  group_by(factor, group, geneId, annotation) %>%
  mutate(n = n()) %>%
  ungroup()

write_rds(factor_targets, 'data/factor_targets.rds')

  
