library(tidyverse)
library(SummarizedExperiment)

transformed_counts <- read_rds('data/transformed_counts.rds')
go_annotation <- read_rds('data/go_annotation.rds')
tf_annotation <- read_rds('data/tf_annotation.rds')
tf <- c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')

(list('All Genes' = rownames(transformed_counts),
     'Autophagy Genes' = rownames(transformed_counts) %in% go_annotation$SYMBOL,
     'Autophagy TF' = rownames(transformed_counts) %in% intersect(go_annotation$SYMBOL, tf_annotation$SYMBOL),
     "Adipogenic TF" = rownames(transformed_counts) %in% tf) %>%
  map(function(x) {
    cmdscale(dist(t(assay(transformed_counts)[x,]))) %>%
      as.data.frame() %>%
      mutate(group = transformed_counts$group)
  }) %>%
  bind_rows(.id = 'type') %>%
  mutate(type = factor(type,
                       levels = c('All Genes', 'Autophagy Genes', 'Autophagy TF', 'Adipogenic TF'))) %>%
  ggplot(aes(x = V1, y = V2, color = group)) +
  geom_point() +
  facet_wrap(~type, nrow = 1, scales = 'free') +
  theme_bw() +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = 'MDS1', y = 'MDS2', color = 'Stage of Differentiation')) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/mds_gene_group.png',
         height = 7, width = 20, units = 'cm')
