library(tidyverse)
library(reshape2)
library(cowplot)
library(SummarizedExperiment)

gene_counts <- read_rds('data/gene_counts.rds')
occupancy <- read_rds('data/factor_occupancy.rds')

go_annotation <- read_rds('data/go_annotation.rds')
tf_annotation <- read_rds('data/tf_annotation.rds')

profiles <- list(Autophagy = intersect(go_annotation$SYMBOL, tf_annotation$SYMBOL),
     Adipogenic = c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')) %>%
  map(function(x) {
    pd <- tibble(time = gene_counts$time,
                 group = gene_counts$group,
                 id = gene_counts$id)
    list('Gene Expression' = assay(gene_counts)[rownames(gene_counts) %in% x,] %>%
           melt() %>%
           setNames(c('geneId', 'id', 'count')) %>%
           left_join(pd),
         'Transcription Activity' = occupancy %>%
           filter(geneId %in% x) %>%
           dplyr::select(-factor)) %>%
      bind_rows(.id = 'count_type')
  }) %>%
  bind_rows(.id = 'tf_type') %>%
  mutate(group = factor(group, levels = c('non', 'early', 'late')))


(map(c('Autophagy', 'Adipogenic'), function(x) {
  profiles %>%
    filter(tf_type == x) %>%
    na.omit() %>%
    mutate(count = log2(count + 1)) %>%
    group_by(count_type, geneId, group) %>%
    mutate(ave = mean(count), sd = sd(count)) %>%
    ggplot(aes(x = group, y = count)) +
    geom_jitter(width = .2, color = 'darkgray', alpha = .5) +
    geom_point(aes(y = ave), color = 'red') +
    geom_linerange(aes(ymin = ave - sd, ymax = ave + sd), color = 'red') +
    facet_grid(count_type~geneId) +
    lims(y = c(0, 15)) +
    theme_bw() + 
    labs(x = 'Stage of Differentiation',
         y = 'Log2 Count') +
    theme(strip.background = element_blank(),
          panel.grid = element_blank())
}) %>%
  plot_grid(plotlist = .,
            ncol = 1,
            labels = 'AUTO',
            label_fontface = 'plain',
            label_size = 10,
            scale = .95)) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/target_tf_profiles.png',
         width = 24, height = 20, units = 'cm')
