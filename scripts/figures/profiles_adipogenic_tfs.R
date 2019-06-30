library(tidyverse)
library(reshape2)
library(cowplot)
library(SummarizedExperiment)

gene_counts <- read_rds('data/gene_counts.rds')
occupancy <- read_rds('data/factor_occupancy.rds')

targets <- c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')

pd <- tibble(time = gene_counts$time,
             group = gene_counts$group,
             id = gene_counts$id)
profiles <- list('Gene Expression' = assay(gene_counts)[rownames(gene_counts) %in% targets,] %>%
                   melt() %>%
                   setNames(c('geneId', 'id', 'count')) %>%
                   left_join(pd),
                 'Transcription Activity' = occupancy %>%
                   filter(geneId %in% targets) %>%
                   dplyr::select(-factor)) %>%
  bind_rows(.id = 'count_type') %>%
  mutate(group = factor(group, levels = c('non', 'early', 'late')))

(profiles %>%
    na.omit() %>%
    mutate(count = log2(count + 1)) %>%
    group_by(count_type, geneId, group) %>%
    mutate(ave = mean(count), sd = sd(count)) %>%
    ggplot(aes(x = group, y = count)) +
    geom_jitter(width = .2, color = 'darkgray', alpha = .5) +
    geom_point(aes(y = ave), color = 'red') +
    geom_linerange(aes(ymin = ave - sd, ymax = ave + sd), color = 'red') +
    facet_grid(count_type~geneId) +
    theme_bw() + 
    labs(x = 'Stage of Differentiation',
         y = 'Log2 Count') +
    theme(strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0,"null"))) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/profiles_adipogenic_tfs.png',
         width = 24, height = 10, units = 'cm')
