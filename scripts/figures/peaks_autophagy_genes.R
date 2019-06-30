library(tidyverse)
library(reshape2)
library(SummarizedExperiment)

binding_data <- read_rds('data/binding_data.rds')

tf <- c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')

targets <- c('Atg4b', 'Ulk1', 'Map1lc3a', 'Map1lc3b', 'Sqstm1', 'Becn1')


df <- binding_data %>%
  map(function(x) {
    pd <- tibble(Var2 = x$id,
                 group = x$group,
                 time = x$time)
    fd <- tibble(Var1 = rownames(x),
                 geneId = mcols(x)$geneId)
    assay(x) %>%
      melt() %>%
      left_join(pd) %>%
      left_join(fd)
  }) %>%
  bind_rows(.id = 'factor')

(df %>%
    filter(geneId %in% targets,
           !factor %in% c('POLR2A'),
           !is.na(group)) %>%
    mutate(group = factor(group, levels = c('non', 'early', 'late'))) %>%
    group_by(group, factor, geneId, Var1) %>%
    summarise(value = log2(mean(value) + 1)) %>%
    group_by(group, factor, geneId) %>%
    mutate(ave = mean(value), sd = sd(value)) %>%
    ggplot(aes(x = group, y = value)) +
    geom_jitter(width = .2, color = 'darkgray', alpha = .5) +
    geom_point(aes(y = ave), color = 'red') +
    geom_linerange(aes(ymin = ave - sd, ymax = ave + sd), color = 'red') +
    facet_grid(factor~geneId) +
    theme_bw() + 
    labs(x = 'Stage of Differentiation',
         y = 'Log2 Count in Peaks') +
    theme(strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0,"null"))) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/peaks_autophagy_genes.png',
         width = 24, height = 20, units = 'cm')
