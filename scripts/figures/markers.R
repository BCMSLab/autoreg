library(tidyverse)
library(SummarizedExperiment)
library(reshape2)
library(cowplot)
library(gghighlight)

markers <- list(Adipogenesis = c('Cebpa', 'Pparg'),
                Lipogenesis = c('Lpl', 'Acly', 'Dgat', 'Elov6', 'Fasn', 'Scd'),
                Autophagy = c('Map1lc3b', 'Sqstm1', 'Becn1'))

gene_counts <- read_rds('data/gene_counts.rds')

pd <- tibble(Var2 = gene_counts$id,
             time = gene_counts$time,
             group = gene_counts$group)

(map(markers, function(x) {
  assay(gene_counts)[rownames(gene_counts) %in% x,] %>%
    melt %>%
    left_join(pd) %>%
    group_by(Var1, time) %>%
    summarise(count = mean(log2(value + 1))) %>%
    ggplot(aes(x = time, y = count, group = Var1)) +
    geom_line() +
    gghighlight(TRUE) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    lims(y = c(0,18)) +
    labs(x = 'Time (hr)',
         y = 'mRNA Level (Log 2)')
}) %>%
  plot_grid(plotlist = .,
            nrow = 1,
            scale = .95,
            labels = 'AUTO',
            label_fontface = 'plain',
            label_size = 10)) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/markers.png',
         width = 20, height = 7, units = 'cm')


# map(markers, function(x) {
#   assay(gene_counts)[rownames(gene_counts) %in% x,] %>%
#     melt %>%
#     left_join(pd) %>%
#     group_by(Var1, group) %>%
#     summarise(count = mean(log2(value + 1))) %>%
#     ggplot(aes(x = Var1, y = count, fill = group)) +
#     geom_col(position = 'dodge') +
#     theme(legend.position = 'top')
# }) %>%
#   plot_grid(plotlist = .,
#             nrow = 1,
#             scale = .95,
#             labels = 'AUTO',
#             label_fontface = 'plain',
#             label_size = 10)
