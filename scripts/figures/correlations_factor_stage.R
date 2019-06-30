library(tidyverse)
library(reshape2)

# dgca <- read_rds('data/dgca.rds')
# 
# (dgca %>%
#     dplyr::select(Gene1, Gene2, cor, group, pValDiff_adj) %>%
#     unique() %>%
#     mutate(group = factor(group, levels = c('non_cor', 'early_cor', 'late_cor'))) %>%
#     filter(pValDiff_adj < .1) %>%
#     ggplot(aes(x = group, y = cor)) +
#     geom_boxplot() +
#     facet_wrap(~Gene2, nrow = 2) +
#     theme_bw() +
#     theme(strip.background = element_blank(),
#           panel.grid = element_blank()) +
#     labs(x = '', y = "Pearson's Correlation") +
#     scale_x_discrete('Stage of Differentiation',
#                      labels = c('non', 'early', 'late'))) %>%
#   ggsave(plot = .,
#          'manuscript/figures/correlations_factor_stage.png',
#          width = 20, height = 16, units = 'cm')

# with deg
deg_res <- read_rds('data/deg_res.rds') %>%
  filter(contrast != 'late_vs_early') %>%
  dplyr::select(group = contrast, Gene1 = row, fc = log2FoldChange) %>%
  mutate(group = ifelse(group == 'early_vs_non', 'early', 'late'),
         dir = case_when(fc > 1 ~ 'UP',
                         fc < -1 ~ 'Down',
                         TRUE ~ 'None'))

ddcor <- read_rds('data/ddcor.rds') %>%
  map(function(x) {
    list(early = filter(x$early) %>% dplyr::select(Gene1, Gene2, cor = early_cor),
         late  = filter(x$late) %>% dplyr::select(Gene1, Gene2, cor = late_cor)) %>%
      bind_rows(.id = 'group')
  }) %>%
  bind_rows()

(left_join(ddcor, deg_res) %>%
  na.omit() %>%
  ggplot(aes(x = group, y = abs(cor), fill = dir)) +
  geom_boxplot() +
  lims(y = c(0, 1)) +
  facet_wrap(~Gene2, nrow = 1) +
  theme_bw() +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null"))+
  labs(y = "Pearson's Correlation",
       x = 'Stage of Differentiation',
       fill = 'Regulation')) %>%
  ggsave(plot = .,
         'manuscript/figures/correlations_factor_stage.png',
         width = 20, height = 8, units = 'cm')
  
