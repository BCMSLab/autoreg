library(tidyverse)

hm_occupancy <- read_rds('data/hm_occupancy.rds')
hm <- c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9me3')

# ddcor <- read_rds('data/dgca.rds')
# (ddcor %>%
#   dplyr::select(Gene2, geneId = Gene1, group, cor, pValDiff_adj) %>%
#   unique() %>%
#   filter(pValDiff_adj < .1) %>%
#   mutate(group = str_split(group, '_', simplify = TRUE)[, 1],
#          group = factor(group, levels = c('non', 'early', 'late')),
#          dir = ifelse(cor > 0, 'Positive', 'Negative')) %>%
#   inner_join(hm_occupancy) %>%
#   filter(factor %in% hm) %>%
#   ggplot(aes(x = group, y = log2(count), fill = dir)) +
#   geom_boxplot() +
#   facet_grid(factor~Gene2) +
#   theme_bw() +
#   theme(strip.background = element_blank(),
#         panel.grid = element_blank(),
#         legend.position = 'top') +
#   labs(x = '',
#        y = "Histone Tag Counts (Log 2)",
#        fill = '') +
#   scale_x_discrete('Stage of Differentiation')) %>%
#   ggsave(plot = .,
#          'manuscript/figures/modification_factor_direction.png',
#          width = 20, height = 20, units = 'cm')

deg_res <- read_rds('data/deg_res.rds') %>%
  filter(contrast != 'late_vs_early') %>%
  dplyr::select(group = contrast, geneId = row, fc = log2FoldChange) %>%
  mutate(group = ifelse(group == 'early_vs_non', 'early', 'late'),
         dir = case_when(fc > 1 ~ 'UP',
                         fc < -1 ~ 'Down',
                         TRUE ~ 'None'))

(hm_occupancy %>%
  inner_join(deg_res) %>%
  ggplot(aes(x = group, y = log2(count + 1), fill = dir)) +
  geom_boxplot() +
  facet_wrap(~factor, nrow = 1, scales = 'free_x') +
  theme_bw() +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null"))+
  labs(y = "HM Marker Tags (Log2)",
       x = 'Stage of Differentiation',
       fill = 'Regulation')) %>%
  ggsave(plot = .,
         'manuscript/figures/modification_factor_direction.png',
         width = 20, height = 7, units = 'cm')
