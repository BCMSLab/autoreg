library(tidyverse)

# dgca <- read_rds('data/dgca.rds')
# occupancy <- read_rds('data/factor_occupancy.rds') %>%
#   filter(factor != 'POLR2A')
# 
# (dgca %>%
#   dplyr::select(factor = Gene2, geneId = Gene1, group, cor, pValDiff_adj) %>%
#   unique() %>%
#   filter(pValDiff_adj < .1) %>%
#   mutate(factor = toupper(factor),
#          group = str_split(group, '_', simplify = TRUE)[, 1],
#          group = factor(group, levels = c('non', 'early', 'late')),
#          dir = ifelse(cor > 0, 'Positive', 'Negative')) %>%
#   inner_join(occupancy) %>%
#   ggplot(aes(x = group, y = log2(count), fill = dir)) +
#   geom_boxplot() +
#   facet_wrap(~factor, nrow = 2) +
#   theme_bw() +
#   theme(strip.background = element_blank(),
#         panel.grid = element_blank(),
#         legend.position = 'top') +
#   labs(x = '',
#        y = "TF Occupancy (Log 2)",
#        fill = '') +
#   scale_x_discrete('Stage of Differentiation')) %>%
#   ggsave(plot = .,
#          'manuscript/figures/occupancy_factor_direction.png',
#          width = 20, height = 16, units = 'cm')

occupancy <- read_rds('data/factor_occupancy.rds') 

deg_res <- read_rds('data/deg_res.rds') %>%
  filter(contrast != 'late_vs_early') %>%
  dplyr::select(group = contrast, geneId = row, fc = log2FoldChange) %>%
  mutate(group = ifelse(group == 'early_vs_non', 'early', 'late'),
         dir = case_when(fc > 2 ~ 'UP',
                         fc < -2 ~ 'Down',
                         TRUE ~ 'None'))

(occupancy %>%
  filter(factor != 'POLR2A') %>%
  group_by(factor, geneId, time, group, time) %>%
  summarise(count = mean(count)) %>%
  inner_join(deg_res) %>%
  ggplot(aes(x = group, y = log2(count+1), fill = dir)) +
  geom_boxplot() +
  facet_grid(~factor, scales = 'free') +
  theme_bw() +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null"))+
  labs(y = "TF Occupancy (Log2)",
       x = 'Stage of Differentiation',
       fill = 'Regulation')) %>%
  ggsave(plot = .,
       'manuscript/figures/occupancy_factor_direction.png',
       width = 20, height = 8, units = 'cm')

# occupancy %>%
#   inner_join(deg_res) %>%
#   ggplot(aes(x = group, y = log2(count+1), fill = dir)) +
#   geom_boxplot() +
#   facet_grid(annotation~factor, scales = 'free') +
#   theme_bw() +
#   theme(legend.position = 'top',
#         panel.grid = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing = unit(0,"null"))+
#   labs(y = "TF Occupancy (Log2)",
#        x = 'Stage of Differentiation',
#        fill = 'Regulation')
#   ggsave(plot = .,
#          'manuscript/figures/occupancy_factor_direction_anno.png',
#          width = 20, height = 22, units = 'cm')
# 
#   