library(tidyverse)

coverage_tracks <- read_rds('data/coverage_tracks.rds')

(coverage_tracks$tf %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c('non', 'early', 'late'))) %>%
  ggplot(aes(x = Var2, y = value, color = group)) +
  geom_line() +
  geom_vline(xintercept = 0, lty = 2, color = 'black') +
  facet_wrap(~factor, nrow = 1, scales = 'free_y') +
  theme_bw() +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = '', y = 'Score', color = 'Stage of Differentiation') +
  scale_x_continuous(breaks = c(-60, 0, 60),
                     labels = c('-3000', 'TSS', '3000'))) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/tf_tracks.png',
         width = 24, height = 7, units = 'cm')
