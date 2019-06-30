library(tidyverse)
library(SummarizedExperiment)
library(cowplot)

binding_data <- read_rds('data/binding_data.rds')
dep_res <- read_rds('data/dep_res.rds')
go_annotation <- read_rds('data/go_annotation.rds')
tfs <- c('CTCF', 'CEBPB', 'PPARG', 'RXRG', 'EP300', 'MED1')

ind <- map(binding_data, function(x){
  mcols(x)$name[mcols(x)$geneId %in% go_annotation$SYMBOL]
}) %>%
  unlist()

p1 <- dep_res %>%
  filter(factor %in% tfs) %>%
  mutate(contrast = case_when(
    contrast == 'early_vs_non' ~ 'Early vs Non',
    contrast == 'late_vs_non' ~ 'Late vs Non',
    contrast == 'late_vs_early' ~ 'Late vs Early'
  ),
  contrast = factor(contrast, levels = c('Early vs Non', 'Late vs Non', 'Late vs Early'))) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(color = 'darkgray', alpha = .5) +
  geom_vline(xintercept = c(-1,1), lty = 2) +
  geom_hline(yintercept = 5, lty = 2) +
  facet_grid(contrast ~ factor, scales = 'free_y') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null")) +
  labs(y = "P-value (-Log 10)",
       x = 'Fold-Change (Log 2)')

p2 <- dep_res %>%
  filter(row %in% ind, factor %in% tfs) %>%
  mutate(contrast = case_when(
    contrast == 'early_vs_non' ~ 'Early vs Non',
    contrast == 'late_vs_non' ~ 'Late vs Non',
    contrast == 'late_vs_early' ~ 'Late vs Early'
  ),
  contrast = factor(contrast, levels = c('Early vs Non', 'Late vs Non', 'Late vs Early'))) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(color = 'darkgray', alpha = .5) +
  geom_vline(xintercept = c(-1,1), lty = 2) +
  geom_hline(yintercept = 5, lty = 2) +
  facet_grid(contrast ~ factor, scales = 'free_y') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null")) +
  labs(y = "P-value (-Log 10)",
       x = 'Fold-Change (Log 2)')

plot_grid(p1, p2,
          ncol = 1,
          scale = .95,
          labels = 'AUTO',
          label_fontface = 'plain',
          label_size = 10) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/volcanos_binding.png',
         width = 20, height = 24, units = 'cm')
