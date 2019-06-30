library(tidyverse)
library(reshape2)
library(SummarizedExperiment)
library(cowplot)

go_annotation <- read_rds('data/go_annotation.rds')
peak_counts <- read_rds('data/peak_counts.rds')

# subset the object
tfs <- c('CTCF', 'CEBPB', 'PPARG', 'RXRG', 'EP300', 'MED1')
hms <- c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9me3')

ind <- mcols(peak_counts)$geneId %in% unique(go_annotation$SYMBOL)

se <- peak_counts[ind, ]
se <- se[, !is.na(se$factor)]
se <- se[, !is.na(se$group)]
se <- se[, se$factor %in% c(tfs, hms)]


pd <- tibble(Var2 = se$id,
             factor = se$factor,
             time = se$time, 
             group = se$group)
fd <- tibble(Var1 = rownames(se),
             gene_id = mcols(se)$geneId,
             anno = str_split(mcols(se)$annotation, ' \\(', simplify = TRUE)[, 1]) %>%
  mutate(anno = ifelse(anno %in% c("3' UTR", "5' UTR", 'Promoter'), anno, 'Other'),
         anno = factor(anno, levels = c('Promoter', "3' UTR", "5' UTR", 'Other')))

table(fd$anno)
counts <- assay(se) %>%
  melt() %>%
  left_join(pd) %>%
  left_join(fd) %>%
  group_by(factor, gene_id, anno, Var1, time, group) %>%
  summarise(count = mean(value)) %>%
  ungroup()

factor_cofactor <- counts %>%
  filter(factor %in% tfs) %>%
  spread(factor, count) %>%
  select(-CTCF) %>%
  gather(factor, count_factor, PPARG, CEBPB) %>%
  gather(cofactor, count_cofactor, EP300, MED1, RXRG) %>%
  na.omit() %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c('non', 'early', 'late')))

p1 <- factor_cofactor %>%
  ggplot(aes(y = scale(log2(count_factor+1)),
             x = scale(log2(count_cofactor+1)))) +
  geom_point(color = 'gray', alpha = .5) +
  geom_smooth(aes(color = group), method = lm, se = FALSE) +
  facet_grid(factor ~ cofactor) +
  theme_bw() +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  lims(x = c(0, 5), y = c(0,5)) +
  labs(x = 'Co-factor Counts (Log 2)',
       y = 'Factor Counts (Log 2)',
       color = 'Stage of Differentiation')

p2 <- factor_cofactor %>%
  ggplot(aes(y = scale(log2(count_factor+1)),
             x = scale(log2(count_cofactor+1)))) +
  geom_point(color = 'gray', alpha = .5) +
  geom_smooth(aes(color = anno), method = lm, se = FALSE) +
  facet_grid(factor ~ cofactor) +
  theme_bw() +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  lims(x = c(0, 5), y = c(0,5)) +
  labs(x = 'Co-factor Counts (Log 2)',
       y = 'Factor Counts (Log 2)',
       color = 'Genomic Region')

factor_hm <- counts %>%
  filter(factor %in% c('CEBPB', 'PPARG', 'H3K27ac', 'H3K4me1', 'H3K4me3')) %>%
  spread(factor, count) %>%
  gather(factor, count_factor, PPARG, CEBPB) %>%
  gather(hm, count_hm, H3K27ac, H3K4me1, H3K4me3) %>%
  na.omit() %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c('non', 'early', 'late')))

p3 <- factor_hm %>%
  ggplot(aes(y = scale(log2(count_factor+1)),
             x = scale(log2(count_hm+1)))) +
  geom_point(color = 'gray', alpha = .5) +
  geom_smooth(aes(color = group), method = lm, se = FALSE) +
  facet_grid(factor ~ hm) +
  theme_bw() +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  lims(x = c(0, 5), y = c(0,5)) +
  labs(x = 'Histone Marker Counts (Log 2)',
       y = 'Factor Counts (Log 2)',
       color = 'Stage of Differentiation')

p4 <- factor_hm %>%
  ggplot(aes(y = scale(log2(count_factor+1)),
             x = scale(log2(count_hm+1)))) +
  geom_point(color = 'gray', alpha = .5) +
  geom_smooth(aes(color = anno), method = lm, se = FALSE) +
  facet_grid(factor ~ hm) +
  theme_bw() +
  theme(legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  lims(x = c(0, 5), y = c(0,5)) +
  labs(x = 'Histone Marker Counts (Log 2)',
       y = 'Factor Counts (Log 2)',
       color = 'Genomic Region')

plot_grid(p1, p3, p2, p4, 
          nrow = 2,
          labels = 'AUTO',
          scale = .95,
          label_fontface = 'plain',
          label_size = 10) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/factor_cofac_hm.png',
         width = 25, height = 22, units = 'cm')
