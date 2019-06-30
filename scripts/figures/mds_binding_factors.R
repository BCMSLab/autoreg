library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

binding_data <- read_rds('data/binding_data.rds')

mds <- binding_data %>%
  map(function(x) {
    dds <- DESeqDataSet(x, design = ~ group - 1)
    dds <- vst(dds)
    cmdscale(dist(t(assay(dds)))) %>%
      as.data.frame() %>%
      mutate(group = x$group)
  }) %>%
  bind_rows(.id = 'type') %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c('non', 'early', 'late')))

(mds %>%
  filter(type != 'POLR2A') %>%
  ggplot(aes(x = V1, y = V2, color = group)) +
  geom_point() +
  facet_wrap(~type, nrow = 1, scales = 'free') +
  theme_bw() +
  theme(axis.text = element_text(size = 7),
        legend.position = 'top',
        panel.grid = element_blank(),
        strip.background = element_blank()) +
  labs(x = 'Dim 1', y = 'Dim 2',
       color = 'Stage of Differentiation')) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/mds_binding_factors.png',
         height = 7, width = 24, units = 'cm')
