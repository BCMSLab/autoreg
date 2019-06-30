# loading required libraries
library(DESeq2)
library(SummarizedExperiment)
library(tidyverse)

# loading data
gene_counts <- read_rds('data/gene_counts.rds')

# normalization
dds <- DESeqDataSet(gene_counts, 
                    design = ~ group-1)

# transformation
dds <- DESeq(dds)

res <- list(early_vs_non = c('group', 'early', 'non'),
     late_vs_non = c('group', 'late', 'non'),
     late_vs_early = c('group', 'late', 'early')) %>%
  map(function(x) {
    results(dds,
            contrast = x,
            tidy = TRUE,
            pAdjustMethod = 'fdr')
  }) %>%
  bind_rows(.id = 'contrast') %>%
  mutate(contrast = factor(contrast, levels = c('early_vs_non', 'late_vs_non', 'late_vs_early')))

write_rds(res, 'data/deg_res.rds')
