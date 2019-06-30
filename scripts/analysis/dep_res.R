library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)

binding_data <- read_rds('data/binding_data.rds')

ddss <- binding_data %>%
  map(function(x) {
    dds <- DESeqDataSet(x, design = ~ group - 1)
    dds <- DESeq(dds)
    dds
  })

res1 <- list(
  early_vs_non = map(ddss, function(x) {
    results(x,
            contrast = c('group', 'early', 'non'),
            tidy = TRUE,
            pAdjustMethod = 'fdr')
  }) %>%
    bind_rows(.id = 'factor'),
  late_vs_non = map(ddss, function(x) {
    results(x,
            contrast = c('group', 'late', 'non'),
            tidy = TRUE,
            pAdjustMethod = 'fdr')
  })  %>%
    bind_rows(.id = 'factor'),
  late_vs_early = map(ddss, function(x) {
    results(x,
            contrast = c('group', 'late', 'early'),
            tidy = TRUE,
            pAdjustMethod = 'fdr')
  })  %>%
    bind_rows(.id = 'factor')
) %>%
  bind_rows(.id = 'contrast')

histone_modification <- read_rds('data/histone_modification.rds')

ddss <- histone_modification %>%
  map(function(x) {
    dds <- DESeqDataSet(x, design = ~ group - 1)
    dds <- DESeq(dds)
    dds
  })

res2 <- list(
  early_vs_non = map(ddss[-5], function(x) {
    results(x,
            contrast = c('group', 'early', 'non'),
            tidy = TRUE,
            pAdjustMethod = 'fdr')
  }) %>%
    bind_rows(.id = 'factor'),
  late_vs_non = map(ddss[-3], function(x) {
    results(x,
            contrast = c('group', 'late', 'non'),
            tidy = TRUE,
            pAdjustMethod = 'fdr')
  })  %>%
    bind_rows(.id = 'factor'),
  late_vs_early = map(ddss[c(-3, -5)], function(x) {
    results(x,
            contrast = c('group', 'late', 'early'),
            tidy = TRUE,
            pAdjustMethod = 'fdr')
  })  %>%
    bind_rows(.id = 'factor')
) %>%
  bind_rows(.id = 'contrast')


write_rds(rbind(res1, res2), 'data/dep_res.rds')
