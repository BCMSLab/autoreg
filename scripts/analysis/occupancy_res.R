library(tidyverse)
library(reshape2)
library(SummarizedExperiment)
library(DESeq2)

occupancy_res <- list()

occupancy <- read_rds('data/factor_occupancy.rds')

occupancy_res$tf <- map(unique(occupancy$factor), function(x) {
  df <- occupancy %>%
    filter(factor == x)
  mat <-  df%>%
    acast(geneId ~ id, value.var = 'count', fill = 0)
  
  pd <- df %>%
    dplyr::select(id, time, group) %>%
    unique()
  
  dds <- DESeqDataSetFromMatrix(mat, pd, design = ~ group)
  dds <- DESeq(dds)
  map(c(early = 'early', late = 'late'), {
    function(y) results(dds,
                        contrast = c('group', y, 'non'),
                        tidy = TRUE,
                        pAdjustMethod = 'fdr')
  }) %>%
    bind_rows(.id = 'contrast')
}) %>%
  set_names(unique(occupancy$factor)) %>%
  bind_rows(.id = 'factor')


occupancy_hm <- read_rds('data/hm_occupancy.rds') %>%
  filter(factor %in% c('H3K27ac', 'H3K4me1', 'H3K4me3'))

occupancy_res$hm <- map(unique(occupancy_hm$factor), function(x) {
  df <- filter(occupancy_hm, factor == x)
  mat <- acast(df, geneId ~ id, value.var = 'count', fill = 0)
  
  pd <- df %>%
    dplyr::select(id, time, group) %>%
    unique()
  
  dds <- DESeqDataSetFromMatrix(mat, pd, design = ~ group)
  dds <- DESeq(dds)
  map(c(early = 'early', late = 'late'), {
    function(y) results(dds,
                        contrast = c('group', y, 'non'),
                        tidy = TRUE,
                        pAdjustMethod = 'fdr')
  }) %>%
    bind_rows(.id = 'contrast')
}) %>%
  set_names(unique(occupancy_hm$factor)) %>%
  bind_rows(.id = 'factor')

write_rds(occupancy_res, 'data/occupancy_res.rds')
