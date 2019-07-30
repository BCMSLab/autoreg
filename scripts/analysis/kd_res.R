library(tidyverse)
library(DESeq2)
library(limma)

se <- read_rds('data/counts_cebpb_kd.rds')
se <- se[rowSums(assay(se)) > 0,]

targets <- list('Adipogenic TF' = c('Cebpa', 'Cebpb', 'Pparg'),
                'Lipogenesis' = c('Lpl', 'Acly', 'Fasn'),
                'Autophagy TF' = c('Foxo1', 'Tfeb', 'Xbp1'),
                'Autophagy Gene' = c('Map1lc3b', 'Becn1', 'Sqstm1')) %>%
  reshape2::melt() %>%
  setNames(c('gene', 'category'))

sum(rownames(se) %in% unlist(targets$gene))

cebpb_kd_res <- map(c('0' = 0, '4' = 4), function(x) {
  e <- se[,se$time == x]
  dds <- DESeqDataSet(e, design = ~group-1)
  dds <- DESeq(dds)
  results(dds, tidy = TRUE, pAdjustMethod = 'fdr', cooksCutoff=FALSE) %>%
    as_tibble()
}) %>%
  bind_rows(.id = 'time') 

write_rds(cebpb_kd_res, 'data/cebpb_kd_res.rds')

# pparg
eset <- read_rds('data/arrays_pparg_kd.rds')
eset <- eset[, !is.na(eset$time)]
exprs(eset) <- normalizeBetweenArrays(exprs(eset))

exprs(eset) <- log2(exprs(eset) + 1)

pparg_kd_res <- map(c('0' = 0, '2' = 2, '4' = 4, '5' = 5, '6' = 6), function(x) {
  e <- eset[, eset$time == x]
  mat <- exprs(e)
  mod <- model.matrix(~group, data = pData(e))
  fit <- lmFit(mat, design = mod, method = 'robust')
  fit <- eBayes(fit)
  
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled
  
  topTable(fit,
           sort.by = 'none',
           number = Inf,
           adjust.method = 'fdr') %>%
    rownames_to_column('gene_id') %>%
    mutate(se = se[, 2]) %>%
    as_tibble()
}) %>%
  bind_rows(.id = 'time')

write_rds(pparg_kd_res, 'data/pparg_kd_res.rds')
