# loading required libraries
library(tidyverse)
library(DGCA)

transformed_counts <- read_rds('data/transformed_counts.rds')
ddcor <- read_rds('data/dgca.rds')

go_annotation <- read_rds('data/go_annotation.rds')
tf <- c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg')

set.seed(12345)
random_sample <- sample(rownames(transformed_counts), length(unique(ddcor$Gene1))-4)
goi <- c(random_sample, tf)

mat <- assay(transformed_counts)[rownames(transformed_counts) %in% goi,]

mod <- model.matrix(~factor(transformed_counts$stage)-1)
colnames(mod) <- paste('stage', 0:3, sep = '_')

ddcor <- map(tf, function(y) {
  map(c(stage_1 = 'stage_1', stage_3 = 'stage_3'),
      function(x) {
        ddcorAll(inputMat = mat,
                 design = mod,
                 compare = c(x, 'stage_0'),
                 adjust = 'fdr',
                 splitSet = y) %>%
          gather(stage, cor, ends_with('cor')) %>%
          gather(stage_pval, pval, ends_with('pval'))
      }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  na.omit() %>%
  unique()

write_rds(ddcor, 'data/random_ddcord.rds')
