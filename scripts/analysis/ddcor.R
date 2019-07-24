# loading required libraries
library(tidyverse)
library(SummarizedExperiment)
library(DGCA)

transformed_counts <- read_rds('data/transformed_counts.rds')

go_annotation <- read_rds('data/go_annotation.rds')
tf <- c('Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')

goi <- c(unique(go_annotation$SYMBOL), tf)

mat <- assay(transformed_counts)[rownames(transformed_counts) %in% goi,]

mod <- model.matrix(~factor(transformed_counts$group)-1)
colnames(mod) <- levels(transformed_counts$group)

ddcor <- map(tf, function(y) {
  map(list(early = 'early', late = 'late'), function(x) {
    ddcorAll(inputMat = mat,
             design = mod,
             compare = c(x, 'non'),
             splitSet = y)
  })
})

write_rds(ddcor, 'data/ddcor.rds')
