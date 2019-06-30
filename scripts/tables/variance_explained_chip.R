library(tidyverse)
library(reshape2)
library(SummarizedExperiment)
library(DESeq2)
library(xtable)

binding_data <- read_rds('data/binding_data.rds')

mds <- binding_data %>%
  map(function(x) {
    dds <- DESeqDataSet(x, design = ~ group - 1)
    dds <- vst(dds)
    mds <- cmdscale(dist(t(assay(dds))),
                    eig = TRUE,
                    add = TRUE,
                    k = 2) 
    (mds$eig*100/sum(mds$eig))[1:2]
  })
mds %>%
  melt(.id = 'Factor') %>%
  setNames(c('value', 'Factor')) %>%
  filter(Factor != 'POLR2A') %>%
  mutate(value = as.integer(round(value)),
         dim = rep(paste('Dim', 1:2), 6)) %>%
  arrange(Factor) %>%
  spread(dim, value) %>%
  xtable(caption = '\\textbf{Percent of binding pattern variance explained by the stage of differentiation.}',
         label = 'tab:variance_explained_chip',
         align = 'clcc') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        sanitize.text.function = identity,
        comment = FALSE,
        file = 'manuscript/tables/variance_explained_chip.tex')
