library(tidyverse)
library(reshape2)
library(SummarizedExperiment)
library(xtable)

go_annotation <- read_rds('data/go_annotation.rds')
tf_annotation <- read_rds('data/tf_annotation.rds')

tf_autophagy <- intersect(go_annotation$SYMBOL, tf_annotation$SYMBOL)
tf <- c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')
targets <- c('Atg4b', 'Ulk1', 'Map1lc3a', 'Map1lc3b', 'Sqstm1', 'Becn1')

cat_gene <- list('Adipogenic TF' = tf,
     'Autophagy TF' = tf_autophagy,
     'Autophagy Gene' = targets) %>%
  melt() %>%
  setNames(c('row', 'cat')) %>%
  mutate(cat = factor(cat, levels = c('Adipogenic TF', 'Autophagy TF', 'Autophagy Gene')))

deg_res <- read_rds('data/deg_res.rds')

header <- paste0("\\multirow{2}{*}{Category} & \\multirow{2}{*}{Gene} &",
                 "\\multicolumn{2}{c}{Early vs Non} & \\multicolumn{2}{c}{Late vs Non} & \\multicolumn{2}{c}{Late vs Early} \\\\",
                 "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
                 "&& FC & SE & FC & SE & FC & SE\\\\")
deg_res %>%
  right_join(cat_gene) %>%
  filter(padj < .2) %>%
  mutate(log2FoldChange = round(log2FoldChange, 2),
         lfcSE = round(lfcSE, 2)) %>%
  unite(values, log2FoldChange, lfcSE, sep = '_') %>%
  dplyr::select(cat, row, contrast, values) %>%
  spread(contrast, values) %>%
  separate(early_vs_non, sep = '_', into = c('fc1', 'se1')) %>%
  separate(late_vs_non, sep = '_', into = c('fc2', 'se2')) %>%
  separate(late_vs_early, sep = '_', into = c('fc3', 'se3')) %>%
  mutate(cat = ifelse(duplicated(cat), '', as.character(cat))) %>%
  xtable(caption = '\\textbf{Significant differentially expressed genes of adipogenic and autophagy transcription factors.}',
         label = 'tab:deg_fc',
         align = 'cllcccccc') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        sanitize.text.function = identity,
        comment = FALSE,
        include.colnames=FALSE,
        add.to.row = list(pos = list(0, 4, 9),
                          command = c(header, rep('\\midrule ', 2))),
        file = 'manuscript/tables/deg_fc.tex')
