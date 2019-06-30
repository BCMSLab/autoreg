library(tidyverse)
library(reshape2)
library(xtable)

markers <- list(Adipogenesis = c('Cebpa', 'Pparg'),
                Lipogenesis = c('Lpl', 'Acly', 'Dgat', 'Elov6', 'Fasn', 'Scd'),
                Autophagy = c('Map1lc3b', 'Sqstm1', 'Becn1'))

deg_res <- read_rds('data/deg_res.rds')

header <- paste0("\\multirow{2}{*}{Category} & \\multirow{2}{*}{Gene} &",
                 "\\multicolumn{3}{c}{Early vs Non} & \\multicolumn{3}{c}{Late vs Non} & \\multicolumn{3}{c}{Late vs Early} \\\\",
                 "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}",
                 "&& FC & SE & FDR & FC & SE & FDR & FC & SE & FDR\\\\")
melt(markers) %>%
  setNames(c('row', 'cat')) %>%
  inner_join(deg_res) %>%
  mutate(log2FoldChange = round(log2FoldChange, 2),
         lfcSE = round(lfcSE, 2),
         padj = ifelse(padj < .2, '$< 0.2$', as.character(round(padj, 2)))) %>%
  unite(values, log2FoldChange, lfcSE, padj, sep = '_') %>%
  select(cat, row, contrast, values) %>%
  spread(contrast, values) %>%
  separate(early_vs_non, sep = '_', into = c('fc1', 'se1', 'padj1')) %>%
  separate(late_vs_non, sep = '_', into = c('fc2', 'se2', 'padj2')) %>%
  separate(late_vs_early, sep = '_', into = c('fc3', 'se3', 'padj3')) %>%
  mutate(cat = ifelse(duplicated(cat), '', cat)) %>%
  xtable(caption = '\\textbf{Differentially expression of gene markers.}',
         label = 'tab:deg_markers',
         align = 'cllccccccccc') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        sanitize.text.function = identity,
        comment = FALSE,
        include.colnames=FALSE,
        add.to.row = list(pos = list(0, 2, 5),
                          command = c(header, rep('\\midrule ', 2))),
        file = 'manuscript/tables/deg_markers.tex')

