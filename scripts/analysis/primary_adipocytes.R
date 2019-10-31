# load required libraries
library(GEOquery)
library(tidyverse)
library(WGCNA)
library(sva)

allowWGCNAThreads(4)

if (!file.exists('data/GSE98680_series_matrix.txt.gz')) {
  getGEO('GSE98680',
         destdir = 'data/')
}

eset <- getGEO('GSE98680', destdir = 'data/')[[1]][, 1:24]

mat <- collapseRows(exprs(eset),
                    rowID = featureNames(eset),
                    rowGroup = fData(eset)$GENE_SYMBOL)[[1]]

newset <- ExpressionSet(mat)

newset$group <- str_split(eset$title, '_', simplify = TRUE)[, 3]
newset$patient <- str_split(eset$title, '_', simplify = TRUE)[, 2]

exprs(newset) <- ComBat(exprs(newset),
                        batch = newset$patient,
                        mod = model.matrix(~group, data = pData(newset)))

write_rds(newset, 'data/primary_adipocyte.rds')

if (!file.exists('data/GSE94752_series_matrix.txt.gz')) {
  getGEO('GSE94752',
         destdir = 'data/')
}

eset <- getGEO('GSE94752', destdir = 'data/')[[1]]

fData(eset)$symbol <- str_split(fData(eset)$gene_assignment,
                                ' // ',
                                simplify = TRUE)[, 2]
eset <- eset[fData(eset)$symbol != '',]

mat <- collapseRows(exprs(eset),
                    rowID = featureNames(eset),
                    rowGroup = fData(eset)$symbol)[[1]]


newset <- ExpressionSet(mat)

newset$group <- str_split(eset$title, ' ', simplify = TRUE)[, 1]
#newset$batch <- str_split(eset$title, '_|\\]', simplify = TRUE)[, 3]
#newset$batch <- str_sub(eset$batch, 1, 1)

write_rds(newset, 'data/primary_adipocyte_insulin.rds')
