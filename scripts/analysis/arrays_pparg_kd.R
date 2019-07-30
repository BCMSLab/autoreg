library(GEOquery)
library(WGCNA)
library(tidyverse)

allowWGCNAThreads(4)

e <- getGEO('GSE12929',
            destdir = 'data/')[[1]]

e$time <- str_split(e$title, '_|rep', simplify = TRUE)[,1]
e$group <- str_split(e$title, '_|rep', simplify = TRUE)[,3]
e <- e[, !is.na(e$time)]
e_sub <- e

e_sub <- e_sub[rowSums(exprs(e_sub)) >= 0,]
exprs(e_sub)[exprs(e_sub) < 0] <- 0

pd <- pData(e_sub)[, c('geo_accession', 'time', 'group')]
fd <- fData(e_sub)
mat <- exprs(e_sub)

mat <- collapseRows(mat, 
                    rowGroup = fd$Symbol,
                    rowID = rownames(mat))[[1]]

eset <- ExpressionSet(mat,
                      phenoData = new('AnnotatedDataFrame', pd))

eset$time <- str_sub(eset$time, start = 4, end = 4) %>% as.numeric()
eset$group <- ifelse(eset$group == 'neg', 'Control', 'Knockdown')

write_rds(eset, 'data/arrays_pparg_kd.rds')
