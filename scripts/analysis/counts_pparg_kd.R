# loading required libraries
library(tidyverse)
library(SummarizedExperiment)
library(org.Mm.eg.db)

# loading data
untar('data/galaxy_pparg_kd/24_ featureCounts on collection 17_ Counts.tar',
      exdir = 'data/galaxy_pparg_kd/')

count_fls <- list.files('data/galaxy_pparg_kd/featureCounts on collection 52: Counts/',
                        full.names = TRUE) %>%
  map(read_tsv) %>%
  map(column_to_rownames, var = 'Geneid')

mat <- do.call('cbind', count_fls) %>% as.matrix()

pd <- read_tsv('data/galaxy_pparg_kd/SraRunTable.txt') %>%
  dplyr::select(id = Sample_Name, srr = Run,
                time = stage_of_adipogenesis, group = source_name) %>%
  mutate(time = rep(c('0h', '1h', '2h', '4h', '2d', '7d'), 2),
         time = factor(time, levels = c('0h', '1h', '2h', '4h', '2d', '7d')),
         group = rep(c('white', 'brown'), each = 6)) %>%
  as.data.frame()

rownames(pd) <- pd$id
colnames(mat) <- pd$id

fd <- AnnotationDbi::select(org.Mm.eg.db,
                            rownames(mat),
                            'SYMBOL',
                            'ENTREZID')

ind <- duplicated(fd$SYMBOL) | is.na(fd$SYMBOL)

mat <- mat[!ind,]
rownames(mat) <- fd$SYMBOL[!ind]

se <- SummarizedExperiment(mat, colData = pd)

write_rds(se, 'data/counts_pparg_kd.rds')
