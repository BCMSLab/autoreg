# loading required libraries
library(tidyverse)
library(SummarizedExperiment)
library(org.Mm.eg.db)

# loading data
untar('data/galaxy_cebpb_kd/70_ featureCounts on collection 61_ Counts.tar',
      exdir = 'data/galaxy_cebpb_kd/')

count_fls <- list.files('data/galaxy_cebpb_kd/featureCounts on collection 61: Counts',
           full.names = TRUE) %>%
  map(read_tsv) %>%
  map(column_to_rownames, var = 'Geneid')

mat <- do.call('cbind', count_fls) %>% as.matrix()

pd <- read_tsv('data/galaxy_cebpb_kd/SraRunTable.txt') %>%
  dplyr::select(id = Sample_Name, srr = Run, time, group = treatment) %>%
  mutate(time = substr(time, start = 1, stop = 1),
         group = str_split(group, ' ', simplify = TRUE)[, 1],
         group = ifelse(group == 'Lentivirus', 'Knockdown', 'Control')) %>%
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

write_rds(se, 'data/counts_cebpb_kd.rds')
