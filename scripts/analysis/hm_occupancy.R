library(tidyverse)
library(reshape2)
library(SummarizedExperiment)

histone_modification <- read_rds('data/histone_modification.rds')

hm_occupancy <- map(histone_modification, function(x) {
  se <- x
  
  pd <- colData(se)[, c('id', 'time', 'group')] %>% as.data.frame()
  fd <- mcols(se)[, c('name', 'geneId', 'annotation')]%>% as.data.frame()
  
  melt(assay(se)) %>%
    setNames(c('name', 'id', 'count')) %>%
    left_join(pd) %>%
    left_join(fd) %>%
    mutate(annotation = str_split(annotation, ' \\(', simplify = TRUE)[, 1]) %>%
    group_by(geneId, id, time, group) %>%
    summarise(count = sum(count)) %>%
    ungroup()
}) %>%
  bind_rows(.id = 'factor')

write_rds(hm_occupancy, 'data/hm_occupancy.rds')
