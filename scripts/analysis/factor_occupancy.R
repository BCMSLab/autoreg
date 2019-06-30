library(tidyverse)
library(reshape2)
library(SummarizedExperiment)

binding_data <- read_rds('data/binding_data.rds')

occupancy <- map(binding_data, function(x) {
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


write_rds(occupancy, 'data/factor_occupancy.rds')

# occupancy_region <- map(binding_data, function(x) {
#   se <- x
#   
#   pd <- colData(se)[, c('id', 'time', 'group')] %>% as.data.frame()
#   fd <- mcols(se)[, c('name', 'geneId', 'annotation')]%>% as.data.frame()
#   
#   melt(assay(se)) %>%
#     setNames(c('name', 'id', 'count')) %>%
#     left_join(pd) %>%
#     left_join(fd) %>%
#     mutate(annotation = str_split(annotation, ' \\(', simplify = TRUE)[, 1]) %>%
#     group_by(geneId, id, time, group, annotation) %>%
#     summarise(count = sum(count)) %>%
#     ungroup()
# }) %>%
#   bind_rows(.id = 'factor')