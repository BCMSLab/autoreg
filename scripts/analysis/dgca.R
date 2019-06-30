# loading required libraries
library(tidyverse)

ddcor <- read_rds('data/ddcor.rds')

dgca <- map(ddcor, function(x) {
  bind_rows(x, .id = 'group') %>%
    gather(group, cor, ends_with('cor')) %>%
    gather(group_pval, pval, ends_with('pval'))
}) %>%
  bind_rows() %>%
  na.omit() %>%
  unique()

write_rds(dgca, 'data/dgca.rds')
