# loading requried libraries
library(tidyverse)

# loading data
fls <- list.files('data/beta_output/direct_targets', full.names = TRUE)
file_md <- str_split(fls, '_|\\(|\\)', simplify = TRUE)[, c(8, 12, 13)] %>%
  as_tibble() %>%
  setNames(c('dir', 'factor', 'time')) %>%
  mutate(file = fls)

names(fls) <- fls

direct_targets <- map(fls, read_tsv) %>%
  bind_rows(.id = 'file') %>%
  left_join(file_md)

# loading data
fls <- list.files('data/beta_output/associated_peaks', full.names = TRUE)
file_md <- str_split(fls, '_|\\(|\\)|\\[', simplify = TRUE)[, c(4, 8, 9)] %>%
  as_tibble() %>%
  setNames(c('dir', 'factor', 'time')) %>%
  mutate(file = fls)

names(fls) <- fls

associated_peaks <- map(fls, read_tsv) %>%
  bind_rows(.id = 'file') %>%
  left_join(file_md)
