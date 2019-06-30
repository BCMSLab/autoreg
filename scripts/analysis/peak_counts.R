# loading requried libraries
library(tidyverse)
library(SummarizedExperiment)

# load data from bioc
# loading data
load('data/peak_counts.rda')

peak_counts$group <- cut(peak_counts$time,
                         breaks = c(-50, 0, 48, 240),
                         labels = c('non', 'early', 'late'))

peak_counts$group <- ifelse(is.na(peak_counts$group) & peak_counts$stage == 0,
                            'non',
                            as.character(peak_counts$group))

peak_counts$group <- ifelse(is.na(peak_counts$group) & peak_counts$stage == 3,
                            'late',
                            as.character(peak_counts$group))

write_rds(peak_counts, 'data/peak_counts.rds')