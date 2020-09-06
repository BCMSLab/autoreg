# loading requried libraries
library(tidyverse)
library(SummarizedExperiment)
library(ExperimentHub)

# load data from bioc
# loading data
# load('data/peak_counts2.rda')
eh <- ExperimentHub()
peak_counts2 <- query(eh, "curatedAdipoChIP")[[1]]

peak_counts2$group <- cut(peak_counts2$time,
                          breaks = c(-50, 0, 48, 240),
                          labels = c('non', 'early', 'late'))

peak_counts2$group <- ifelse(is.na(peak_counts2$group) & peak_counts2$stage == 0,
                            'non',
                            as.character(peak_counts2$group))

peak_counts2$group <- ifelse(is.na(peak_counts2$group) & peak_counts2$stage == 3,
                            'late',
                            as.character(peak_counts2$group))

write_rds(peak_counts2, 'data/peak_counts.rds')