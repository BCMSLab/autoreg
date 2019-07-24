# loading requried libraries
library(tidyverse)
library(SummarizedExperiment)

# loading data
peak_counts <- read_rds('data/peak_counts.rds')

# subset the object
tfs <- c('CEBPB', 'PPARG', 'POLR2A', 'RXRG', 'EP300', 'MED1')

se <- map(tfs, function(x) {
  ind <- is.na(peak_counts$factor)
  e <- peak_counts[, !ind]
  sample_ind <- (e$factor == x)
  sample_ids <- colnames(e)[sample_ind]
  
  # select peaks from the selected samples
  peak_ind <- lapply(mcols(e)$peak, function(x) sum(sample_ids %in% x))
  peak_ind <- unlist(peak_ind) > 2
  
  # subset the object
  se <- e[peak_ind, sample_ind]
  low_counts <- apply(assay(se), 1, function(x) length(x[x>5])>=2)
  se <- se[low_counts, ]
  se
})
names(se) <- tfs


write_rds(se, 'data/binding_data.rds')
