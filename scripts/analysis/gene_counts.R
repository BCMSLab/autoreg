# loading required libraries
library(curatedAdipoRNA)
library(SummarizedExperiment)
library(tidyverse)

# loading data
data("adipo_counts")

# subset by time points
time <- c(0, 2, 4, 24, 48, 144, 168, 192, 240)
se <- adipo_counts[, adipo_counts$time %in% time]

# subset by low counts
low_counts <- apply(assay(se), 1, function(x) length(x[x>10])>=2)
se <- se[low_counts, ]

# subset by low quality
qc <- unlist(se$qc, recursive = FALSE)
per_base <- lapply(qc, function(x) {
  df <- x[['per_base_sequence_quality']]
  df %>%
    dplyr::select(Base, Mean) %>%
    transform(Base = strsplit(as.character(Base), '-')) %>%
    unnest(Base) %>%
    mutate(Base = as.numeric(Base))
}) %>%
  bind_rows(.id = 'run') %>%
  group_by(run) %>%
  mutate(length = max(Base) > 150,
         run_average = mean(Mean) > 35)

per_base %>%
  ggplot(aes(x = Base, y = Mean, group = run, color = run_average)) +
  geom_line() +
  facet_wrap(~length, scales = 'free_x')

bad_samples <- data.frame(samples = unique(per_base$run[per_base$run_average == FALSE]))
bad_samples <- separate(bad_samples, col = samples, into = c('id', 'run'), sep = '\\.')

se <- se[, !se$id %in% bad_samples$id]

se$group <- group <- cut(se$time,
    breaks = c(-50, 0, 48, 240),
    labels = c('non', 'early', 'late'))

write_rds(se, 'data/gene_counts.rds')
