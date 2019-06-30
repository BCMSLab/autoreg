library(tidyverse)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# if(!file.exists('data/peaks.tar.gz') {
#   download.file('',
#                 destfile = 'data/peaks.tar.gz')
#   untar('data/peaks.tar.gz')
# })

hms <- c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9me3')
tfs <- c('CTCF', 'CEBPB', 'PPARG', 'POLR2A', 'RXRG', 'EP300', 'MED1')

peak_counts <- read_rds('data/peak_counts.rds')

pd <- tibble(
  factor = peak_counts$factor,
  group = peak_counts$group,
  time = peak_counts$time,
  id = peak_counts$id
) %>%
  filter(!is.na(factor)) %>%
  filter(factor %in% c(tfs, hms)) %>%
  mutate(file = paste0('data/peaks/', id, '_peaks.xls'))

peak_overlaps <- map(pd$file, function(x) {
  enrichPeakOverlap(queryPeak = x,
                    targetPeak = pd$file,
                    TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                    pAdjustMethod = 'fdr',
                    nShuffle = 50,
                    mc.cores = 3)
}) %>%
  bind_rows() %>%
  mutate_at(vars(1:2), function(x) str_split(x, '_', simplify = TRUE)[, 1])

write_rds(peak_overlaps, 'data/peak_overlaps.rds')
