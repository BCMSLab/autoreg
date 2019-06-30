library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)

peak_counts <- read_rds('data/peak_counts.rds')

pd <- tibble(factor = peak_counts$factor,
             sample = peak_counts$id,
             group = peak_counts$group,
             time = peak_counts$time)

data_tracks <- read_rds('data/factor_track_signal.rds') %>%
  as.list() %>%
  map(as_tibble) %>%
  bind_rows(.id = 'sample') %>%
  inner_join(pd) %>%
  with(split(., .$factor)) %>%
  map(function(x) {
    makeGRangesListFromDataFrame(x,
                                 split.field = 'time',
                                 keep.extra.columns = TRUE) %>%
      as.list() %>%
      imap(function(f, .y) {
        DataTrack(f,
                  name = .y,
                  background.title = 'white',
                  col.axis = 'black',
                  col.title = 'black')
      })
  })

write_rds(data_tracks, 'data/data_tracks.rds')
