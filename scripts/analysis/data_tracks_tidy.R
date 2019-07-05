library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(org.Mm.eg.db)

peak_counts <- read_rds('data/peak_counts.rds')
factor_track_signal <- read_rds('data/factor_track_signal.rds')
go_annotation <- read_rds('data/go_annotation.rds')

pd <- tibble(factor = peak_counts$factor,
             sample = peak_counts$id,
             group = peak_counts$group,
             time = peak_counts$time)

tf <- c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')


goi <- AnnotationDbi::select(org.Mm.eg.db,
                             c(unique(go_annotation$SYMBOL), tf, 'Cidec', 'Klf5'),
                             'ENTREZID', 'SYMBOL')

prom <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,  
                     filter = list(gene_id = goi$ENTREZID),
                     upstream = 3000,
                     downstream = 3000,
                     use.names = FALSE,
                     columns=c('gene_id',"tx_id")) %>%
  with(split(., as.integer(.$gene_id))) %>%
  as.list() %>%
  map(function(x) {
    tibble(seqnames = unique(as.character(seqnames(x))),
           start = min(start(x)),
           end = max(end(x)),
           strand = unique(as.character(strand(x))))
  }) %>%
  bind_rows(.id = 'ENTREZID') %>%
  left_join(goi) %>%
  dplyr::select(-ENTREZID, -strand)

data_tracks <- factor_track_signal %>%
  as.list() %>%
  map(as_tibble) %>%
  bind_rows(.id = 'sample') %>%
  inner_join(pd) %>%
  with(split(., .$factor)) 


data_tracks_tidy <- split(prom, prom$SYMBOL) %>%
  map(function(x) {
    map(data_tracks, function(f) {
      f %>%
        filter(seqnames == x$seqnames,
               start >= x$start,
               end <= x$end) %>%
        mutate(pos = (start + end)/2,
               pos = pos - min(pos)-3000) %>%
        na.omit()
    }) %>%
      bind_rows(.id = 'factor')
    })%>%
  bind_rows(.id = 'gene_id')

write_rds(data_tracks_tidy, 'data/data_tracks_tidy.rds')
