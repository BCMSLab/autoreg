library(tidyverse)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(rtracklayer)

go_annotation <- read_rds('data/go_annotation.rds')
tf <- c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')
goi <- select(org.Mm.eg.db,
              c(unique(go_annotation$SYMBOL), tf, 'Cidec', 'Klf5'),
              'ENTREZID', 'SYMBOL')

prom_gr <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,  
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
  bind_rows(.id = 'gene_id') %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

bw <- list.files('data/bws_tissue', full.names = TRUE) %>%
  map(function(x) {
    import.bw(x,
              selection = BigWigSelection(prom_gr))
  })

names(bw) <- str_split(list.files('data/bws_tissue/'), '\\[|\\.bw', simplify = TRUE)[,2]

data_tracks <- bw %>%
  as.list() %>%
  map(as_tibble) %>%
  bind_rows(.id = 'sample') %>%
  mutate(source = sample) %>%
  separate(source, into = c('tissue', 'factor', 'rep')) %>%
  with(split(., .$sample))

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

data_tracks_tissue_tidy <- split(prom, prom$SYMBOL) %>%
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

write_rds(data_tracks_tissue_tidy, 'data/data_tracks_tissue_tidy.rds')
