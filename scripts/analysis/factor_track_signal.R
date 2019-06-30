library(tidyverse)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicRanges)
library(org.Mm.eg.db)

# if(!file.exists('data/bws.tar.gz') {
#   download.file('',
#                 destfile = 'data/bws.tar.gz')
#   untar('data/bws.tar.gz')
# })
  
peak_counts <- read_rds('data/peak_counts.rds')
tf <- c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')

go_annotation <- read_rds('data/go_annotation.rds')
goi <- AnnotationDbi::select(org.Mm.eg.db,
              c(unique(go_annotation$SYMBOL), tf, 'Cidec', 'Klf5'),
              'ENTREZID', 'SYMBOL')

pd <- tibble(id = peak_counts$id,
             factor = peak_counts$factor) %>%
  filter(!is.na(factor),
         factor %in% toupper(tf))

fls <- paste0('data/bws/', pd$id, '_EF.bw')
names(fls) <- pd$id

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
  
factor_track_signal <- map(fls, function(x) {
  import.bw(x,
            selection = BigWigSelection(prom_gr))
})

write_rds(factor_track_signal, 'data/factor_track_signal.rds')
