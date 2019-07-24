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
tf <- c('Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')

go_annotation <- read_rds('data/go_annotation.rds')
goi <- select(org.Mm.eg.db,
              c(unique(go_annotation$SYMBOL), tf),
              'ENTREZID', 'SYMBOL')

pd <- tibble(id = peak_counts$id,
             factor = peak_counts$factor) %>%
  filter(!is.na(factor),
         factor == 'POLR2A')

fls <- paste0('data/bws/', pd$id, '_EF.bw')
names(fls) <- pd$id

gene_gr <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene,  
                     filter = list(gene_id = goi$ENTREZID),
                     columns=c('gene_id',"tx_id")) %>%
  reduce()

polr2a_track_signal <- map(fls, function(x) {
  import.bw(x,
            selection = BigWigSelection(gene_gr))
})

write_rds(polr2a_track_signal, 'data/polr2a_track_signal.rds')
