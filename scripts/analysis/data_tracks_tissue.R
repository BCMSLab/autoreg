library(tidyverse)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(rtracklayer)
library(Gviz)

go_annotation <- read_rds('data/go_annotation.rds')
tf <- c('Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')
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

nms <- str_split(names(bw), '_', simplify = TRUE)[, 1]

data_tracks_tissue <- list()


for(i in 1:length(bw)) {
  onms <- names(bw)[i]
  x <- bw[[i]]
  y <- nms[i]

  data_tracks_tissue[onms] <- DataTrack(x,
                        name = y,
                        background.title = 'white',
                        col.axis = 'black',
                        col.title = 'black')
}

write_rds(data_tracks_tissue, 'data/data_tracks_tissue.rds')
