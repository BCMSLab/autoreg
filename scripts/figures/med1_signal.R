library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)

gene_id <- list(Foxo1 = 56458,
                Tfeb = 21425,
                Trp53 =  22059,
                Xbp1 =  22433)

track_signal <- read_rds('data/data_tracks.rds')

## MED1
png(filename = paste0('manuscript/figures/med1_signal.png'),
    height = 8, width = 24, units = 'cm', res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 4)))

imap(gene_id, function(x, .y) {
  # get promoter region
  prom_gr <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene,
                       filter = list(gene_id = x),
                       upstream = 3000,
                       downstream = 3000) %>%
    GenomicRanges::reduce()
  pushViewport(viewport(layout.pos.col = unname(which(gene_id == x)),
                        layout.pos.row = 1))
  # extract gr info
  gr_genome <- unique(as.character(genome(prom_gr)))
  gr_chrom <- unique(as.character(seqnames(prom_gr)))
  gr_strand <- unique(as.character(strand(prom_gr)))
  gr_start <- min(start(prom_gr))
  gr_end <- max(end(prom_gr))
  
  title_add = ifelse(.y == 'Foxo1', TRUE, FALSE)
  
  # MED1
  plotTracks(list(track_signal$MED1$`0`,
                  track_signal$MED1$`4`,
                  track_signal$MED1$`48`,
                  track_signal$MED1$`96`,
                  track_signal$MED1$`168`),
             type = 'h',
             chromosome = gr_chrom,
             from = gr_start,
             to = gr_end,
             add = TRUE,
             main = .y,
             showTitle = title_add,
             cex.main = 1)
  popViewport(1)
  return(NULL)
})
dev.off()
