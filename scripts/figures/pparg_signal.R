library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)

gene_id <- list(Foxo1 = 56458,
                Tfeb = 21425,
                Trp53 =  22059,
                Xbp1 =  22433,
                Cidec = 14311)

track_signal <- read_rds('data/data_tracks.rds')

# generate figures
## PPARG
png(filename = paste0('manuscript/figures/pparg_signal.png'),
    height = 10, width = 24, units = 'cm', res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 5)))

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
  
  # PPARG
  plotTracks(list(track_signal$PPARG$`0`,
                  track_signal$PPARG$`24`,
                  track_signal$PPARG$`48`,
                  track_signal$PPARG$`72`,
                  track_signal$PPARG$`144`,
                  track_signal$PPARG$`168`),
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
