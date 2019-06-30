library(Gviz)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
library(org.Mm.eg.db)

targets <- c('Atg4b', 'Ulk1', 'Map1lc3a', 'Map1lc3b', 'Sqstm1', 'Becn1')

symbol_entrez <- AnnotationDbi::select(org.Mm.eg.db,
                        targets,
                        'ENTREZID', 'SYMBOL')
gene_id <- symbol_entrez$ENTREZID
names(gene_id) <- symbol_entrez$SYMBOL

track_signal <- read_rds('data/data_tracks.rds')

## CEBPB
png(filename = paste0('manuscript/figures/cebpb_signal_direct.png'),
    height = 8, width = 30, units = 'cm', res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, length(targets))))

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
  
  title_add = ifelse(.y == 'Atg4b', TRUE, FALSE)
  
  # CEBPB
  plotTracks(list(track_signal$CEBPB$`0`,
                  track_signal$CEBPB$`2`,
                  track_signal$CEBPB$`4`,
                  track_signal$CEBPB$`6`,
                  track_signal$CEBPB$`48`),
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
