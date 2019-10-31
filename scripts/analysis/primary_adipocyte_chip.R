# # PPARG
# GSE25836
# GSE59703 # also has RNA-seq
# GSE41578
# 
# # CEBPB hMSC + induction
# GSE68864

# GSE59703  25504365  GSM1443809  SRR1523496  PPARG hMADS Day 19
#                     GSM1443811  SRR1523498  PPARG hMADS Day 19  
#                     GSM1443817  SRR1523504  Input hMADS Day 19
#                     GSM1443818  SRR1523505  Input hMADS Day 19 
# GSE68864  26111340  GSM1684636	SRR2025078  CEBPB hMSC  Hour 6  
#                     GSM1684651	SRR2025101  Input hMSC  Hour 6

library(rtracklayer)
library(ChIPseeker)
library(tidyverse)

peaks <- c('data/CEBPB_hMSC_6h.bed.gz', 'data/PPARG_hMADS_19d.bed.gz') %>%
  map(function(x) {
    import.bed(x) %>% 
      annotatePeak(level = 'gene',
                   verbose = FALSE,
                   TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                   annoDb = 'org.Hs.eg.db') %>%
      as.GRanges()
  })

write_rds(peaks, 'data/primary_adipocyte_chip.rds')
