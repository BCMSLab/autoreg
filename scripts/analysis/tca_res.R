library(tidyverse)
library(TCseq)
library(SummarizedExperiment)

go_annotation <- read_rds('data/go_annotation.rds')

transformed_counts <- read_rds('data/transformed_counts.rds')
se <- transformed_counts[, transformed_counts$time %in% c(0, 4, 24, 48, 168, 192)]

gene_counts <- read_rds('data/gene_counts.rds')
se <- gene_counts[rownames(gene_counts) %in% c(unique(go_annotation$SYMBOL), 'Cebpb', 'Pparg')]

table(se$time)
design <- data.frame(
  sampleid = se$id,
  timepoint = as.factor(se$time),
  group = se$stage,
  stringsAsFactors = FALSE)

gf <- as.data.frame(rowRanges(se), row.names = NULL)[, c(1:3, 6)]
names(gf) <- c('chr', 'start', 'end', 'id')

counts <- round(assay(se))

tca <- TCA(design = design,
           counts = counts,
           genomicFeature = gf)

# run analysis
tca <- DBanalysis(tca,
                  filter.type = "raw",
                  filter.value = 10,
                  samplePassfilter = 2)

tca_res <- DBresult(tca,
                    contrasts = colnames(tca@contrasts))

tca_res %>%
  as.list() %>%
  map(function(x) mcols(x) %>% as.data.frame) %>%
  bind_rows(.id = 'contrast') %>%
  write_rds('data/tca_res.rds')

tca <- timecourseTable(tca,
                       value = 'expression',
                       norm.method = 'rpkm',
                       filter = TRUE)
tca <- timeclust(tca, 'cm', 5)

write_rds(tca, 'data/tca.rds')
