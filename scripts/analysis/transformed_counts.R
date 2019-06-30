# loading required libraries
library(DESeq2)
library(SummarizedExperiment)
library(tidyverse)

# loading data
gene_counts <- read_rds('data/gene_counts.rds')

# normalization
dds <- DESeqDataSet(gene_counts, 
                    design = ~ group)

# transformation
dds <- DESeq(dds)
transformed_counts <- vst(dds)

write_rds(transformed_counts, 'data/transformed_counts.rds')
