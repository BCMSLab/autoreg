library(tidyverse)
library(stringapi)

go_annotation <- read_rds('data/go_annotation.rds')
tf_annotation <- read_rds('data/tf_annotation.rds')

tf_autophagy <- intersect(go_annotation$SYMBOL, tf_annotation$SYMBOL)
tf <- c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')

ids <- c(tf_autophagy, tf)

interactions <- list()
interactions$tf_tf <- network(identifiers = ids,
                              species = 10090)

targets <- c('Atg4b', 'Ulk1', 'Map1lc3a', 'Map1lc3b', 'Sqstm1', 'Becn1')

ids <- c(tf, targets)

interactions$tf_gene <- network(identifiers = ids,
                              species = 10090)

write_rds(interactions, 'data/interactions.rds')
