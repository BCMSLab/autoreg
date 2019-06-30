# loading required libraries
library(tidyverse)
library(org.Mm.eg.db)

tf_annotation <- select(org.Mm.eg.db,
                        keys = 'GO:0003700',
                        columns = 'SYMBOL',
                        keytype = 'GO')

write_rds(tf_annotation, 'data/tf_annotation.rds')