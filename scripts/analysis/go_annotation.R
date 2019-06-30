# loading required libraries
library(tidyverse)
library(org.Mm.eg.db)

go_annotation <- select(org.Mm.eg.db,
                        keys = 'GO:0006914',
                        columns = 'SYMBOL',
                        keytype = 'GO')

write_rds(go_annotation, 'data/go_annotation.rds')

tf_annotation <- select(org.Mm.eg.db,
                        keys = 'GO:0003700',
                        columns = 'SYMBOL',
                        keytype = 'GO')

write_rds(tf_annotation, 'data/tf_annotation.rds')
 