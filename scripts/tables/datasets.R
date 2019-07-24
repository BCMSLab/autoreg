library(tidyverse)
library(SummarizedExperiment)
library(xtable)

peak_counts <- read_rds('data/peak_counts.rds')
gene_counts <- read_rds('data/gene_counts.rds')

tfs <- c('CEBPB', 'PPARG', 'POLR2A', 'RXRG', 'EP300', 'MED1')
hms <- c('H3K27ac', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K9me3')

pd <- colData(peak_counts)

t1 <- pd[, !colnames(pd) %in% 'qc'] %>%
  as_tibble() %>%
  filter(!is.na(factor),
         factor %in% c(tfs, hms)) %>%
  group_by(study) %>%
  summarise(n_sample = n(),
            factor = paste(unique(factor), collapse = '/'),
            time = paste(unique(time), collapse = '/'),
            group = paste(unique(group), collapse = '/'),
            pmid = as.character(unique(pmid)),
            bibtexkey = paste0('\\cite{', unique(bibtexkey), '}')) %>%
  arrange(desc(n_sample))

pd <- colData(gene_counts)

t2 <- pd[, !colnames(pd) %in% 'qc'] %>%
  as_tibble() %>%
  group_by(study) %>%
  summarise(n_sample = n(),
            factor = '',
            time = paste(unique(time), collapse = '/'),
            group = paste(unique(group), collapse = '/'),
            pmid = as.character(unique(pmid)),
            bibtexkey = paste0('\\cite{', unique(bibtexkey), '}')) %>%
  arrange(desc(n_sample))

bind_rows(list('ChIP Seq' = t1,
               'RNA Seq' = t2),
          .id = 'type') %>%
  mutate(type = ifelse(duplicated(type), '', type)) %>%
  setNames(c('Type', 'Study', '(N)', 'Factor', 'Time (hr)', 'Stage', 'PMID', 'Ref')) %>%
  xtable(caption = '\\textbf{Datasets of RNA and ChIP Seq.}',
         label = 'tab:datasets',
         align = 'cllclllll') %>%
  print(floating = TRUE,
        include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        sanitize.text.function = identity,
        comment = FALSE,
        add.to.row = list(pos = list(14), command = '\\midrule '),
        file = 'manuscript/tables/datasets.tex')
