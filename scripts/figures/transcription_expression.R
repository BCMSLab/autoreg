library(tidyverse)
library(reshape2)
library(GGally)
library(SummarizedExperiment)

go_annotation <- read_rds('data/go_annotation.rds')
occupancy <- read_rds('data/factor_occupancy.rds')
transformed_counts <- read_rds('data/transformed_counts.rds')


df1 <- occupancy %>%
  filter(factor == 'POLR2A',
         geneId %in% go_annotation$SYMBOL) %>%
  group_by(geneId, time, group) %>%
  summarise(count = log2(mean(count)+1)) %>%
  ungroup()

df2 <- assay(transformed_counts)[rownames(transformed_counts) %in% go_annotation$SYMBOL,] %>%
  melt() %>%
  setNames(c('geneId', 'id', 'count')) %>%
  left_join(tibble(id = transformed_counts$id,
                   time = transformed_counts$time,
                   group = transformed_counts$group)) %>%
  group_by(geneId, time, group) %>%
  summarise(count = mean(count)) %>%
  ungroup()

df <- list(polr2a = df1, rna = df2) %>%
  bind_rows(.id = 'type') %>%
  filter(time %in% intersect(df1$time, df2$time)) %>%
  dplyr::select(-group) %>%
  unique() %>%
  unite(col, type, time) %>%
  spread(col, count) %>%
  na.omit() %>%
  dplyr::select(-geneId)

xcol <- grep('polr2a_', names(df))
nxcol <- str_split(names(df)[xcol], '_', simplify = TRUE)[, 2] %>% as.numeric()
ind <- order(nxcol)

xrow <- grep('rna_', names(df))
nxrow <- str_split(names(df)[xrow], '_', simplify = TRUE)[, 2] %>% as.numeric()
ind2 <- order(nxrow)

(ggduo(df, 
      columnsX = xcol[ind],
      columnLabelsX = nxcol[ind],
      columnsY = xrow[ind2],
      columnLabelsY = nxrow[ind2],
            types = list(continuous = wrap('smooth_lm',
                                           color = 'gray',
                                           se = FALSE,
                                           alpha = .7))) +
  labs(x = 'Transcriptional Activity (POLR2A)',
       y = 'Gene Expression (mRNA)') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0,"null"))) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/transcription_expression.png',
         width = 20, height = 20, units = 'cm')
