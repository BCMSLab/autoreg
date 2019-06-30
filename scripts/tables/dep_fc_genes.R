library(tidyverse)
library(xtable)
library(SummarizedExperiment)

binding_data <- read_rds('data/binding_data.rds')
dep_res <- read_rds('data/dep_res.rds')

ind <- c('Atg4b', 'Ulk1', 'Map1lc3a', 'Map1lc3b', 'Sqstm1', 'Becn1')

peak_symbol <- map(binding_data, function(x){
  tibble(row = mcols(x)$name,
         gene_id = mcols(x)$geneId)
}) %>%
  bind_rows() %>%
  unique() %>%
  filter(gene_id %in% ind)

header <- paste0("\\multirow{2}[3]{*}{Category} & \\multirow{2}[3]{*}{Factor} & \\multirow{2}[3]{*}{Gene} &",
                 " \\multicolumn{4}{c}{Early vs Non} & \\multicolumn{4}{c}{Late vs Non} &  \\multicolumn{4}{c}{Late vs Early} \\\\",
                 " \\cmidrule(lr){4-7} \\cmidrule(lr){8-11} \\cmidrule(lr){12-15}",
                 "&&& (N) & Range & Ave & SD & (N) & Range & Ave & SD & (N) & Range & Ave & SD \\\\")

cat_fac <- list(factor = c('CTCF', 'CEBPB', 'PPARG'),
                co_factor = c('EP300', 'MED1', 'RXRG'),
                hm = c('H3K27ac', 'H3K4me3'))
fac <- tibble(cat = factor(rep(c('Factor', 'Cofactor', 'Histone Marker'), times = c(3,3,2)),
                           levels = c('Factor', 'Cofactor', 'Histone Marker')),
              factor = factor(unlist(cat_fac, use.names = FALSE),
                              levels = unlist(cat_fac, use.names = FALSE)))

peak_symbol %>%
  inner_join(dep_res) %>%
  filter(padj < .2) %>%
  group_by(factor, contrast, gene_id) %>%
  summarise(n = n(),
            range = ifelse(n() == 1, as.character(round(log2FoldChange, 2)),
                           paste(round(min(log2FoldChange), 2), round(max(log2FoldChange), 2), sep = '/')),
            ave = round(mean(log2FoldChange), 2),
            sd = ifelse(is.na(sd(log2FoldChange)), '', as.character(round(sd(log2FoldChange), 2)))) %>%
  unite(values, n, range, ave, sd, sep = '_') %>%
  spread(contrast, values) %>%
  ungroup() %>%
  inner_join(fac) %>%
  arrange(cat) %>%
  dplyr::select(cat, factor, gene = gene_id, early_vs_non, late_vs_non, late_vs_early) %>%
  separate(early_vs_non, sep = '_', into = c('n1', 'range1', 'ave1', 'sd1')) %>%
  separate(late_vs_non, sep = '_', into = c('n2', 'range2', 'ave2', 'sd2')) %>%
  separate(late_vs_early, sep = '_', into = c('n3', 'range3', 'ave3', 'sd3')) %>%
  mutate(factor = ifelse(duplicated(factor), '', factor)) %>%
  mutate(cat = ifelse(duplicated(cat), '', as.character(cat))) %>%
  xtable(caption = '\\textbf{Significant peaks of adipogenic factors on autophagy genes.}',
         label = 'tab:dep_fc_genes',
         align = 'clllcccccccccccc') %>%
  print(floating = TRUE,
        include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        sanitize.text.function = identity,
        comment = FALSE,
        include.colnames=FALSE,
        add.to.row = list(pos = list(0, 2, 5, 10, 6, 11),
                          command = c(header, rep('\\cmidrule{2-15} ', 3), rep('\\midrule ', 2))),
        file = 'manuscript/tables/dep_fc_genes.tex')
