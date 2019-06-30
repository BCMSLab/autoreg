ddcor <- read_rds('data/dgca.rds')

gene_lists <- ddcor %>%
  filter(stage %in% c('stage_1_cor', 'stage_3_cor')) %>%
  group_by(Gene2, stage, Classes) %>%
  summarise(genes = list(unique(Gene1))) %>%
  ungroup() %>%
  mutate(id = as.character(row_number()))

go_children_id <- unlist(as.list(GOBPCHILDREN['GO:0006914']), use.names = FALSE) %>%
  map(go2term) %>%
  melt() %>%
  dplyr::select(ID = go_id, Term)

go_children <- AnnotationDbi::select(
  org.Mm.eg.db,
  go_children_id$ID,
  'SYMBOL', 'GO') %>%
  select(term = GO, gene = SYMBOL)

go_ids <- tibble(ID = c('GO:0010507', 'GO:0010508'),
                 term = c('negative regulation of autophagy', 'positive regulation of autophagy'))

go_genes <- AnnotationDbi::select(
  org.Mm.eg.db,
  go_ids$ID,
  'SYMBOL', 'GO') %>%
  dplyr::select(term = GO, gene = SYMBOL) %>%
  unique() %>%
  filter(gene %in% go_annotation$SYMBOL)

names(gene_lists$Gene2) <- 1:length(gene_lists$genes)

comp <- map(gene_lists$genes, function(x) {
  ob <- enricher(x, TERM2GENE = (go_genes))
  if(!is_null(ob)) {
    ob@result
  }
}) %>%
  bind_rows(.id = 'id') %>%
  left_join(go_ids) %>%
  left_join(gene_lists) %>%
  filter(p.adjust < .2)

View(comp)
table(go_genes$gene %in% go_annotation$SYMBOL)

ddcor %>%
  left_join(go_genes, by = c('Gene1'='gene')) %>%
  na.omit() %>%
  filter(pValDiff_adj < .1) %>%
  group_by(Gene2, stage, Classes, term) %>%
  summarise(n()) %>%
  View
