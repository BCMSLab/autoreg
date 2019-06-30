library(tidyverse)
library(reshape2)
library(tidygraph)
library(ggraph)

interactions <- read_rds('data/interactions.rds')
interactions <- interactions$tf_tf
factor_targets <- read_rds('data/factor_targets.rds')

go_annotation <- read_rds('data/go_annotation.rds')
tf_annotation <- read_rds('data/tf_annotation.rds')

nds <- list(Autophagy = intersect(go_annotation$SYMBOL, tf_annotation$SYMBOL),
            Adipogenic = c('Ctcf', 'Cebpb', 'Pparg', 'Rxrg', 'Ep300', 'Med1')) %>%
  melt() %>%
  setNames(c('name', 'category')) %>%
  mutate(name = toupper(name))

string <- interactions %>%
  dplyr::select(from = preferredName_A, to = preferredName_B) %>%
  unique()

chip <- factor_targets %>%
  filter(geneId %in% unique(unlist(string))) %>%
  mutate(annotation = str_split(annotation, ' ', simplify = TRUE)[, 1]) %>%
  filter(annotation %in% c('Promoter')) %>%
  dplyr::select(from = factor, to = geneId) %>%
  unique()

set.seed(12345)
(list(string = string, chip = chip) %>%
  bind_rows(.id = 'type') %>%
  dplyr::select(from, to, type) %>%
  mutate_all(toupper) %>%
  as_tbl_graph() %>%
  left_join(nds) %>%
  ggraph() +
  geom_edge_link(color = 'darkgray') +
  geom_node_point(size = 10, aes(color = category)) +
  geom_node_text(aes(label = name)) +
  facet_wrap(~type, scales = 'free') +
  theme_graph() +
  labs(color = '') +
  theme(legend.position = 'top')) %>%
  ggsave(plot = .,
         filename = 'manuscript/figures/string_chip_networks.png',
         width = 20, height = 20, units = 'cm')

  