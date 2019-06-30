factor_targets <- read_rds('data/factor_targets.rds')
go_annotation <- read_rds('data/go_annotation.rds')

map(list(all = factor_targets,
         autophagy = factor_targets %>% filter(geneId %in% unique(go_annotation$SYMBOL))),
    function(df) {
      df %>%
      with(split(., stage)) %>%
        map(function(x) {
          targets <- map(c('PPARG', 'CEBPB'), function(y) {
            filter(x, factor == y) %>% pull(geneId) %>% unique()
          })
          tot <- length(unlist(targets))
          list(PPARG = length(setdiff(targets[[1]], targets[[2]]))/tot,
               CEBPB = length(setdiff(targets[[2]], targets[[1]]))/tot,
               COMMON = length(intersect(targets[[1]], targets[[2]]))/tot)
        }) %>%
        bind_rows(.id = 'stage')
    }) %>%
  bind_rows(.id = 'type')
