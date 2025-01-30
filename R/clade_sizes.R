
getCDF = function(
    values,
    cdf_range = range(values),
    resolution = diff(cdf_range)/1000
) {

  sample_points = seq(
    cdf_range[[1]],
    resolution*ceiling(cdf_range[[2]]/resolution),
    resolution
  )
  cdf = c()
  for (sp in sample_points){
    cdf[as.character(sp)] = mean(values<sp)
  }
  cdf
}


#' Distribution of clade sizes compared to synonymous expectation
#'
#' @param tree_and_sequences ...
#' @param tree_tibble_for_node_search Perhaps you want to restrict node search to a subset of the tree in tree_and_sequences. Because tree_and_sequences$tree is used to count descendant tips, the "intact" tree_and_sequences (i.e. without tree_and_sequences being subsetted, meaning `tree` corresponds exactly to `tree_tibble`) must also be passed
#' @param aa_substitutions ...
#' @param nt_mutations ...
#' @param n_bootstraps ...
#' @param max_log2_clade_size ...
#' @param log2_clade_size_resolution ...
#'
#' @return ...
#' @export
getCladeSizeCDFs = function(
    tree_and_sequences,
    tree_tibble_for_node_search = tree_and_sequences$tree_tibble,
    aa_substitutions,
    nt_mutations,
    n_bootstraps = 100,
    max_log2_clade_size = 10,
    log2_clade_size_resolution = 0.1
){


  tree_and_sequences$tree_tibble$tips_per_node = c(
    rep(1, Ntip(tree_and_sequences$tree)),
    castor::count_tips_per_node(tree_and_sequences$tree)
  )

  tips_per_node = setNames(
    tree_and_sequences$tree_tibble$tips_per_node,
    as.character(tree_and_sequences$tree_tibble$node)
  )

  all_synonymous_nodes = tree_tibble_for_node_search %>%
    unnest(nt_mutations_syn) %>%
    pluck("node")

  all_synonymous_clade_sizes = tips_per_node[as.character(all_synonymous_nodes)]

  aa_substitution_clade_sizes = list()

  for (s in substitutions){
    aa_substitution_nodes = tree_tibble_for_node_search %>%
      filter(map_lgl(aa_mutations, ~s %in% .x)) %>%
      pluck("node")

    aa_substitution_clade_sizes[[s]] = tips_per_node[as.character(aa_substitution_nodes)]
  }

  nt_mutation_clade_sizes = list()

  for (m in nt_mutations){
    nt_mutation_nodes = tree_tibble_for_node_search %>%
      filter(map_lgl(nt_mutations, ~m %in% .x)) %>%
      pluck("node")

    nt_mutation_clade_sizes[[m]] = tips_per_node[as.character(nt_mutation_nodes)]
  }

  focal_cdfs = tibble(
    mutation = character(),
    aa_or_nt = character(),
    cdf = numeric(),
    real_or_bootstrap = character(),
    bootstrap_replicate_id = integer()
  )

  synonymous_nt_cdfs = tibble(
    mutation = character(),
    aa_or_nt = character(),
    cdf = numeric(),
    real_or_bootstrap = character(),
    bootstrap_replicate_id = integer()
  )

  synonymous_nt_cdfs = add_row(
    synonymous_nt_cdfs,
    mutation = "Any synonymous mutation",
    aa_or_nt = "nt",
    cdf = getCDF(
      log2(all_synonymous_clade_sizes),
      resolution = log2_clade_size_resolution,
      cdf_range = c(0, max_log2_clade_size)
    ),
    real_or_bootstrap = "real",
    bootstrap_replicate_id = NA
  )

  for (s in aa_substitutions){
    focal_cdfs = add_row(
      focal_cdfs,
      mutation = s,
      aa_or_nt = "aa",
      cdf = getCDF(
        log2(aa_substitution_clade_sizes[[s]]),
        resolution = log2_clade_size_resolution,
        cdf_range = c(0, max_log2_clade_size)
      ),
      real_or_bootstrap = "real",
      bootstrap_replicate_id = NA
    )
  }

  for (m in nt_mutations){
    focal_cdfs = add_row(
      focal_cdfs,
      mutation = m,
      aa_or_nt = "nt",
      cdf = getCDF(
        log2(nt_mutation_clade_sizes[[m]]),
        resolution = log2_clade_size_resolution,
        cdf_range = c(0, max_log2_clade_size)
      ),
      real_or_bootstrap = "real",
      bootstrap_replicate_id = NA
    )
  }

  for (i in 1:n_bootstraps){

    synonymous_nt_cdfs = add_row(
      synonymous_nt_cdfs,
      mutation = "Any synonymous mutation",
      aa_or_nt = "nt",
      cdf = getCDF(
        log2(sample(all_synonymous_clade_sizes, replace = T)),
        resolution = log2_clade_size_resolution,
        cdf_range = c(0, max_log2_clade_size)
      ),
      real_or_bootstrap = "bootstrap",
      bootstrap_replicate_id = i
    )

    for (s in substitutions){
      focal_cdfs = add_row(
        focal_cdfs,
        mutation = s,
        aa_or_nt = "aa",
        cdf = getCDF(
          log2(sample(aa_substitution_clade_sizes[[s]], replace = T)),
          resolution = log2_clade_size_resolution,
          cdf_range = c(0, max_log2_clade_size)
        ),
        real_or_bootstrap = "bootstrap",
        bootstrap_replicate_id = i
      )
    }

    for (m in nt_mutations){
      focal_cdfs = add_row(
        focal_cdfs,
        mutation = m,
        aa_or_nt = "nt",
        cdf = getCDF(
          log2(sample(nt_mutation_clade_sizes[[m]], replace = T)),
          resolution = log2_clade_size_resolution
        ),
        real_or_bootstrap = "bootstrap",
        bootstrap_replicate_id = i
      )
    }
  }

  focal_cdfs$mutation = as_factor(focal_cdfs$mutation)

  focal_cdfs = focal_cdfs %>%
    mutate(
      clade_size = as.numeric(names(cdf)),
      p = unname(cdf)
    ) %>%
    select(-cdf)

  synonymous_nt_cdfs = synonymous_nt_cdfs %>%
    mutate(
      clade_size = as.numeric(names(cdf)),
      p = unname(cdf)
    ) %>%
    select(-cdf)

  to_intervals = function(cdfs){
    intervals =  cdfs %>%
      filter(real_or_bootstrap == "bootstrap") %>%
      group_by(
        mutation,
        clade_size
      ) %>%
      summarise(
        p_lower = quantile(p, c(0.05, 0.95))[1],
        p_upper = quantile(p, c(0.05, 0.95))[2],
        .groups = "drop"
      )

    left_join(
      cdfs %>%
        filter(real_or_bootstrap == "real") %>%
        select(-real_or_bootstrap, - bootstrap_replicate_id),
      intervals,
      by = c("mutation", "clade_size")
    )
  }

  raw = list(
    all_synonymous_mutations = synonymous_nt_cdfs,
    focal_substitutions_and_mutations = focal_cdfs
  )

  intervals = list(
    all_synonymous_mutations = to_intervals(synonymous_nt_cdfs),
    focal_substitutions_and_mutations = to_intervals(focal_cdfs)
  )

  list(
    raw = raw,
    intervals = intervals
  )
}
