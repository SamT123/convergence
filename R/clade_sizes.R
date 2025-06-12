
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

getAllCladeSizes = function(
    tree_and_sequences,
    tree_tibble_for_node_search = tree_and_sequences$tree_tibble,
    ps = F
){
  aa_clade_sizes = rep(
    list(NULL),
    length(unique(unlist(tree_tibble_for_node_search$aa_mutations_nonsyn)))
  ) %>%
    setNames(unique(unlist(tree_tibble_for_node_search$aa_mutations_nonsyn)))

  syn_nt_clade_sizes = rep(
    list(NULL),
    length(unique(unlist(tree_tibble_for_node_search$nt_mutations_syn)))
  ) %>%
    setNames(unique(unlist(tree_tibble_for_node_search$nt_mutations_syn)))

  tree_and_sequences$tree_tibble$tips_per_node = c(
    rep(1, ape::Ntip(tree_and_sequences$tree)),
    castor::count_tips_per_node(tree_and_sequences$tree)
  )

  tips_per_node = setNames(
    tree_and_sequences$tree_tibble$tips_per_node,
    as.character(tree_and_sequences$tree_tibble$node)
  )

  tips_per_node_tree_tibble_for_node_search = tips_per_node[
    tree_tibble_for_node_search$node
  ]

  for (i in seq_along(tree_tibble_for_node_search$aa_mutations_nonsyn)){

    mutations = unlist(tree_tibble_for_node_search$aa_mutations_nonsyn[[i]])

    aa_clade_sizes[mutations] = purrr::map(
      aa_clade_sizes[mutations],
      append,
      tips_per_node_tree_tibble_for_node_search[[i]]
    )
  }

  for (i in seq_along(tree_tibble_for_node_search$nt_mutations_syn)){

    mutations = unlist(tree_tibble_for_node_search$nt_mutations_syn[[i]])

    syn_nt_clade_sizes[mutations] = purrr::map(
      syn_nt_clade_sizes[mutations],
      append,
      tips_per_node_tree_tibble_for_node_search[[i]]
    )
  }

  all_syn_nt_clade_sizes = unlist(syn_nt_clade_sizes)
  mean_log2_syn_clade_size = mean(log2(all_syn_nt_clade_sizes))
  q0.8_log2_syn_clade_size = quantile(log2(all_syn_nt_clade_sizes), 0.8)

  makeTibble = function(clade_sizes, ps){
    tib = tibble(
      mutation = names(clade_sizes),
      clade_sizes = clade_sizes
    ) %>%
      mutate(
        n_clades = purrr::map_int(clade_sizes, length),

        mean_log2_clade_size = purrr::map_dbl(
          clade_sizes,
          ~mean(log2(.x))
        ),
        mean_log2_clade_size_diff = mean_log2_clade_size - mean_log2_syn_clade_size,

        q0.8_log2_clade_size = purrr::map_dbl(
          clade_sizes,
          ~quantile(log2(.x), 0.8)
        ),
        q0.8_log2_clade_size_diff = q0.8_log2_clade_size - q0.8_log2_syn_clade_size
      )


    if (ps){

      quantile_f = function(x) {quantile(x, 0.8)}
      mean_f = function(x) {mean(x)}
      n_resamples = 1e3

      bootPop = function(sample_size, null_pop, n_boots, summary_fun){
        sample(
          x = null_pop,
          size = sample_size*n_boots,
          replace = T
        ) %>%
          matrix(
            nrow = sample_size,
            ncol = n_boots
          ) %>%
          apply(2, summary_fun)
      }

      efficient_map = function(v, f, verbose = T, ...){
        v_unq = unique(v)
        if (verbose) message('length = ', length(v_unq))
        args = list(...)
        o = list()
        for (i in seq_along(v_unq)){
          o[[i]] = f(v_unq[[i]], ...)
        }
        o[match(v, v_unq)]
      }


      tib = tib %>%
        mutate(
          q0.8_log2_clade_size_diff_p = efficient_map(
            n_clades,
            bootPop,
            verbose = T,
            null_pop = log2(all_syn_nt_clade_sizes),
            n_boots = n_resamples,
            summary_fun = quantile_f
          ) %>%
            purrr::map2_dbl(
              q0.8_log2_clade_size,
              ~mean(.x>=.y)
            ),

          mean_log2_clade_size_diff_p = efficient_map(
            n_clades,
            bootPop,
            verbose = T,
            null_pop = log2(all_syn_nt_clade_sizes),
            n_boots = n_resamples,
            summary_fun = mean_f
          ) %>%
            purrr::map2_dbl(
              mean_log2_clade_size,
              ~mean(.x>=.y)
            )
        )
    }
    tib
  }


  list(
    aas = makeTibble(aa_clade_sizes, ps = ps),
    syn_nucs = makeTibble(syn_nt_clade_sizes, ps = ps)
  )
}


#' Distribution of clade sizes compared to synonymous expectation
#'
#' @param tree_and_sequences ...
#' @param tree_tibble_for_node_search Perhaps you want to restrict node search to a subset of the tree in tree_and_sequences. Because tree_and_sequences$tree is used to count descendant tips, the "intact" tree_and_sequences (i.e. without tree_and_sequences being subsetted, meaning `tree` corresponds exactly to `tree_tibble`) must also be passed
#' @param aa_substitutions ...
#' @param syn_nuc_mutations ...
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
    syn_nuc_mutations,
    n_bootstraps = 100,
    max_log2_clade_size = 10,
    log2_clade_size_resolution = 0.1
){

  all_clade_sizes = getAllCladeSizes(
    tree_and_sequences,
    tree_tibble_for_node_search
  )

  aa_substitution_clade_sizes = all_clade_sizes$aas %>%
    {setNames(.$clade_sizes, .$mutation)}



  syn_nuc_mutation_clade_sizes = all_clade_sizes$syn_nucs %>%
    {setNames(.$clade_sizes, .$mutation)}

  all_synonymous_clade_sizes = unlist(all_clade_sizes$syn_nucs$clade_sizes)


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

  for (m in syn_nuc_mutations){
    focal_cdfs = add_row(
      focal_cdfs,
      mutation = m,
      aa_or_nt = "nt",
      cdf = getCDF(
        log2(syn_nuc_mutation_clade_sizes[[m]]),
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

    for (s in aa_substitutions){
      focal_cdfs = add_row(
        focal_cdfs,
        mutation = s,
        aa_or_nt = "aa",
        cdf = getCDF(
          log2(
            sample(
              all_synonymous_clade_sizes,
              size = length(aa_substitution_clade_sizes[[s]]),
              replace = T
            )
          ),
          resolution = log2_clade_size_resolution,
          cdf_range = c(0, max_log2_clade_size)
        ),
        real_or_bootstrap = "bootstrap",
        bootstrap_replicate_id = i
      )
    }

    for (m in syn_nuc_mutations){
      focal_cdfs = add_row(
        focal_cdfs,
        mutation = m,
        aa_or_nt = "nt",
        cdf = getCDF(
          log2(
            sample(
              all_synonymous_clade_sizes,
              size = length(syn_nuc_mutation_clade_sizes[[m]]),
              replace = T
            )
          ),
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
