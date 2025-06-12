#' Get information about a tree required for calculating convergence scores
#'
#' @param tree_and_sequences a list containing `tree` and `tree_tibble`
#' @param noise_model one of the 'noise_models'
#' @param n_iter iterations for rstan
#' @param n_chains chains for rstan
#' @param excluded_positions positions to exclude for synonymous rate variation calculation
#' @param use_fitted_rates should the fitted mutation rates be used, rather than simply taking the mean of the observed values?
#' @param seed seed for rstan
#'
#'@export
getTreeSizeAndNucRates = function(
    tree_and_sequences,
    noise_model,
    n_iter = 2000,
    n_chains = 3,
    excluded_positions = c(),
    use_fitted_rates = T,
    seed = NULL
) {

  nuc_counts = getNucCounts(
    tree_and_sequences,
    excluded_positions=excluded_positions
  )

  nuc_counts_with_model = addMutationRateModel(
    nuc_counts,
    model = noise_model,
    n_iter = n_iter,
    n_chains = n_chains,
    seed = seed
  )
  nuc_counts = nuc_counts_with_model[["nuc_counts"]]
  model_fit = nuc_counts_with_model[["model_fit"]]

  if (use_fitted_rates){

    nuc_counts$mutations_per_site = purrr::pmap_dbl(
        list(
          .x = nuc_counts$parameters,
          .y = nuc_counts$mutations_per_site_measured
        ),
        \(.x,.y) {
          2**(log2(.y)+.x$centrality)
        }
      )

    nuc_counts$parameters = purrr::map(
      nuc_counts$parameters,
      \(p) {
        p$centrality = 0
        p
      }
    )

  } else {
    nuc_counts$mutations_per_site = nuc_counts$mutations_per_site_measured
  }

  tree_size_ratios = getTreeSizeNucWeightings(nuc_counts)

  tree_size_fn = function(
    nuc_counts
    # tree_size_ratios
  ) {
    tree_size_components = nuc_counts$mutations_per_site_measured *
      tree_size_ratios$ratio[match(nuc_counts$from_to, tree_size_ratios$from_to)]

    sum(tree_size_components)
  }

  tree_size = tree_size_fn(nuc_counts)

  nuc_rates = getNucRates(nuc_counts, tree_size)

  rm(list = c("tree_and_sequences"))
  list(
    tree_size = tree_size,
    nuc_rates = nuc_rates,
    tree_size_fn = tree_size_fn,
    noise_model_fit = model_fit,
    sampling_function = noise_model$sampling_function
  )
}

getNucCounts = function(
    t_and_s,
    four_fold_syn_nuc_positions = NULL,
    excluded_positions = c()
) {

  rates = tibble(
    from = character(),
    to = character(),
    from_to = character(),
    n_mutations_individual_positions = list(),
    n_mutations = integer(),
    mutations_per_site_measured = numeric(),
    n_sites = integer(),
    n_seqs = integer(),
    n_nodes = integer()
  )

  fixations = t_and_s$tree_tibble %>%
    select(nt_mutations_syn) %>%
    tidyr::unnest(nt_mutations_syn) %>%
    mutate(
      from = stringr::str_sub(nt_mutations_syn, 1, 1),
      at = stringr::str_sub(nt_mutations_syn, 2, -2) %>% as.integer(),
      to = factor(
        stringr::str_sub(nt_mutations_syn, -1, -1),
        levels = c("A", "T", "C", "G")
      )
    )

  if (missing(four_fold_syn_nuc_positions)){
    four_fold_syn_nuc_positions = getSplitFourFoldSynPositionsFromTreeAndSequences(
      t_and_s,
      proportion_syn_threshold = 0.95,
      identity_threshold = 0.9
    )
  }

  n_s = sum(!is.na(t_and_s$tree_tibble$dna_sequence))
  n_n = nrow(t_and_s$tree_tibble)

  for (fr in names(four_fold_syn_nuc_positions)){

    rates_nt = fixations %>%
      filter(
        at %in% four_fold_syn_nuc_positions[[fr]],
        !at %in% excluded_positions
      ) %>%
      mutate(at = factor(at, levels = four_fold_syn_nuc_positions[[fr]])) %>%
      group_by(to, .drop = F) %>%
      filter(as.character(to) != fr) %>%
      summarise(
        n_mutations_individual_positions = list(table(at)),
        n_mutations = n()
      ) %>%
      complete(to, fill = list(n_mutations = 0)) %>%
      filter(as.character(to) != fr) %>%
      mutate(
        n_sites = length(four_fold_syn_nuc_positions[[fr]]),
        mutations_per_site_measured = n_mutations / n_sites
      )

    rates = add_row(
      rates,
      from = rep(fr, 3),
      to = as.character(rates_nt$to),
      from_to = paste0(rep(fr, 3), as.character(rates_nt$to)),
      n_mutations_individual_positions =
        rates_nt$n_mutations_individual_positions,
      n_mutations = rates_nt$n_mutations,
      mutations_per_site_measured = rates_nt$mutations_per_site_measured,
      n_sites = rates_nt$n_sites,
      n_seqs = rep(n_s, 3),
      n_nodes = rep(n_n, 3)
    )
  }

  rates = rates %>%
    filter(from != to)

  rates
}


addMutationRateModel = function(
    nuc_counts,
    model,
    n_chains,
    n_iter,
    seed = NULL
) {
  nuc_counts_long = nuc_counts %>%
    mutate(
      n = n_mutations_individual_positions,
      mutation = purrr::pmap(
        list(from, to, n),
        \(f,t,x) paste0(f, names(x), t)
      ),
      mutation_class = paste0(from, to),
      expected_n = mutations_per_site_measured
    ) %>%
    tidyr::unnest(c(mutation, n)) %>%
    select(
      from,
      to,
      from_to,
      mutation_class,
      mutation,
      expected_n,
      n
    )

  model_fit = fitNoiseModel(
    nuc_counts_long,
    model$specification,
    n_chains = n_chains,
    n_iter = n_iter,
    seed = seed
  )

  nuc_counts$parameters =
    model$extract_parameters_function(
      model_fit
      )[nuc_counts$from_to]

  list(
    nuc_counts = nuc_counts,
    model_fit = model_fit
  )
}


getFourFoldSynPositions = function(t_and_s, threshold){

  four_fold_syn_codons = names(Biostrings::GENETIC_CODE)[purrr::map_lgl(
    names(Biostrings::GENETIC_CODE),
    function(cod){
      length(unique(Biostrings::GENETIC_CODE[
        paste0(stringr::str_sub(cod, 1, 2), c("A", "T", "C", "G"))])) == 1
    }
  )] %>%
    unique()

  t_and_s$tree_tibble$codons = purrr::map(
    t_and_s$tree_tibble$reconstructed_dna_sequence,
    ~substring(.x, seq(1, nchar(.x), 3), seq(3, nchar(.x), 3))
  )

  t_and_s$tree_tibble$four_fold_syn_codons = purrr::map(
    t_and_s$tree_tibble$codons,
    ~which(.x %in% four_fold_syn_codons)
  )

  four_fold_syn_frequencies = table(unlist(
    t_and_s$tree_tibble$four_fold_syn_codons
  )) / nrow(t_and_s$tree_tibble)

  majority_four_fold_syn = names(four_fold_syn_frequencies)[
    four_fold_syn_frequencies > threshold
  ]

  majority_four_fold_syn = as.integer(majority_four_fold_syn)
}


getSplitFourFoldSynPositionsFromTreeAndSequences = function(
    t_and_s,
    proportion_syn_threshold,
    identity_threshold
){
  four_fold_syn_nuc_positions = getFourFoldSynPositions(
    t_and_s,
    proportion_syn_threshold
  ) * 3

  four_fold_syn_nuc_position_aas = purrr::map(
    four_fold_syn_nuc_positions,
    function(pos){
      nts = substr(
        t_and_s$tree_tibble$dna_sequence,
        pos,
        pos
      )
      nts = nts[!is.na(nts)]
      nts = sort(table(nts), decreasing = T) /
        sum(nts %in% c("A", "T", "C", "G"))
      list(nt = names(nts)[[1]],
           prop = nts[[1]])
    }
  )

  four_fold_syn_nuc_positions = split(
    four_fold_syn_nuc_positions[
      purrr::map_dbl(four_fold_syn_nuc_position_aas, "prop") > identity_threshold],
    purrr::map_chr(four_fold_syn_nuc_position_aas, "nt")[
      purrr::map_dbl(four_fold_syn_nuc_position_aas, "prop") > identity_threshold]
  )

  four_fold_syn_nuc_positions
}


getTreeSizeNucWeightings = function(nuc_counts){
  nuc_counts %>%
    transmute(
      from = from,
      to = to,
      from_to = from_to,
      n_sites = n_sites,
      ratio = n_sites
    )
}

# getTreeSize = function(
#   nuc_counts,
#   tree_size_ratios
# ) {
#   tree_size_components = nuc_counts$mutations_per_site_measured *
#     tree_size_ratios$ratio[match(nuc_counts$from_to, tree_size_ratios$from_to)]
#
#   sum(tree_size_components)
# }


getNucRates = function(nuc_counts, tree_size){
  nuc_counts$mutations_per_site_rate = nuc_counts$mutations_per_site / tree_size
  nuc_counts
}

orderNucRateColumns = function(df) {
  select(
    df,
    cluster,
    from,
    to,
    from_to,
    n_mutations_individual_positions,
    n_mutations,
    n_sites,
    mutations_per_site,
    mutations_per_site_rate,
    mutations_per_site_interp,
    n_seqs,
    n_nodes
  ) %>%
    mutate(
      from_to = factor(from_to, levels = unlist(nt_substitutions)),
      from = factor(from, levels = c("A", "T", "C", "G")),
      to = factor(to, levels = c("A", "T", "C", "G"))
    )
}
