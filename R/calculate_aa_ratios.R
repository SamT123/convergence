#' Calculate log convergence ratios
#'
#' @param t_and_s_filtered a list containing `tree` and `tree_tibble`
#' @param nuc_rates from `getTreeSizeAndNucRates`
#' @param tree_size_fn from `getTreeSizeAndNucRates`
#' @param sampling_function from `getTreeSizeAndNucRates`
#' @param positions aa positions for which to calculate log convergence
#' @param calculate_p_values controls p-value calculation. Can be TRUE (calculate for all positions),
#'   FALSE (skip p-values), or a list with optional elements: `aas` (positions for amino acid p-values),
#'   `nucs` (positions for nucleotide p-values), `min_LCR` (minimum log convergence ratio threshold),
#'   `max_LCR` (maximum log convergence ratio threshold). Default TRUE
#' @param p_alternative alternative hypothesis for p-value calculation. One of "greater" (default),
#'   "less", or "two.sided"
#' @param verbose report on progress? default TRUE
#'
#'@export
getSubstitutionRatios = function(
  t_and_s_filtered,
  nuc_rates,
  tree_size_fn,
  sampling_function,
  positions,
  calculate_p_values = TRUE,
  p_alternative = "greater",
  verbose = TRUE
) {
  # fmt: skip
  if (verbose) message("Counting codons...")
  codon_table = getCodonTable(
    t_and_s_filtered$tree_tibble,
    positions
  )

  codon_table_filtered = filter(codon_table, at %in% positions)

  # fmt: skip
  if (verbose) message("Calculating expected mutation counts...")
  codon_table_with_expected_n = addExpectedMutationsToCodonTable(
    t_and_s_filtered$tree_tibble,
    codon_table_filtered,
    nuc_rates,
    tree_size_fn
  )

  # fmt: skip
  if (verbose) message("Counting observed mutation counts...")
  mutation_table = addObservedMutationCounts(
    codon_table_with_expected_n,
    t_and_s_filtered$tree_tibble
  )

  # fmt: skip
  if (verbose) message("Summarising...")
  mutation_table = summariseMutationTableToAAsAndSyns(mutation_table)

  # fmt: skip
  if (isTRUE(calculate_p_values)) {
    p_value_positions = list(
      aas = unique(mutation_table[["aas"]][["at"]]),
      nucs = unique(mutation_table[["nucs"]][["at"]])
    )
    min_LCR = -Inf
    max_LCR = +Inf
  } else if (is.list(calculate_p_values)){

    if (!all(names(calculate_p_values) %in% c("aas", "nucs", "min_LCR", "max_LCR"))){
      stop("")
    }

    p_value_positions = list()

    if (!"aas" %in% names(calculate_p_values)){
      p_value_positions[["aas"]] = positions
    } else {
      p_value_positions[["aas"]] = calculate_p_values[["aas"]]
    }

    if (!"nucs" %in% names(calculate_p_values)){
      p_value_positions[["nucs"]] = positions
    } else {
      p_value_positions[["nucs"]] = calculate_p_values[["nucs"]]
    }

    if (!"min_LCR" %in% names(calculate_p_values)){
      min_LCR = -Inf
    } else {
      min_LCR = calculate_p_values[["min_LCR"]]
    }

    if (!"max_LCR" %in% names(calculate_p_values)){
      max_LCR = Inf
    } else {
      max_LCR = calculate_p_values[["max_LCR"]]
    }
    
  } else if (isFALSE(calculate_p_values)){
    return(mutation_table)
  } else{
    stop("calculate_p_values should be TRUE, FALSE, or list")
  }

  # fmt: skip
  if (verbose) message("Calculating p-values...")

  mutation_table = purrr::map2(
    mutation_table[c("aas", "nucs")],
    p_value_positions[c("aas", "nucs")],
    function(mutation_table, positions) {
      bind_rows(
        addPValuesToMutationTable(
          filter(
            mutation_table,
            at %in% positions | dplyr::between(log2(ratio), min_LCR, max_LCR)
          ),
          nuc_rates,
          sampling_function,
          alternative = p_alternative
        ),
        filter(mutation_table, !at %in% positions) %>% mutate(p = NA)
      )
    }
  )

  mutation_table
}


# gets all codons at each amino acid position
getCodonTable = function(tree_tibble, positions, min_prop = 0, min_count = 1) {
  codons = purrr::map(
    positions,
    function(aa_pos) {
      codons = substr(
        tree_tibble[["reconstructed_dna_sequence"]],
        3 * (aa_pos - 1) + 1,
        aa_pos * 3
      )
      counts = sort(table(codons))
      counts = counts[!stringr::str_detect(names(counts), "N")]
      props = counts / sum(counts)
      props = props[props >= min_prop & counts >= min_count]
      whiches = purrr::map(
        names(props),
        ~ tree_tibble$node[which(codons == .x)]
      )

      list(
        codons = names(props),
        props = as.numeric(unname(props)),
        nodes = whiches
      )
    }
  )

  codon_table = tibble::tibble(
    at = positions,
    codon = purrr::map(codons, "codons"),
    props = purrr::map(codons, "props"),
    nodes = purrr::map(codons, "nodes")
  ) %>%
    tidyr::unnest(c(codon, props, nodes))

  codon_table$tip_nodes = purrr::map(
    codon_table$nodes,
    ~ .x[!is.na(tree_tibble$dna_sequence[match(.x, tree_tibble$node)])]
  )

  codon_table
}


splitFourFoldSynPositionsFromSequenceTable = function(
  sequences_at_four_fold_syn_sites,
  tip_nodes
) {
  sequences_at_four_fold_syn_sites_idx = sequences_at_four_fold_syn_sites[
    sequences_at_four_fold_syn_sites$node %in% tip_nodes,
    -ncol(sequences_at_four_fold_syn_sites)
  ]

  sequences_at_four_fold_syn_sites_idx = setNames(
    apply(
      sequences_at_four_fold_syn_sites_idx,
      2,
      \(x) names(sort(table(x), decreasing = T))[[1]]
    ),
    colnames(sequences_at_four_fold_syn_sites_idx)
  )

  four_fold_syn_sites_idx = split(
    as.integer(stringr::str_sub(
      names(sequences_at_four_fold_syn_sites_idx),
      4
    )),
    sequences_at_four_fold_syn_sites_idx
  )

  four_fold_syn_sites_idx = four_fold_syn_sites_idx[c("A", "T", "C", "G")]

  four_fold_syn_sites_idx
}


addExpectedMutationsToCodonTable = function(
  tree_tibble,
  codon_table,
  nuc_rates,
  tree_size_fn
) {
  # make sequence table for four fold synonymous sites in passed tree_tibble
  four_fold_syn_sites = getFourFoldSynPositions(
    list(tree_tibble = filter(tree_tibble, !is.na(reconstructed_dna_sequence))),
    threshold = 0.95
  ) *
    3

  sequences_at_four_fold_syn_sites = purrr::map(
    four_fold_syn_sites,
    ~ substr(tree_tibble$reconstructed_dna_sequence, .x, .x)
  ) %>%
    setNames(paste0("pos", four_fold_syn_sites)) %>%
    bind_cols() %>%
    cbind(node = tree_tibble$node)

  # initialise columns
  codon_table$tree_size = NA

  codon_table$mutations_per_site_interp =
    rep(list(NULL), nrow(codon_table))

  for (idx in seq_len(nrow(codon_table))) {
    # get sites which are four fold syn at the nodes where the codon being considered is present
    four_fold_syn_sites_idx = splitFourFoldSynPositionsFromSequenceTable(
      sequences_at_four_fold_syn_sites,
      unique(
        c(
          codon_table[["nodes"]][[idx]],
          codon_table[["tip_nodes"]][[idx]]
        )
      )
    )

    # get tree size estimate + sample for the codon being considered, and interpolate nuc rates for the codon
    nuc_counts_idx = getNucCounts(
      list(
        tree_tibble = filter(
          tree_tibble,
          node %in% codon_table[["nodes"]][[idx]]
        )
      ),
      four_fold_syn_nuc_positions = four_fold_syn_sites_idx
    )

    tree_size_idx = tree_size_fn(nuc_counts = nuc_counts_idx)

    thiscodon_nuc_rates = nuc_rates

    thiscodon_nuc_rates$mutations_per_site =
      thiscodon_nuc_rates$mutations_per_site_rate * tree_size_idx

    codon_table$tree_size[[idx]] = tree_size_idx

    codon_table$mutations_per_site_interp[[idx]] = setNames(
      thiscodon_nuc_rates$mutations_per_site,
      thiscodon_nuc_rates$from_to
    )
  }

  codon_table
}

addObservedMutationCounts = function(codon_table, tree_tibble) {
  mutation_table = codon_table %>%
    ungroup() %>%
    mutate(at = purrr::map(at, ~ seq((.x - 1) * 3 + 1, .x * 3))) %>%
    tidyr::unnest(at) %>%
    mutate(
      at_m3 = 1 + ((at - 1) %% 3),
      from = substr(codon, at_m3, at_m3),
      to = list(c("A", "T", "C", "G")),
    ) %>%
    tidyr::unnest(to) %>%
    filter(from != to) %>%
    mutate(
      from_to = paste0(from, to),
      nt_mutation = paste0(from, at, to),

      mutations_per_site_interp = purrr::map2_dbl(
        mutations_per_site_interp,
        from_to,
        ~ unlist(.x[[.y]])
      ),
    ) %>%
    select(nt_mutation, from_to, from, at, to, everything(), -at_m3) %>%
    ungroup()

  observed_mutations = tree_tibble %>%
    select(
      nt_mutations,
      nt_mutations_is_single,
      codon_changes
    ) %>%
    tidyr::unnest(c(
      nt_mutations,
      nt_mutations_is_single,
      codon_changes
    )) %>%
    filter(nt_mutations_is_single) %>%
    group_by(nt_mutations, codon_changes) %>%
    reframe(
      n = n(),
      codon_from = stringr::str_sub(unique(codon_changes), 1, 3),
      .groups = "drop"
    ) %>%
    ungroup() %>%
    mutate(
      from = stringr::str_sub(nt_mutations, 1, 1),
      at = stringr::str_sub(nt_mutations, 2, -2),
      to = stringr::str_sub(nt_mutations, -1, -1),
      nt_mutation = nt_mutations
    )

  observed_mutations$nt_mutation = as.character(observed_mutations$nt_mutation)

  mutation_table = left_join(
    mutation_table,
    select(observed_mutations, nt_mutation, codon_from, n),
    by = c("nt_mutation" = "nt_mutation", "codon" = "codon_from")
  )

  mutation_table$n[is.na(mutation_table$n)] = 0

  mutation_table = mutation_table %>%
    mutate(
      expected_n = mutations_per_site_interp,
      aa_mutation = paste0(
        Biostrings::GENETIC_CODE[codon],
        ceiling(at / 3),
        Biostrings::GENETIC_CODE[
          purrr::pmap_chr(list(codon, at, to), function(codon, at, to) {
            applySubstitutions(
              codon,
              paste0("X", 1 + (unique(at) - 1) %% 3, unique(to))
            )
          })
        ]
      ),
      ratio = n / mutations_per_site_interp
    )

  mutation_table
}

summariseMutationTableToAAsAndSyns = function(mutation_table) {
  # fmt: skip
  if (any(duplicated(paste(mutation_table$nt_mutation, mutation_table$codon)))) {
    stop("duplicated codon - mutation pairs in mutation table!")
  }

  mutation_table = mutation_table %>%
    group_by(nt_mutation, codon) %>%
    summarise(
      expected_n = sum(mutations_per_site_interp),
      n = sum(n),
      aa_mutation = unique(aa_mutation),
      is_syn = Biostrings::GENETIC_CODE[codon] ==
        Biostrings::GENETIC_CODE[
          applySubstitutions(
            codon,
            paste0("X", 1 + (unique(at) - 1) %% 3, unique(to))
          )
        ],
      tree_size = unique(tree_size),
      .groups = "drop"
    ) %>%
    group_by(aa_mutation, nt_mutation, codon) %>%
    mutate(
      components = list(list(
        nt_mutation = (nt_mutation),
        nt_mutation_class = paste0(
          stringr::str_sub(nt_mutation, 1, 1),
          stringr::str_sub(nt_mutation, -1, -1)
        ),
        tree_size = (tree_size),
        expected_n = (expected_n),
        n = (n)
      ))
    ) %>%
    ungroup()

  aa_mutation_table = mutation_table %>%
    filter(!is_syn) %>%
    group_by(aa_mutation) %>%
    summarise(
      components = list(components),
      expected_n = sum(expected_n),
      n = sum(n),
      tree_size = sum(tree_size),
      ratio = n / expected_n,
      .groups = "drop"
    ) %>%
    arrange(-ratio) %>%
    mutate(
      aa_mutation = factor(aa_mutation, levels = unique(aa_mutation)),
      from = stringr::str_sub(aa_mutation, 1, 1),
      at = stringr::str_sub(aa_mutation, 2, -2),
      to = stringr::str_sub(aa_mutation, -1, -1)
    )

  nt_mutation_table = mutation_table %>%
    mutate(
      components = list(components), # to match aa_mutation_table
      ratio = n / expected_n,
      from = stringr::str_sub(nt_mutation, 1, 1),
      at = stringr::str_sub(nt_mutation, 2, -2),
      to = stringr::str_sub(nt_mutation, -1, -1)
    )

  list(
    aas = aa_mutation_table,
    nucs = nt_mutation_table
  )
}
