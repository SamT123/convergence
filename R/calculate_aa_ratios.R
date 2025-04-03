#'@export
getSubstitutionRatios = function(
    t_and_s_filtered,
    nuc_rates,
    tree_size_ratios,
    tree_size_fn,
    sampling_function,
    positions
){

  codon_table = getCodonTable(
    t_and_s_filtered$tree_tibble,
    positions
  )

  codon_table_filtered = filter(codon_table, at %in% positions)

  codon_table_with_expected_n = addExpectedMutationsToCodonTable(
    t_and_s_filtered$tree_tibble,
    codon_table_filtered,
    nuc_rates,
    tree_size_ratios,
    tree_size_fn
  )

  mutation_table = addObservedMutationCounts(
    codon_table_with_expected_n,
    t_and_s_filtered$tree_tibble
  )

  mutation_table = summariseMutationTableToAAsAndSyns(mutation_table)

  mutation_table = map(
    mutation_table,
    ~addPValuesToMutationTable(.x, nuc_rates, sampling_function)

  )

  mutation_table
}


# gets all codons with frequency >1% at each amino acid position
getCodonTable = function(tree_tibble, positions, min_prop = 0, min_count = 1){
  codons = map(
    positions,
    function(aa_pos){
      codons = substr(
        tree_tibble[["reconstructed_dna_sequence"]],
        3*(aa_pos-1)+1,
        aa_pos*3
      )
      counts = sort(table(codons))
      counts = counts[!str_detect(names(counts), "N")]
      props = counts / sum(counts)
      props = props[props >= min_prop & counts >= min_count] # remove this?
      whiches = map(names(props), ~tree_tibble$node[which(codons == .x)])

      list(
        codons = names(props),
        props = as.numeric(unname(props)),
        nodes = whiches
      )
    }
  )

  codon_table = tibble(
    at = positions,
    codon = map(codons, "codons"),
    props = map(codons, "props"),
    nodes = map(codons, "nodes")
  ) %>%
    unnest(c(codon, props, nodes))

  codon_table$tip_nodes = map(
    codon_table$nodes,
    ~.x[!is.na(tree_tibble$dna_sequence[match(.x, tree_tibble$node)])]
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
    apply(sequences_at_four_fold_syn_sites_idx, 2, \(x){names(sort(table(x)))[[1]]}),
    colnames(sequences_at_four_fold_syn_sites_idx)
  )

  four_fold_syn_sites_idx = split(
    as.integer(str_sub(names(sequences_at_four_fold_syn_sites_idx), 4)),
    sequences_at_four_fold_syn_sites_idx
  )

  four_fold_syn_sites_idx = four_fold_syn_sites_idx[c("A", "T", "C", "G")]

  four_fold_syn_sites_idx
}


addExpectedMutationsToCodonTable = function(
    tree_tibble,
    codon_table,
    nuc_rates,
    tree_size_ratios,
    tree_size_fn
){

  # make sequence table for four fold synonymous sites in passed tree_tibble
  four_fold_syn_sites = getFourFoldSynPositions(
    list(tree_tibble = filter(tree_tibble, !is.na(reconstructed_dna_sequence))),
    threshold = 0.95
  ) * 3

  sequences_at_four_fold_syn_sites = map(
    four_fold_syn_sites,
    ~substr(tree_tibble$reconstructed_dna_sequence, .x, .x)
  ) %>%
    setNames(paste0("pos", four_fold_syn_sites)) %>%
    bind_cols() %>%
    cbind(node = tree_tibble$node)

  # initialise columns
  codon_table$tree_size = NA

  codon_table$mutations_per_site_interp =
    rep(list(NULL), nrow(codon_table))

  for (idx in seq_len(nrow(codon_table))){
    # get sites which are four fold syn at the tips with the codon being considered
    four_fold_syn_sites_idx = splitFourFoldSynPositionsFromSequenceTable(
      sequences_at_four_fold_syn_sites,
      unique(
        c(
          codon_table[["nodes"]][[idx]],
          codon_table[["tip_nodes"]][[idx]]
        )
      )
    )

    # get tree size estimate + sample for the codon being considered, and interpolate mutation nuc rates for the codon
    nuc_counts_idx = getNucCounts(
      list(
        tree_tibble = filter(
          tree_tibble,
          node %in% codon_table[["nodes"]][[idx]]
        )
      ),
      four_fold_syn_nuc_positions = four_fold_syn_sites_idx
    )

    tree_size_idx = tree_size_fn(
      nuc_counts = nuc_counts_idx,
      tree_size_ratios = tree_size_ratios
    )

    thiscodon_nuc_rates = interpolateMutationRates(
      nuc_counts_idx,
      tree_size_idx,
      nuc_rates
    )

    codon_table$tree_size[[idx]] = tree_size_idx

    codon_table$mutations_per_site_interp[[idx]] = setNames(
      thiscodon_nuc_rates$mutations_per_site_interp,
      thiscodon_nuc_rates$from_to)

  }

  codon_table
}

addObservedMutationCounts = function(codon_table, tree_tibble){

  mutation_table = codon_table %>%
    ungroup() %>%
    mutate(at = map(at, ~seq((.x-1)*3+1, .x*3))) %>%
    unnest(at) %>%
    mutate(
      at_m3 = 1+((at-1)%%3),
      from = substr(codon, at_m3, at_m3),
      to = list(c("A", "T", "C", "G")),
    ) %>%
    unnest(to) %>%
    filter(from != to) %>%
    mutate(
      from_to = paste0(from, to),
      nt_mutation = paste0(from, at, to),

      mutations_per_site_interp = map2_dbl(
        mutations_per_site_interp,
        from_to,
        ~unlist(.x[[.y]])
      ),
    ) %>%
    select(nt_mutation, from_to, from, at, to, everything(), -at_m3) %>%
    ungroup()


  observed_mutations = tree_tibble %>%
    select(nt_mutations, nt_mutations_is_single, codon_changes, aa_mutations) %>%
    unnest(c(nt_mutations, nt_mutations_is_single, codon_changes, aa_mutations)) %>%
    filter(nt_mutations_is_single) %>%
    group_by(nt_mutations, codon_changes) %>%
    reframe(
      n = n(),
      aa_mutations = list(unique(aa_mutations)),
      codon_from = str_sub(unique(codon_changes), 1, 3),
      .groups = "drop"
    ) %>%
    ungroup() %>%
    mutate(
      from = str_sub(nt_mutations, 1, 1),
      at = str_sub(nt_mutations, 2, -2),
      to = str_sub(nt_mutations, -1, -1),
      aa_mutation = unlist(aa_mutations),
      nt_mutation = nt_mutations
    )

  if (!"aa_mutation" %in% colnames(observed_mutations)){
    observed_mutations$aa_mutation = character()
  }

  observed_mutations$nt_mutation = as.character(observed_mutations$nt_mutation)

  mutation_table = left_join(
    mutation_table,
    select(observed_mutations, nt_mutation, codon_from, n, aa_mutation),
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
          pmap_chr(list(codon, at, to), function(codon,at,to){applySubstitutions(
            codon,
            paste0("X", 1+(unique(at)-1) %% 3, unique(to))
          )})]
      ),
      ratio = n / mutations_per_site_interp
    )

  mutation_table
}

summariseMutationTableToAAsAndSyns  = function(mutation_table){

  syn_nt_mutations = mutation_table %>%
    mutate(
      is_syn = (
        Biostrings::GENETIC_CODE[codon] ==
          Biostrings::GENETIC_CODE[
            pmap_chr(
              list(codon, at, to),
              function(codon,at,to){convergence:::applySubstitutions(
                codon,
                paste0("X", 1+(unique(at)-1) %% 3, unique(to))
              )}
            )]
      )) %>%
    filter(is_syn)

  if (nrow(syn_nt_mutations) > 0){

    syn_nt_mutations = syn_nt_mutations%>%
      group_by(nt_mutation, codon) %>%
      summarise(
        expected_n = sum(mutations_per_site_interp),
        n = sum(n),
        tree_size = unique(tree_size),
        components = list(list(list(
          nt_mutation = (nt_mutation),
          nt_mutation_class = paste0(
            str_sub(nt_mutation, 1, 1),
            str_sub(nt_mutation, -1, -1)
          ),
          tree_size = (tree_size),
          expected_n = (expected_n)
        ))),
        .groups = "drop"
      ) %>%
      ungroup() %>%
      mutate(
        ratio = n / expected_n,
      ) %>%
      arrange(-ratio) %>%
      mutate(
        from = str_sub(nt_mutation, 1, 1),
        at = str_sub(nt_mutation, 2, -2),
        to = str_sub(nt_mutation, -1, -1)
      )
  }

  aa_mutation_table = mutation_table %>%
    distinct() %>%
    group_by(nt_mutation, codon) %>%
    summarise(
      expected_n = sum(mutations_per_site_interp),
      n = sum(n),
      aa_mutation = unique(aa_mutation),
      is_syn = Biostrings::GENETIC_CODE[codon] == Biostrings::GENETIC_CODE[
        convergence:::applySubstitutions(
          codon,
          paste0("X", 1+(unique(at)-1) %% 3, unique(to))
        )],
      tree_size = unique(tree_size),
      .groups = "drop"
    ) %>%
    filter(!is_syn) %>%
    distinct() %>%
    group_by(aa_mutation, nt_mutation, codon) %>%
    mutate(
      components = list(list(
        nt_mutation = (nt_mutation),
        nt_mutation_class = paste0(
          str_sub(nt_mutation, 1, 1),
          str_sub(nt_mutation, -1, -1)
        ),
        tree_size = (tree_size),
        expected_n = (expected_n)
      ))
    ) %>%
    group_by(aa_mutation) %>%
    summarise(
      components = list(components),
      expected_n = sum(expected_n),
      n = sum(n),
      tree_size = sum(tree_size),
      .groups = "drop"
    ) %>%
    mutate(
      ratio = n / expected_n
    ) %>%
    arrange(-ratio) %>%
    mutate(
      aa_mutation = factor(aa_mutation, levels = unique(aa_mutation)),
      from = str_sub(aa_mutation, 1, 1),
      at = str_sub(aa_mutation, 2, -2),
      to = str_sub(aa_mutation, -1, -1)
    )

  list(aas = aa_mutation_table, syn_nucs = syn_nt_mutations)
}

interpolateMutationRates = function(
    nuc_rates_by_cluster,
    cluster_tree_size,
    nuc_rates_overall
) {

  nuc_rates_by_cluster$mutations_per_site_interp = rep(
    NA,
    nrow(nuc_rates_by_cluster)
  )

  for (from in c("A", "T", "C", "G")){
    for (to in c("A", "T", "C", "G")){
      if (from == to) next
      row = which(
        nuc_rates_by_cluster$from == from & nuc_rates_by_cluster$to == to
      )

      nuc_rates_by_cluster$mutations_per_site_interp[[row]] =
        cluster_tree_size *
        nuc_rates_overall$mutations_per_site_rate[[
          which(nuc_rates_overall$from == from &
                  nuc_rates_overall$to == to)]]

    }
  }


  nuc_rates_by_cluster
}
