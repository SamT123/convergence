# USHER ------------------------------------------------------------------------




# ------------------------------------------------------------------------------
#'@export
makeOverallTreeInfo = function(t_and_s){

  overall_nuc_counts = getNucCounts(t_and_s)
  tree_size_ratios = getTreeSizeRatios(overall_nuc_counts)
  overall_tree_size = getTreeSize(
    overall_nuc_counts,
    tree_size_ratios
  )
  overall_nuc_rates = getNucRates(overall_nuc_counts, overall_tree_size)

  tree_size_fn = function(
    nuc_counts,
    tree_size_ratios
  ) {

    size = 0

    for (nt in unique(nuc_counts$from_to)){
      size_increment = nuc_counts$mutations_per_site[nuc_counts$from_to == nt] *
        tree_size_ratios$ratio[tree_size_ratios$from_to == nt]
      size = size + size_increment
    }

    list(point = size)
  }

  rm(list = c("t_and_s"))
  list(
    tree_size_ratios = tree_size_ratios,
    overall_tree_size = overall_tree_size,
    overall_nuc_rates = overall_nuc_rates,
    fast_tree_size_fn = tree_size_fn
  )
}

getFourFoldSynPositions = function(t_and_s, threshold){

  four_fold_syn_codons = names(Biostrings::GENETIC_CODE)[map_lgl(
    names(Biostrings::GENETIC_CODE),
    function(cod){
      length(unique(Biostrings::GENETIC_CODE[
        paste0(str_sub(cod, 1, 2), c("A", "T", "C", "G"))])) == 1
    }
  )] %>%
    unique()

  t_and_s$tree_tibble$codons = map(
    t_and_s$tree_tibble$reconstructed_dna_sequence,
    ~substring(.x, seq(1, nchar(.x), 3), seq(3, nchar(.x), 3))
  )

  t_and_s$tree_tibble$four_fold_syn_codons = map(
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


getSplitFourFoldSynNucPositions = function(
    t_and_s,
    proportion_syn_threshold,
    identity_threshold
){
  four_fold_syn_nuc_positions = getFourFoldSynPositions(
    t_and_s,
    proportion_syn_threshold
  ) * 3

  four_fold_syn_nuc_position_aas = map(
    four_fold_syn_nuc_positions,
    function(pos){
      nts = substr(
        t_and_s$sequences$dna_sequence,
        pos, pos)
      nts = sort(table(nts), decreasing = T) /
        sum(nts %in% c("A", "T", "C", "G"))
      list(nt = names(nts)[[1]],
           prop = nts[[1]])
    }
  )

  four_fold_syn_nuc_positions = split(
    four_fold_syn_nuc_positions[
      map_dbl(four_fold_syn_nuc_position_aas, "prop") > identity_threshold],
    map_chr(four_fold_syn_nuc_position_aas, "nt")[
      map_dbl(four_fold_syn_nuc_position_aas, "prop") > identity_threshold]
  )

  four_fold_syn_nuc_positions
}

getNucCounts = function(
    t_and_s,
    four_fold_syn_nuc_positions = NULL
) {

  rates = tibble(
    from = character(),
    to = character(),
    from_to = character(),
    n_mutations_individual_positions = list(),
    n_mutations = integer(),
    mutations_per_site = numeric(),
    n_sites = integer(),
    n_seqs = integer(),
    n_nodes = integer()
  )

  fixations = t_and_s$tree_tibble %>%
    select(nt_mutations_syn) %>%
    unnest(nt_mutations_syn) %>%
    mutate(
      from = str_sub(nt_mutations_syn, 1, 1),
      at = str_sub(nt_mutations_syn, 2, -2) %>% as.integer(),
      to = factor(
        str_sub(nt_mutations_syn, -1, -1),
        levels = c("A", "T", "C", "G")
      )
    )

  if (missing(four_fold_syn_nuc_positions)){
    four_fold_syn_nuc_positions = getSplitFourFoldSynNucPositions(
      t_and_s,
      proportion_syn_threshold = 0.95,
      identity_threshold = 0.9
    )
  }

  n_s = sum(!is.na(t_and_s$tree_tibble$dna_sequence))
  n_n = nrow(t_and_s$tree_tibble)

  for (fr in names(four_fold_syn_nuc_positions)){

    rates_nt = fixations %>%
      filter(at %in% four_fold_syn_nuc_positions[[fr]]) %>%
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
        mutations_per_site = n_mutations / n_sites
      )

    rates = add_row(
      rates,
      from = rep(fr, 3),
      to = as.character(rates_nt$to),
      from_to = paste0(rep(fr, 3), as.character(rates_nt$to)),
      n_mutations_individual_positions =
        rates_nt$n_mutations_individual_positions,
      n_mutations = rates_nt$n_mutations,
      mutations_per_site = rates_nt$mutations_per_site,
      n_sites = rates_nt$n_sites,
      n_seqs = rep(n_s, 3),
      n_nodes = rep(n_n, 3)
    )
  }

  rates = rates %>%
    filter(from != to)

  rates
}


getTreeSizeRatios = function(overall_nuc_counts){
  overall_nuc_counts %>%
    transmute(
      from = from,
      to = to,
      from_to = from_to,
      n_sites = n_sites,
      ratio = n_sites/sum(n_sites)
    )
}

getTreeSize = function(nuc_counts,
                       tree_size_ratios) {

  size = 0
  for (nt in unique(nuc_counts$from_to)){
    size = size +
      nuc_counts$mutations_per_site[nuc_counts$from_to == nt] *
      tree_size_ratios$ratio[tree_size_ratios$from_to == nt]
  }

  list(point = size)
}

getNucRates = function(nuc_counts, tree_size){

  nuc_counts$mutations_per_site_rate = nuc_counts$mutations_per_site / tree_size[["point"]]
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




# ------------------------------------------------------------------------------

#'@export
addEstimatedNodeDates = function(tree_and_sequences, steps = 10000){

  myUtils::add_to_PATH(path = "/Users/samturner/miniforge3/envs/treebuild/bin/")

  dir = paste0("results/trees/", paste0(sample(LETTERS, 10), collapse = ""), "/")
  dir.create(dir)
  on.exit(unlink(dir, recursive = T))

  name_to_hash = setNames(
    openssl::md4(tree_and_sequences$tree$tip.label),
    tree_and_sequences$tree$tip.label
  )

  hash_to_name = setNames(
    names(name_to_hash),
    unname(name_to_hash)
  )

  tree_and_sequences$tree_tibble$label = name_to_hash[
    tree_and_sequences$tree_tibble$label]

  tree_and_sequences$tree$tip.label = name_to_hash[
    tree_and_sequences$tree$tip.label]

  tree_and_sequences$tree$node.label = 1:Nnode(tree_and_sequences$tree)
  tree_and_sequences$tree_tibble$label[
    Ntip(tree_and_sequences$tree) + (1:Nnode(tree_and_sequences$tree)) ] =
    1:Nnode(tree_and_sequences$tree)

  tree_and_sequences$tree_tibble$Collection_date_clean = tree_and_sequences$tree_tibble$Collection_date

  tree_and_sequences$tree_tibble$Collection_date_clean[
    tree_and_sequences$tree_tibble$Collection_date_clean < as_date("1968-01-01")
  ] = NA

  tree_and_sequences$tree_tibble$Collection_date_clean[
    year(tree_and_sequences$tree_tibble$Collection_date_clean) != tree_and_sequences$tree_tibble$year
  ] = NA

  ref = (tree_and_sequences$tree_tibble %>%
           filter(str_sub(Collection_date_clean, 5) != "-01-01") %>%
           arrange(Collection_date_clean) %>%
           pluck("label"))[1]

  tip_dates = tree_and_sequences$tree_tibble$Collection_date_clean[
    match(tree_and_sequences$tree$tip.label, tree_and_sequences$tree_tibble$label)
  ]

  # tip_dates[tip_dates < as_date("1967-01-01")] = NA

  tip_dates = as.character(tip_dates)

  tip_dates[map_lgl(substr(tip_dates, 6, 10) == "01-01", isTRUE)] = substr(
    tip_dates[map_lgl(substr(tip_dates, 6, 10) == "01-01", isTRUE)],
    1, 4
  )

  tip_dates[map_lgl(substr(tip_dates, 9, 10) == "01", isTRUE)] = substr(
    tip_dates[map_lgl(substr(tip_dates, 9, 10) == "01", isTRUE)],
    1, 7
  )


  branch_length_factor = 1e7
  tree_and_sequences$tree$edge.length = tree_and_sequences$tree$edge.length * branch_length_factor

  if (steps > 0){

    ape::write.tree(
      tree_and_sequences$tree,
      paste0(dir, "tree.nw")
    )

    cbind(
      strain = tree_and_sequences$tree$tip.label[1:Ntip(tree_and_sequences$tree)],
      date = as.character(tip_dates[1:Ntip(tree_and_sequences$tree)])
    ) %>%
      write.csv(
        paste0(dir, "tree_dates.csv")
      )

    system(
      paste0(
        "chronumental ",
        "--tree ", dir, "tree.nw ",
        "--tree_out ", dir, "tree_out_no_var.nw ",
        "--dates ", dir, "tree_dates.csv ",
        "--dates_out ", dir, "tree_dates_no_var.csv ",
        "--steps ", steps
      )
    )


    dates_no_var = read_delim(paste0(dir, "tree_dates_no_var.csv")) %>%
      {setNames(.$predicted_date, .$strain)}

    tree_and_sequences$tree_tibble$estimated_date_chronumental =
      dates_no_var[tree_and_sequences$tree_tibble$label]

    tree_and_sequences$tree_dated = ape::read.tree(paste0(dir, "tree_out_no_var.nw"))

  } else {
    tree_and_sequences$tree_tibble$estimated_date_chronumental = NA

    # tree_and_sequences$tree_tibble$estimated_date_var = NA
  }

  # nearest tip dates
  nearest_tips = castor::find_nearest_tips(
    tree = tree_and_sequences$tree,
    target_tips = which(!is.na(tip_dates))
  )

  nearest_desc_tips = castor::find_nearest_tips(
    tree = tree_and_sequences$tree,
    target_tips = which(!is.na(tip_dates)),
    only_descending_tips = T
  )

  node_dates_from_nearest_tips = c(
    tip_dates[nearest_tips$nearest_tip_per_tip],
    tip_dates[nearest_tips$nearest_tip_per_node]
  )

  node_dates_from_nearest_desc_tips = c(
    tip_dates[nearest_tips$nearest_tip_per_tip], # no desc from tip
    tip_dates[nearest_desc_tips$nearest_tip_per_node]
  )

  tree_and_sequences$tree_tibble$estimated_date_nearest =
    node_dates_from_nearest_tips

  tree_and_sequences$tree_tibble$estimated_date_nearest_desc =
    node_dates_from_nearest_desc_tips

  # tree_and_sequences$tree_tibble$estimated_date_nearest_except_self = tip_dates[
  # find_nearest_tips_except_self(tree_and_sequences$tree,which(!is.na(tip_dates)))]


  tree_and_sequences$tree_tibble$distance_from_root = castor::get_all_distances_to_root(
    tree_and_sequences$tree
  ) / branch_length_factor

  tree_and_sequences$tree_tibble$node_height = ape::node.height(
    tree_and_sequences$tree
  )

  tree_and_sequences$tree_tibble$label = hash_to_name[
    tree_and_sequences$tree_tibble$label]

  tree_and_sequences$tree$tip.label = hash_to_name[
    tree_and_sequences$tree$tip.label]

  list(tree_and_sequences, ref)
}

find_nearest_tips_except_self = function(tree, considered_tips = 1:Ntip(tree)){

  nearest_tips = c()

  for (t in 1:Ntip(tree)){
    if (t %% 1000 == 0) message(t)

    nearest_tips[[t]] = castor::find_nearest_tips(
      tree,
      target_tips = considered_tips[considered_tips != t]
    )$nearest_tip_per_tip[[t]]
  }

  c(
    nearest_tips,
    unlist(
      castor::find_nearest_tips(
        tree,
        target_tips = considered_tips
      )$nearest_tip_per_node
    )
  ) %>% unlist()
}

# ------------------------------------------------------------------------------

#'@export
getAASubstitutionRatios = function(
    t_and_s_filtered,
    overall_nuc_rates,
    tree_size_ratios,
    tree_size_fn,
    positions
){
  message("\tCodon table")
  codon_table = getCodonTable(
    t_and_s_filtered$tree_tibble,
    positions
  )

  codon_table_filtered = filter(codon_table, at %in% positions)
  # codon_table_filtered = filter(codon_table, at == 156)
  message("\tExpected counts")
  codon_table_with_expected_n = addExpectedMutationsToCodonTable(
    t_and_s_filtered$tree_tibble,
    codon_table_filtered,
    overall_nuc_rates,
    tree_size_ratios,
    tree_size_fn
  )


  message("\tRatios")
  mutation_table = addObservedMutationCounts(
    codon_table_with_expected_n,
    t_and_s_filtered$tree_tibble
  )

  mutation_table = summariseMutationTableToAAsAndSyns(mutation_table)

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

get_ffs = function(sequences_at_four_fold_syn_sites, tip_nodes){
  sequences_at_four_fold_syn_sites_idx = sequences_at_four_fold_syn_sites[
    sequences_at_four_fold_syn_sites$node %in% tip_nodes,
    -ncol(sequences_at_four_fold_syn_sites)
  ]

  sequences_at_four_fold_syn_sites_idx = setNames(
    apply(sequences_at_four_fold_syn_sites_idx, 2, \(x) DescTools::Mode(x)[[1]]),
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
    overall_nuc_rates,
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
    four_fold_syn_sites_idx = get_ffs(
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
      overall_nuc_rates
    )

    codon_table$tree_size[[idx]] = tree_size_idx$point

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
  acc = function(b){
    if (length(b) == 0) return(list(NULL))
    list(purrr::reduce(
      .x = b,
      .f = function(a,b){a+b},
      .init = rep(0, length(b[[1]]))
    ))
  }

  syn_nt_mutations = mutation_table %>%
    mutate(
      is_syn = (
        Biostrings::GENETIC_CODE[codon] ==
          Biostrings::GENETIC_CODE[
            pmap_chr(list(codon, at, to), function(codon,at,to){applySubstitutions(
              codon,
              paste0("X", 1+(unique(at)-1) %% 3, unique(to))
            )})]
      )) %>%
    filter(is_syn)

  if (nrow(syn_nt_mutations) > 0){

    syn_nt_mutations = syn_nt_mutations%>%
      group_by(nt_mutation, codon) %>%
      summarise(
        mutations_per_site_interp = sum(mutations_per_site_interp),
        n = sum(n),
        # n_MAP = n,
        tree_size = unique(tree_size),
        .groups = "drop"
      ) %>%
      ungroup() %>%
      mutate(
        ratio = n / mutations_per_site_interp,
      ) %>%
      arrange(-ratio) %>%
      mutate(
        from = str_sub(nt_mutation, 1, 1),
        at = str_sub(nt_mutation, 2, -2),
        to = str_sub(nt_mutation, -1, -1),
        expected_n = mutations_per_site_interp
      )
  }

  # browser()
  aa_mutation_table = mutation_table %>%
    distinct() %>%
    group_by(nt_mutation, codon) %>%
    summarise(
      mutations_per_site_interp = sum(mutations_per_site_interp),
      n = sum(n),
      aa_mutation = unique(aa_mutation),
      is_syn = Biostrings::GENETIC_CODE[codon] == Biostrings::GENETIC_CODE[
        applySubstitutions(
          codon,
          paste0("X", 1+(unique(at)-1) %% 3, unique(to))
        )],
      tree_size = unique(tree_size),
      .groups = "drop"
    ) %>%
    filter(!is_syn) %>%
    group_by(aa_mutation) %>%
    summarise(
      mutations_per_site_interp = sum(mutations_per_site_interp),
      n = sum(n),
      tree_size = sum(tree_size),
      .groups = "drop"
    ) %>%
    mutate(
      ratio = n / mutations_per_site_interp
    ) %>%
    arrange(-ratio) %>%
    mutate(
      aa_mutation = factor(aa_mutation, levels = unique(aa_mutation)),
      from = str_sub(aa_mutation, 1, 1),
      at = str_sub(aa_mutation, 2, -2),
      to = str_sub(aa_mutation, -1, -1),
      expected_n = mutations_per_site_interp
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
        cluster_tree_size[["point"]] *
        nuc_rates_overall$mutations_per_site_rate[[
          which(nuc_rates_overall$from == from &
                  nuc_rates_overall$to == to)]]

    }
  }


  nuc_rates_by_cluster
}
