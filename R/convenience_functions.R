#' Title
#'
#' @param tree_and_sequences tree_and_sequences list
#' @param tree_info tree_info from `getTreeSizeAndNucRates`
#' @param time_boundaries Boundaries between time intervals. The final interval will have boundaries (last_time_interval, Inf)
#' @param positions amino acid positions to analyse
#' @param date_column column to use for date filtering
#'
#' @export
getTimeIntervalSubstitutionRatios = function(
  tree_and_sequences,
  tree_info,
  time_boundaries,
  positions = 1:550,
  date_column = "estimated_date_nearest"
) {
  time_boundaries = c(time_boundaries, "2100-01-01")

  stopifnot(date_column %in% colnames(tree_and_sequences$tree_tibble))
  tree_and_sequences$tree_tibble[["DATECOLUMN"]] =
    tree_and_sequences$tree_tibble[[date_column]]

  ratio_list = list()
  for (i in seq_along(time_boundaries[-1])) {
    tree_tibble_filtered = filter(
      tree_and_sequences$tree_tibble,
      lubridate::as_date(
        DATECOLUMN,
        format = c("%Y-%m-%d", "%Y-%m", "%Y")
      ) >=
        lubridate::as_date(time_boundaries[[i]]),
      lubridate::as_date(
        DATECOLUMN,
        format = c("%Y-%m-%d", "%Y-%m", "%Y")
      ) <
        lubridate::as_date(time_boundaries[[i + 1]])
    )

    if (nrow(tree_tibble_filtered) == 0) {
      next
    }

    ratios = getSubstitutionRatios(
      list(
        tree_tibble = tree_tibble_filtered
      ),
      tree_info$nuc_rates,
      tree_info$tree_size_fn,
      tree_info$sampling_function,
      positions = positions
    )

    time_interval_name = paste0(
      as.character(time_boundaries[[i]]),
      " to ",
      ifelse(
        i == length(time_boundaries) - 1,
        "Inf",
        as.character(time_boundaries[[i + 1]])
      )
    )

    ratio_list[[time_interval_name]] = ratios
  }

  ratio_list
}

#' Helper function to get information about nodes in the tree
#'
#' Useful for learning more about occurrences of substitutions of interest - e.g. genetic context
#'
#' @param tree_and_sequences ...
#' @param nodes nodes to get info about
#' @param comparison_sequence for genetic context information
#'
#'@export
getInfoAboutNodes = function(
  tree_and_sequences,
  nodes,
  comparison_sequence
) {
  recent_occurrences = tree_and_sequences$tree_tibble %>%
    filter(node %in% nodes)

  recent_occurrence_descs = purrr::map(
    recent_occurrences$node,
    ~ c(.x, fastGetDescendants(tree_and_sequences$tree, .x))
  )

  recent_occurrence_descs_info = purrr::map(
    recent_occurrence_descs,
    function(desc_nodes) {
      node_root_sequence_dna = filter(
        tree_and_sequences$tree_tibble,
        node == desc_nodes[[1]]
      )$reconstructed_dna_sequence

      node_root_sequence_aa = node_root_sequence_dna %>%
        seqUtils::translate(seqUtils::alaska_232_2015_aas)

      root_diff_to_comparison_comparison_sequence = seqUtils::get_substitutions(
        comparison_sequence,
        node_root_sequence_aa
      ) %>%
        paste(collapse = " ")

      sequences = filter(
        tree_and_sequences$tree_tibble,
        node %in% desc_nodes
      ) %>%
        filter(!is.na(dna_sequence)) %>%
        mutate(
          aa_diff_from_root = seqUtils::get_substitutions(
            node_root_sequence_aa,
            aa_sequence,
            simplify = F
          ),

          dna_diff_from_root = seqUtils::get_substitutions(
            node_root_sequence_dna,
            dna_sequence,
            simplify = F
          ),
          syn_diff_from_root = purrr::map2(
            dna_diff_from_root,
            aa_diff_from_root,
            function(dna, aa) {
              dna_at = ceiling(as.integer(stringr::str_sub(dna, 2, -2)) / 3)
              aa_at = as.integer(stringr::str_sub(aa, 2, -2))

              dna[!dna_at %in% aa_at]
            }
          ),

          aa_diff_from_root = aa_diff_from_root %>%
            map_chr(~ paste(.x, collapse = " ")),

          syn_diff_from_root = syn_diff_from_root %>%
            map_chr(~ paste(.x, collapse = " ")),
        ) %>%
        select(
          Isolate_name_standard,
          Collection_date,
          aa_diff_from_root,
          syn_diff_from_root
        )

      list(
        root_diff_to_comparator = root_diff_to_comparison_comparison_sequence,
        sequences = sequences
      )
    }
  )

  recent_occurrence_descs_info
}

#' Get the consensus sequence from a vector of sequences
#' @param sequences a vector of sequences
#'
#' @export
getConsensus = function(sequences) {
  Biostrings::consensusMatrix(sequences) %>%
    {
      apply(., 2, \(x) {
        rownames(.)[which.max(x)]
      })
    } %>%
    paste(collapse = "")
}

#' Find branch with mutations
#'
#' Finds the branch where all of the target_aa_substitutions and target_nt_mutations are present, and at leaast occurs on the branch. Errors if the number of such branches is not exactly one.
#'
#' @param tree_tibble tree_tibble
#' @param target_aa_substitutions aa substitutions which the branch should have
#' @param target_nt_mutations nt mutations which the branch should have
#'
#' @export
findBranch = function(
  tree_tibble,
  target_aa_substitutions,
  target_nt_mutations
) {
  possible_ancestors = tree_tibble %>%
    as_tibble()

  # must have specified aa substitutions
  for (s in target_aa_substitutions) {
    at = as.integer(stringr::str_sub(s, 1, -2))
    to = stringr::str_sub(s, -1, -1)
    possible_ancestors = filter(
      possible_ancestors,
      stringr::str_sub(reconstructed_aa_sequence, at, at) == to
    )
  }

  # must have specified nt mutations
  for (t in target_nt_mutations) {
    at = as.integer(stringr::str_sub(t, 1, -2))
    to = stringr::str_sub(t, -1, -1)
    possible_ancestors = filter(
      possible_ancestors,
      stringr::str_sub(reconstructed_dna_sequence, at, at) == to
    )
  }

  # one of the specified nt mutations or aa substitutions must have occurred on the branch

  possible_ancestors = filter(
    possible_ancestors,
    purrr::map_lgl(
      aa_mutations_nonsyn,
      ~ any(target_aa_substitutions %in% stringr::str_sub(.x, 2))
    ) |
      purrr::map_lgl(
        nt_mutations,
        ~ any(target_nt_mutations %in% stringr::str_sub(.x, 2))
      )
  )

  if (nrow(possible_ancestors) == 0) {
    stop("no branches found")
  } else if (nrow(possible_ancestors) > 1) {
    stop("more than one branch found")
  }

  possible_ancestors$node[[1]]
}

#' Trim tree and sequences to descendants of a branch
#'
#' @param tree_and_sequences tree_and_sequences list
#' @param branch branch to keep descendants of
#'
#' @export
trimTreeAndSequences = function(
  tree_and_sequences,
  branch
) {
  treedata = tidytree::as.treedata(
    tree_and_sequences$tree_tibble,
    label = "label"
  )

  subset_treedata = tidytree::keep.tip(
    treedata,
    treeio::offspring(
      treedata,
      branch,
      type = "tips"
    )
  )

  subset_tree_and_sequences = list(
    tree = subset_treedata@phylo,
    tree_tibble = tidytree::as_tibble(subset_treedata),
    sequences = filter(
      tree_and_sequences$sequences,
      Isolate_unique_identifier %in% subset_treedata@phylo$tip.label
    )
  )

  lad = ladderizeTreeAndTib(
    subset_tree_and_sequences$tree,
    subset_tree_and_sequences$tree_tibble
  )

  list(
    tree = lad$tree,
    tree_tibble = lad$tib,
    sequences = subset_tree_and_sequences$sequences
  )
}
