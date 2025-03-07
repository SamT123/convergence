
#' Title
#'
#' @param tree_and_sequences tree_and_sequences list
#' @param tree_info tree_info from `getTreeSizeAndNucRates`
#' @param time_boundaries Boundaries between time intervals. The final interval will have boundaries [last_time_interval, Inf]
#' @param positions amino acid positions to analyse
#'
#' @export
getTimeIntervalSubstitutionRatios = function(
    tree_and_sequences,
    tree_info,
    time_boundaries,
    positions = 1:550
){

  time_boundaries = c(time_boundaries, "2100-01-01")

  ratio_list = list()
  for (i in seq_along(time_boundaries[-1])){

    tree_tibble_filtered = filter(
      tree_and_sequences$tree_tibble,
      as_date(
        estimated_date_nearest,
        format = c("%Y-%m-%d", "%Y-%m", "%Y")
      ) >= as_date(time_boundaries[[i]]),
      as_date(
        estimated_date_nearest,
        format = c("%Y-%m-%d", "%Y-%m", "%Y")
      ) < as_date(time_boundaries[[i+1]])
    )

    if (nrow(tree_tibble_filtered) == 0) next

    ratios = getSubstitutionRatios(
      list(
        tree_tibble = tree_tibble_filtered
      ),
      tree_info$nuc_rates,
      tree_info$tree_size_ratios,
      tree_info$tree_size_fn,
      tree_info$sampling_function,
      positions = positions
    )


    time_interval_name = paste0(
      as.character(time_boundaries[[i]]),
      " to ",
      ifelse(i == length(time_boundaries)-1, "Inf", as.character(time_boundaries[[i+1]]))
    )

    ratio_list[[time_interval_name]] = ratios
  }

  ratio_list
}

#'@export
getInfoAboutNodes = function(
    tree_and_sequences,
    nodes,
    comparison_sequence
) {

  recent_occurrences = tree_and_sequences$tree_tibble %>%
    filter(node %in% nodes)


  recent_occurrence_descs = map(
    recent_occurrences$node,
    ~c(.x, convergence:::getDescendants(tree_and_sequences$tree, .x))
  )

  recent_occurrence_descs_info = map(
    recent_occurrence_descs,
    function(desc_nodes){

      node_root_sequence_dna = filter(
        tree_and_sequences$tree_tibble,
        node == desc_nodes[[1]]
      )$reconstructed_dna_sequence

      node_root_sequence_aa = node_root_sequence_dna %>%
        seqUtils::translate(seqUtils::alaska_232_2015_aas)

      root_diff_to_comparison_comparison_sequence = seqUtils::get_substitutions(
        comparison_sequence,
        node_root_sequence_aa
      ) %>% paste(collapse = " ")


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
          syn_diff_from_root = map2(
            dna_diff_from_root,
            aa_diff_from_root,
            function(dna,aa){
              dna_at = ceiling(as.integer(str_sub(dna, 2, -2))/3)
              aa_at = as.integer(str_sub(aa, 2, -2))

              dna[!dna_at %in% aa_at]
            }
          ),

          aa_diff_from_root = aa_diff_from_root %>%
            map_chr(~paste(.x, collapse = " ")),

          syn_diff_from_root = syn_diff_from_root %>%
            map_chr(~paste(.x, collapse = " ")),

        ) %>%
        select(Isolate_name_standard, Collection_date, aa_diff_from_root, syn_diff_from_root)

      list(
        root_diff_to_comparator = root_diff_to_comparison_comparison_sequence,
        sequences = sequences
      )
    }
  )

  recent_occurrence_descs_info
}

#' @export
getConsensus = function(sequences){
  Biostrings::consensusMatrix(sequences) %>%
    {apply(., 2, \(x){rownames(.)[which.max(x)]})} %>%
    paste(collapse = "")
}
