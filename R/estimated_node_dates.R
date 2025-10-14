#' Add estimated node dates to tree_tibble
#'
#' @param tree_and_sequences a list containing `tree` and `tree_tibble`
#'
#' Two types of node date are added: the collection date of the nearest tip, and the collection date of the nearest _descendant_ tip.
#'@export
addNearestDescendantNodeDates = function(tree_and_sequences) {
  # clean collection dates
  tree_and_sequences$tree_tibble$Collection_date_clean = tree_and_sequences$tree_tibble$Collection_date

  tree_and_sequences$tree_tibble$Collection_date_clean[
    tree_and_sequences$tree_tibble$Collection_date_clean <
      lubridate::as_date("1968-01-01")
  ] = NA

  tree_and_sequences$tree_tibble$Collection_date_clean[
    year(tree_and_sequences$tree_tibble$Collection_date_clean) !=
      tree_and_sequences$tree_tibble$year
  ] = NA

  tip_dates = tree_and_sequences$tree_tibble$Collection_date_clean[
    match(
      tree_and_sequences$tree$tip.label,
      tree_and_sequences$tree_tibble$label
    )
  ]

  tip_dates = as.character(tip_dates)

  tip_dates[purrr::map_lgl(
    substr(tip_dates, 6, 10) == "01-01",
    isTRUE
  )] = substr(
    tip_dates[purrr::map_lgl(substr(tip_dates, 6, 10) == "01-01", isTRUE)],
    1,
    4
  )

  tip_dates[purrr::map_lgl(substr(tip_dates, 9, 10) == "01", isTRUE)] = substr(
    tip_dates[purrr::map_lgl(substr(tip_dates, 9, 10) == "01", isTRUE)],
    1,
    7
  )

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

  tree_and_sequences
}

#' Convert tree_and_sequences to a timetree using Chronumental
#'
#' @param tree_and_sequences a list containing `tree` and `tree_tibble`
#' @param reference_strain name of the strain to use as date reference (t=0). NA -> earliest date
#' @param chronumental_path path to chronumental binary
#' @param n_steps number of chronumental iterations
#' @param genome_size nt sequence length when branch lengths are in units of substitutions/site
#'
#' Two types of node date are added: the collection date of the nearest tip, and the collection date of the nearest _descendant_ tip.
#'@export
toChronumentalTree = function(
  tree_and_sequences,
  reference_strain = NA,
  chronumental_path = NULL,
  n_steps = 10000,
  genome_size = NA
) {
  if ("nt_mutations" %in% colnames(tree_and_sequences$tree_tibble)) {
    stop(
      "toChronumentalTree should be run BEFORE addASRusher, as toChronumentalTree changes the tree rooting"
    )
  }

  tree = treeio::as.phylo(
    tree_and_sequences$tree_tibble,
    label = "label",
    branch.length = "branch.length"
  )

  dates = tree_and_sequences$sequences %>%
    select(Isolate_unique_identifier, Collection_date) %>%
    rename(strain = Isolate_unique_identifier, date = Collection_date)

  if (is.na(reference_strain)) {
    reference_strain = dates %>%
      arrange(dates) %>%
      slice(1) %>%
      pluck("Isolate_unique_identifier")
  }

  tree_and_dates = makeChronumentalTree(
    tree,
    dates,
    reference_strain,
    chronumental_path,
    n_steps,
    genome_size
  )

  tree_and_sequences$original_tree = tree_and_sequences$tree
  tree_and_sequences$original_tree_tibble = tree_and_sequences$tree_tibble

  tree_and_sequences$tree = ape::ladderize(tree_and_dates$tree)
  tree_and_sequences$tree_tibble = tree_and_sequences$tree %>%
    treeio::as_tibble() %>%
    left_join(
      select(
        tree_and_sequences$original_tree_tibble[
          seq_len(ape::Ntip(tree_and_sequences$original_tree)),
        ],
        -parent,
        -node,
        -branch.length
      ),
      by = "label"
    ) %>%
    left_join(
      tree_and_dates$dates,
      by = c("label" = "strain")
    )

  tree_and_sequences
}


makeChronumentalTree = function(
  tree,
  dates,
  reference_strain,
  chronumental_path,
  n_steps = 10000,
  genome_size
) {
  if (!is.null(chronumental_path)) {
    convergence::addASRusher(chronumental_path)
  }

  ### write tree --------------------
  tree_file = fs::file_temp(ext = ".nwk")
  castor::write_tree(
    tree,
    file = tree_file
  )

  ### write dates --------------------
  dates_file = fs::file_temp(ext = ".csv")
  readr::write_csv(
    dates,
    file = dates_file
  )

  ### outfiles --------------------
  tree_out_file = fs::file_temp(ext = ".nwk")
  dates_out_file = fs::file_temp(ext = ".csv")

  # fmt: skip
  chronumental_call = c(
    "chronumental",
    "--tree", tree_file,
    "--dates", dates_file,
    "--tree_out", tree_out_file,
    "--dates_out", dates_out_file,
    "--reference_node", paste0("'", reference_strain, "'"),
    "--steps", n_steps
)

  if (!is.na(genome_size)) {
    chronumental_call = c(
      chronumental_call,
      "--treat_mutation_units_as_normalised_to_genome_size",
      genome_size
    )
  }

  system(paste(chronumental_call, collapse = " "))

  tree_out = castor::read_tree(file = tree_out_file)
  dates_out = readr::read_tsv(dates_out_file)
  dates_out$predicted_date = lubridate::as_date(dates_out$predicted_date)

  list(
    tree = tree_out,
    dates = dates_out
  )
}
