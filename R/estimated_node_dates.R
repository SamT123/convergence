#' Add estiamted node dates to tree_tibble
#'
#' @param tree_and_sequences a list containing `tree` and `tree_tibble`
#'
#' Two types of node date are added: the collection date of the nearest tip, and the collection date of the nearest _descendant_ tip.
#'
#' To do:
#' - Node dating with TimeTree and/or Chronumental
#' - Earliest date of any descendant tip (but this is very sensitive to misplaced sequences)
#'
#'@export
addEstimatedNodeDates = function(tree_and_sequences){

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

  tree_and_sequences$tree$node.label = 1:ape::Nnode(tree_and_sequences$tree)
  tree_and_sequences$tree_tibble$label[
    ape::Ntip(tree_and_sequences$tree) + (1:ape::Nnode(tree_and_sequences$tree)) ] =
    1:ape::Nnode(tree_and_sequences$tree)


  tree_and_sequences$tree_tibble$Collection_date_clean = tree_and_sequences$tree_tibble$Collection_date

  tree_and_sequences$tree_tibble$Collection_date_clean[
    tree_and_sequences$tree_tibble$Collection_date_clean < as_date("1968-01-01")
  ] = NA

  tree_and_sequences$tree_tibble$Collection_date_clean[
    year(tree_and_sequences$tree_tibble$Collection_date_clean) != tree_and_sequences$tree_tibble$year
  ] = NA


  tip_dates = tree_and_sequences$tree_tibble$Collection_date_clean[
    match(tree_and_sequences$tree$tip.label, tree_and_sequences$tree_tibble$label)
  ]

  tip_dates = as.character(tip_dates)

  tip_dates[map_lgl(substr(tip_dates, 6, 10) == "01-01", isTRUE)] = substr(
    tip_dates[map_lgl(substr(tip_dates, 6, 10) == "01-01", isTRUE)],
    1, 4
  )

  tip_dates[map_lgl(substr(tip_dates, 9, 10) == "01", isTRUE)] = substr(
    tip_dates[map_lgl(substr(tip_dates, 9, 10) == "01", isTRUE)],
    1, 7
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

  tree_and_sequences$tree_tibble$label = hash_to_name[
    tree_and_sequences$tree_tibble$label]

  tree_and_sequences$tree$tip.label = hash_to_name[
    tree_and_sequences$tree$tip.label]

  tree_and_sequences
}
