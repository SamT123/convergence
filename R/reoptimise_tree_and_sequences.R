#' Reoptimise the tree in a tree_and_sequences object
#'
#' @param tree_and_sequences tree_and_sequences list
#' @param f function that takes a named vector of DNA sequences and returns a phylo object
#'
#' @examples
#' \dontrun{
#'   reoptimiseTreeAndSequences(
#'     tree_and_sequences,
#'     function(seqs) {
#'       t = seqUtils::make_cmaple_tree(
#'         sequences = c(seqs, outgroup = seqUtils::alaska_232_2015_nts),
#'         tree_path = "tmp.nwk",
#'         cmaple_path = "/path/to/bin",
#'         keep_files = c(),
#'         out_sequence = "outgroup"
#'       )
#'       ape::drop.tip(t, "outgroup")
#'     }
#'   )
#' }
#'
#' @export
reoptimiseTreeAndSequences = function(tree_and_sequences, f) {
  stopifnot(is.function(f))

  sequence_vector = setNames(
    tree_and_sequences$sequences$dna_sequence,
    tree_and_sequences$sequences$Isolate_unique_identifier
  )

  stopifnot(
    identical(
      unname(sort(names(sequence_vector))),
      unname(sort(tree_and_sequences$tree$tip.label))
    )
  )

  new_tree = f(sequence_vector)

  # Validate the new tree
  stopifnot("phylo" %in% class(new_tree))

  stopifnot(
    identical(
      unname(sort(new_tree$tip.label)),
      unname(sort(names(sequence_vector)))
    )
  )

  makeTreeAndSequences(new_tree, tree_and_sequences$sequences)
}
