#' @title Wrapper for CMAPLE
#'
#' @param sequences named character vector of DNA/RNA sequences
#' @param tree_path path to save tree
#' @param multifurcating is the tree allowed to contain multifucations? default TRUE.
#' @param starting_tree_path path to starting tree containing a subset of the sequences
#' @param freeze_starting_tree should the relationships in the starting tree be preserved?
#' @param cmaple_path path to your cmaple installation
#' @param out_sequence name of outgroup for tree rooting
#'
#'@export
make_cmaple_tree = function(
  sequences,
  tree_path,
  multifurcating = T,
  starting_tree_path = NULL,
  freeze_starting_tree = F,
  cmaple_path = NULL,
  out_sequence = NULL
) {
  if (!is.null(cmaple_path)) {
    add_to_PATH(cmaple_path)
  }

  if (!is.null(starting_tree_path)) {
    starting_tree = castor::read_tree(file = starting_tree_path)
    message(
      sum(names(sequences) %in% starting_tree$tip.label),
      " / ",
      length(sequences),
      " sequences already in starting tree"
    )

    message(
      sum(starting_tree$tip.label %in% names(sequences)),
      " / ",
      length(starting_tree$tip.label),
      " starting tree tips in sequence list"
    )
  }

  fasta_path = fs::path_ext_set(tree_path, ".fasta")

  seqUtils::write_fast_fasta(
    sequences,
    names(sequences),
    path = fasta_path
  )

  # fmt: skip
  cmaple_call = c(
    "cmaple",
    "-aln", fasta_path
  )

  if (multifurcating) {
    # fmt: skip
    cmaple_call = c(
      cmaple_call,
      "--out-mul-tree"
    )
  }

  if (!is.null(starting_tree_path)) {
    # fmt: skip
    cmaple_call = c(
      cmaple_call,
      "-t ", starting_tree_path
    )
  }

  if (!freeze_starting_tree) {
    # fmt: skip
    cmaple_call = c(
      cmaple_call,
      "-search", "EXHAUSTIVE"
    )
  } else {
    # fmt: skip
    cmaple_call = c(
      cmaple_call,
      "-search", "NORMAL"
    )
  }

  # fmt: skip
  cmaple_call = c(
      cmaple_call,
      "-nt", "AUTO",
      "-m", "GTR"
    )

  system(paste(cmaple_call, collapse = " "))

  system(paste0("mv ", fasta_path, ".treefile ", tree_path))
  system(paste0("mv ", fasta_path, ".log ", tree_path, ".log"))

  tree = castor::read_tree(file = tree_path)
  ladderizeAndMaybeRoot(tree, out_sequence)
  castor::write_tree(tree, tree_path)

  tree_path
}


ladderizeAndMaybeRoot = function(tree, out_sequence = NULL) {
  if (!is.null(out_sequence)) {
    stopifnot(out_sequence %in% tree$tip.label)

    tree = tree %>%
      ape::unroot() %>%
      ape::root(outgroup = out_sequence)
  }

  tree = ape::ladderize(tree)

  f = fs::file_temp()
  castor::write_tree(tree, file = f)
  tree_read = castor::read_tree(file = f)

  tree_read
}
