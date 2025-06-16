#' @title Wrapper for CMAPLE
#'
#' @param sequences named character vector of DNA/RNA sequences
#' @param tree_path path to save tree
#' @param starting_tree_path path to starting tree containing a subset of the sequences
#' @param freeze_starting_tree should the relationships in the starting tree be preserved?
#'
#'@export
make_cmaple_tree = function(
  sequences,
  tree_path,
  starting_tree_path = NULL,
  freeze_starting_tree = F
) {
  fasta_path = paste0(
    tools::file_path_sans_ext(tree_path),
    ".fasta"
  )

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

  old_path <- Sys.getenv("PATH")
  Sys.setenv(
    PATH = paste(old_path, "~/miniforge3/envs/treebuild/bin/", sep = ":")
  )

  seqUtils::write_fast_fasta(
    sequences,
    names(sequences),
    path = fasta_path
  )

  if (!is.null(starting_tree_path)) {
    if (!freeze_starting_tree) {
      system(
        paste0(
          "cmaple -aln ",
          fasta_path,
          " -t ",
          starting_tree_path,
          " -search EXHAUSTIVE -nt AUTO -m GTR"
        )
      )
    } else {
      system(
        paste0(
          "cmaple -aln ",
          fasta_path,
          " -t ",
          starting_tree_path,
          " -nt AUTO -m GTR"
        )
      )
    }
  } else {
    system(
      paste0(
        "cmaple -aln ",
        fasta_path,
        " -search EXHAUSTIVE -nt AUTO -m GTR"
      )
    )
  }

  system(paste0("mv ", fasta_path, ".treefile ", tree_path))
  system(paste0("mv ", fasta_path, ".log ", tree_path, ".log"))

  tree_path
}
