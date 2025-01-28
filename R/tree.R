make_fn_safe = function(txt){
  str_replace_all(txt, "[ /]", "_")
}

#'@export
make_cmaple_tree = function(sequences, tree_path, starting_tree_path = NULL, freeze_starting_tree = F){

  fasta_path = paste0(
    tools::file_path_sans_ext(tree_path),
    ".fasta"
  )

  if (!is.null(starting_tree_path)){
    starting_tree = castor::read_tree(file = starting_tree_path)
    message(
      sum(names(sequences) %in% starting_tree$tip.label), " / ", length(sequences),
      " sequences already in starting tree"
    )

    message(
      sum(starting_tree$tip.label %in% names(sequences)), " / ", length(starting_tree$tip.label),
      " starting tree tips in sequence list"
    )
  }

  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(old_path, "~/miniforge3/envs/treebuild/bin/", sep = ":"))

  myUtils::fast_fasta.write(
    sequences,
    names(sequences),
    path = fasta_path
  )

  if (!is.null(starting_tree_path)){

    if (!freeze_starting_tree){
      system(
        paste0(
          "cmaple -aln ", fasta_path, " -t ", starting_tree_path, " -search EXHAUSTIVE -nt AUTO -m GTR"
        )
      )
    } else {
      system(
        paste0(
          "cmaple -aln ", fasta_path, " -t ", starting_tree_path, " -nt AUTO -m GTR"
        )
      )
    }
  } else {
    system(
      paste0(
        "cmaple -aln ", fasta_path, " -search EXHAUSTIVE -nt AUTO -m GTR"
      )
    )
  }

  system(paste0("mv ", fasta_path, ".treefile ", tree_path))
  system(paste0("mv ", fasta_path, ".log ", tree_path, ".log"))

  tree_path
}




#' Run IQ-TREE
#'
#' For details, see http://www.iqtree.org/doc/
#'
#' @param alignment DNA alignment to use for phylogenetic analysis. Must be matrix
#' (i.e., aligned sequences) of class DNAbin
#' @param wd Path to working directory. The alignment and IQ-TREE intermediate files
#' and results will be written here.
#' @param bb Optional; number of ultrafast bootstrap replicates to run.
#' @param nt Optional; number of cores to use. Set to "AUTO" to determine automatically.
#' @param alrt Optional; number of SH-aLRT tests to run.
#' @param m Optional; specify model. If no model is given, ModelTest will be run
#' to identify the best model for the data.
#' @param redo Logical; should the analysis be redone from scratch if output from
#' previous runs is present?
#' @param echo Logical; should STDERR be written to the screen?
#' @param ... Other arguments not used by this function but used by
#' drake for tracking.
#'
#'
#' @return Phylogenetic tree (list of class "phylo")
#'
#' @examples
#' data(woodmouse)
#' # Rapid boot-strap tree with 1000 replicates on best-fitting model
#' tree <- iq_tree(woodmouse, tempdir(), bb = 1000, echo = TRUE)
#' plot(tree)
#' # Check the optimum number of cores to use for GTR+I+G model
#' iq_tree(tempdir(), woodmouse, m = "GTR+I+G", nt = "AUTO", echo = TRUE, redo = TRUE)
make_iq_tree <- function(alignment,
                          wd,
                          bb = NULL,
                          nt = NULL,
                          alrt = NULL,
                          m = NULL,
                          redo = FALSE,
                          echo = FALSE,
                          ...) {

  assertthat::assert_that(inherits(alignment, "DNAbin"),
                          msg = "alignment must be of class 'DNAbin'")
  assertthat::assert_that(is.matrix(alignment),
                          msg = "alignment must be a matrix (not a list of unaligned sequences)")

  assertthat::assert_that(assertthat::is.dir(wd))

  assertthat::assert_that(is.logical(echo))

  assertthat::assert_that(is.logical(redo))

  if(!is.null(bb))
    assertthat::assert_that(assertthat::is.number(bb))

  if(!is.null(alrt))
    assertthat::assert_that(assertthat::is.number(alrt))

  if(!is.null(nt))
    assertthat::assert_that(assertthat::is.number(nt) | assertthat::is.string(nt))

  if(!is.null(m))
    assertthat::assert_that(assertthat::is.string(m))

  wd <- fs::path_norm(wd)

  # check that iqtree is installed and on the PATH
  tryCatch({
    processx::run("iqtree2", "-h", echo = FALSE)
  }, warning = function(w) {
    stop("iqtree not installed and on path")
  }, error = function(e) {
    stop("iqtree not installed and on path")
  }, finally = {
    TRUE
  })

  # Write alignment to working directory in phylip format
  aln_path <- "aln.fa"

  ape::write.FASTA(x = alignment, file = paste0(wd, "/", aln_path))

  # Set up arguments
  iqtree_arguments <- c(
    "-s", aln_path,
    if(!is.null(bb)) "-bb",
    bb,
    if(!is.null(alrt)) "-alrt",
    alrt,
    if(!is.null(nt)) "-nt",
    nt,
    if(!is.null(m)) "-m",
    m,
    if(isTRUE(redo)) "-redo"
  )

  iqtree_arguments_p = paste(iqtree_arguments, collapse = " ")
  print(iqtree_arguments_p)

  with_dir(
    wd,
    code = {system(paste0("iqtree2 ", iqtree_arguments_p))}
  )


  # Run iqtree command
  # processx::run(
  #   "iqtree2",
  #   iqtree_arguments, wd = wd, echo = echo)

  # Read in resulting consensus tree
  tree_path <- paste0(wd, "/", aln_path, ".treefile")

  ape::read.tree(tree_path)

}
