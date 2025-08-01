#' Add ancestral state reconstruction to `tree_and_sequences`
#' @param tree_and_sequences list with a sequence dataframe (see details) and a phylogenetic tree
#'
#' @param aa_ref amino acid reference sequence for alignment
#' @param nuc_ref nucleotide reference sequence for alignment
#'
#' @details The sequence dataframe in `tree_and_sequences` must contain columns `Isolate_unique_identifier` and `dna_sequence`. Other columns containing inforamtion about the sequences can also be present, and wil be preserved. Tip names in the tree in `tree_and_sequences` should correspond to names in the `Isolate_unique_identifier` column.
#'
#'@export
addASRusher = function(tree_and_sequences, aa_ref, nuc_ref) {
  tree_and_sequences = remakeTreeAndSequencesWithUsherASR(
    tree_and_sequences,
    aa_ref = aa_ref,
    nuc_ref = nuc_ref
  )

  tree_and_sequences_ladderized = ladderizeTreeAndTib(
    tree_and_sequences$tree,
    tree_and_sequences$tree_tibble
  )

  tree_and_sequences$tree = tree_and_sequences_ladderized$tree
  tree_and_sequences$tree_tibble = tree_and_sequences_ladderized$tib

  tree_and_sequences = addAaSequence(tree_and_sequences, aa_ref)
  tree_and_sequences = addSynonymousInfo(tree_and_sequences)

  tree_and_sequences
}

addSynonymousInfo = function(tree_and_sequences) {
  tree_and_sequences$tree_tibble$aa_mutations_syn = purrr::map2(
    tree_and_sequences$tree_tibble$aa_mutations,
    tree_and_sequences$tree_tibble$nt_mutations_is_syn,
    ~ .x[.y]
  )

  tree_and_sequences$tree_tibble$aa_mutations_nonsyn = purrr::map2(
    tree_and_sequences$tree_tibble$aa_mutations,
    tree_and_sequences$tree_tibble$nt_mutations_is_syn,
    ~ .x[!.y]
  )

  tree_and_sequences$tree_tibble$nt_mutations_syn = purrr::map2(
    tree_and_sequences$tree_tibble$nt_mutations,
    tree_and_sequences$tree_tibble$nt_mutations_is_syn,
    ~ .x[.y]
  )

  tree_and_sequences$tree_tibble$nt_mutations_nonsyn = purrr::map2(
    tree_and_sequences$tree_tibble$nt_mutations,
    tree_and_sequences$tree_tibble$nt_mutations_is_syn,
    ~ .x[!.y]
  )

  tree_and_sequences
}

addAaSequence = function(tree_and_sequences, aa_ref) {
  tree_and_sequences$tree_tibble$aa_sequence = NA
  tree_and_sequences$tree_tibble$aa_sequence[
    1:ape::Ntip(tree_and_sequences$tree)
  ] =
    suppressWarnings(
      seqUtils::translate(
        tree_and_sequences$tree_tibble$dna_sequence[
          1:ape::Ntip(tree_and_sequences$tree)
        ],
        reference_aas = substr(
          aa_ref,
          1,
          nchar(tree_and_sequences$tree_tibble$dna_sequence[[1]]) / 3
        )
      )
    )

  tree_and_sequences
}

do_ASR_usher = function(
  tree,
  sequences,
  nuc_ref,
  usher_path = NULL
) {
  if (!is.null(usher_path)) {
    add_to_PATH(usher_path)
  }

  dir = paste0(
    "tempdir",
    paste(sample(LETTERS, 10, replace = T), collapse = ""),
    "/"
  )
  dir.create(dir)
  on.exit(unlink(dir, recursive = T))

  castor::write_tree(
    tree = tree,
    file = paste0(dir, "/tree.nw")
  )

  seqUtils::write_fast_fasta(
    sequences,
    names(sequences),
    paste0(dir, "/sequences.fasta")
  )

  seqUtils::write_fast_fasta(
    nuc_ref[[1]],
    "ref",
    paste0(dir, "/nuc_ref.fasta")
  )

  system(
    paste0(
      "faToVcf ",
      paste0(dir, "/sequences.fasta"),
      " ",
      paste0(dir, "/sequences.vcf")
    )
  )

  system(
    paste0(
      "usher -t ",
      paste0(dir, "/tree.nw"),
      " -v ",
      paste0(dir, "/sequences.vcf"),
      " -o ",
      paste0(dir, "/MAT.pb"),
      " -l "
    )
  )

  cat(
    paste0(
      'NC_045512v2	ncbiGenes.genePred	CDS	1	',
      nchar(sequences[[1]]),
      '	.	+	0	gene_id "HA"; transcript_id "HA"; exon_number "1"; exon_id "HA";'
    ),
    file = paste0(dir, "/GTF.gtf")
  )

  system(
    paste0(
      "matUtils summary --translate ",
      paste0(dir, "/out.tsv"),
      " --input-mat ",
      paste0(dir, "/MAT.pb"),
      "  --input-gtf ",
      paste0(dir, "/GTF.gtf"),
      " --input-fasta ",
      paste0(dir, "/nuc_ref.fasta")
    )
  )

  system(
    paste0(
      "matUtils extract --input-mat ",
      paste0(dir, "/MAT.pb"),
      " --write-tree ",
      paste0(dir, "/new_tree.nw")
    )
  )

  list(
    tree = castor::read_tree(file = paste0(dir, "/new_tree.nw")),
    asr = readr::read_tsv(paste0(dir, "/out.tsv"))
  )
}

add_to_PATH = function(path) {
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(old_path, path, sep = ":"))
}

ladderizeTreeAndTib = function(tree, tib) {
  recover_tip_labels = setNames(
    tree$tip.label,
    paste0("t", seq_along(tree$tip.label))
  )

  node_labels_present = !is.null(tree$node.labels)
  if (node_labels_present) {
    recover_node_labels = setNames(
      tree$node.label,
      paste0("n", seq_len(ape::Nnode(tree)))
    )
  }

  tree$tip.label = paste0(
    "t",
    seq_along(tree$tip.label)
  )

  tree$node.label = paste0("n", seq_len(ape::Nnode(tree)))

  tib$nms = c(
    tree$tip.label,
    tree$node.label
  )

  tree = ape::ladderize(tree)

  tib = tib[match(c(tree$tip.label, tree$node.label), tib$nms), ]

  old2new = setNames(tib$node, stringr::str_sub(tib$nms, 2))
  tib$node = old2new[tib$node]
  tib$parent = old2new[tib$parent]
  tree$tip.label = recover_tip_labels[tree$tip.label]
  if (node_labels_present) {
    tree$node.label = recover_node_labels[tree$node.label]
  } else {
    tree$node.label = NULL
  }

  list(tree = tree, tib = tib)
}

applySubstitutions = function(seq, subs, backwards = F) {
  for (s in subs) {
    at = as.integer(stringr::str_sub(s, 2, -2))
    if (!backwards) {
      substr(seq, at, at) = stringr::str_sub(s, -1, -1)
    }
    if (backwards) substr(seq, at, at) = stringr::str_sub(s, 1, 1)
  }

  seq
}

getParents = function(tree_tibble, node) {
  while (tree_tibble[["parent"]][[node[[1]]]] != node[[1]]) {
    node = c(tree_tibble[["parent"]][[node[[1]]]], node)
  }
  node
}

getParentsFromNodelist = function(parents, node) {
  while (parents[[node[[1]]]] != node[[1]]) {
    node = c(parents[[node[[1]]]], node)
  }
  node
}

fastGetDescendants = function(phy, node) {
  if (node <= ape::Ntip(phy)) {
    return(NULL)
  }

  st = castor::get_subtree_at_node(
    phy,
    node = node - ape::Ntip(phy)
  )

  o = c(st$new2old_tip, st$new2old_node + ape::Ntip(phy))
  o[o != node]
}


reconstructNodeSequences = function(tree_tibble) {
  tree_tibble$nt_mutations[tree_tibble$parent == tree_tibble$node] = list(NULL)
  tree_tibble$reconstructed_dna_sequence = NA

  random_tip = sample(
    filter(
      tree_tibble,
      !node %in% parent,
      stringr::str_count(dna_sequence, "N") ==
        min(stringr::str_count(dna_sequence, "N"), na.rm = T)
    )$node,
    1
  )

  nts_to_root = tree_tibble$nt_mutations[getParents(tree_tibble, random_tip)]
  nts_to_root = rev(unlist(nts_to_root))

  root_sequence = applySubstitutions(
    filter(tree_tibble, node == random_tip)$dna_sequence,
    nts_to_root,
    backwards = T
  )

  tree_tibble$reconstructed_dna_sequence[
    tree_tibble$node == tree_tibble$parent
  ] = root_sequence

  nt_mutations = tree_tibble$nt_mutations
  nodes = tree_tibble$node
  parents = tree_tibble$parent
  reconstructed_dna_sequences = tree_tibble$reconstructed_dna_sequence

  for (tip in filter(tree_tibble, !node %in% parent)$node) {
    descs = getParentsFromNodelist(parents, tip)
    nts_to_tip = nt_mutations[descs]
    first_empty_slot = min(which(is.na(reconstructed_dna_sequences[descs])))
    sequence = reconstructed_dna_sequences[descs[first_empty_slot - 1]]

    for (p_i in first_empty_slot:length(descs)) {
      sequence = applySubstitutions(
        sequence,
        nts_to_tip[[p_i]]
      )
      reconstructed_dna_sequences[
        nodes == descs[[p_i]]
      ] = sequence
    }
  }

  tree_tibble$reconstructed_dna_sequence = reconstructed_dna_sequences

  tree_tibble
}


remakeTreeAndSequencesWithUsherASR = function(
  tree_and_sequences,
  aa_ref,
  nuc_ref
) {
  asr_and_tree = do_ASR_usher(
    tree_and_sequences$tree,
    setNames(
      tree_and_sequences$sequences$dna_sequence,
      tree_and_sequences$sequences$Isolate_unique_identifier
    ),
    nuc_ref
  )

  tree_and_sequences$tree = asr_and_tree$tree

  tree_and_sequences$tree_tibble = treeio::as_tibble(tree_and_sequences$tree)
  tree_and_sequences$tree_tibble = dplyr::left_join(
    tree_and_sequences$tree_tibble,
    tree_and_sequences$sequences,
    by = c("label" = "Isolate_unique_identifier")
  )

  tree_and_sequences$tree_tibble = cbind(
    tree_and_sequences$tree_tibble,
    asr_and_tree$asr[
      match(tree_and_sequences$tree_tibble$label, asr_and_tree$asr$node_id),
      -1
    ]
  )

  tree_and_sequences$tree_tibble$nt_mutations =
    tree_and_sequences$tree_tibble$nt_mutations %>%
    stringr::str_split(";|,") %>%
    purrr::map(~ .x[!is.na(.x)])

  tree_and_sequences$tree_tibble$nt_positions = purrr::map(
    tree_and_sequences$tree_tibble$nt_mutations,
    ~ as.integer(stringr::str_sub(.x, 2, -2))
  )

  tree_and_sequences$tree_tibble = reconstructNodeSequences(
    tree_and_sequences$tree_tibble
  )

  tree_and_sequences$tree_tibble = tree_and_sequences$tree_tibble %>%
    mutate(
      reconstructed_aa_sequence = seqUtils::translate(
        reconstructed_dna_sequence,
        aa_ref
      )
    )

  tree_and_sequences$tree_tibble$codon_changes = purrr::pmap(
    list(
      positions = tree_and_sequences$tree_tibble$nt_positions,
      seq = tree_and_sequences$tree_tibble$reconstructed_dna_sequence,
      parental_seq = tree_and_sequences$tree_tibble$reconstructed_dna_sequence[
        tree_and_sequences$tree_tibble$parent
      ]
    ),
    function(positions, seq, parental_seq) {
      if (length(positions) == 0) {
        return(character())
      }
      aa_positions = ceiling(positions / 3)
      paste0(
        stringr::str_sub(
          parental_seq,
          3 * (aa_positions - 1) + 1,
          3 * aa_positions
        ),
        ">",
        stringr::str_sub(seq, 3 * (aa_positions - 1) + 1, 3 * aa_positions)
      )
    }
  )
  tree_and_sequences$tree_tibble$aa_mutations = purrr::map2(
    tree_and_sequences$tree_tibble$codon_changes,
    tree_and_sequences$tree_tibble$nt_positions,
    function(codon, pos) {
      from = Biostrings::GENETIC_CODE[substr(codon, 1, 3)]
      from[is.na(from)] = "X"
      to = Biostrings::GENETIC_CODE[substr(codon, 5, 7)]
      to[is.na(to)] = "X"

      paste0(
        from,
        ceiling(pos / 3),
        to
      )
    }
  )

  parent_row = which(
    tree_and_sequences$tree_tibble$parent == tree_and_sequences$tree_tibble$node
  )

  tree_and_sequences$tree_tibble$nt_mutations[parent_row] = list(NULL)
  tree_and_sequences$tree_tibble$nt_positions[parent_row] = list(NULL)
  tree_and_sequences$tree_tibble$aa_mutations[parent_row] = list(NULL)
  tree_and_sequences$tree_tibble$codon_changes[parent_row] = list(NULL)

  tree_and_sequences$tree_tibble$nt_mutations_is_single = purrr::map(
    tree_and_sequences$tree_tibble$nt_positions,
    function(positions) {
      aa_positions = ceiling(positions / 3)
      !(duplicated(aa_positions, fromLast = F) |
        duplicated(aa_positions, fromLast = T))
    }
  )

  tree_and_sequences$tree_tibble$nt_mutations_is_syn = purrr::map(
    tree_and_sequences$tree_tibble$codon_changes,
    function(codon_changes) {
      froms = stringr::str_sub(codon_changes, 1, 3)
      tos = stringr::str_sub(codon_changes, -3, -1)

      froms = Biostrings::GENETIC_CODE[froms]
      froms[is.na(froms)] = "X"
      tos = Biostrings::GENETIC_CODE[tos]
      tos[is.na(tos)] = "X"

      froms == tos
    }
  )

  tree_and_sequences$tree_tibble = as_tibble(tree_and_sequences$tree_tibble)

  tree_and_sequences
}
