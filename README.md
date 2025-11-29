# convergence

<!-- badges: start -->

[![R-CMD-check](https://github.com/SamT123/convergence/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SamT123/convergence/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

R package for detecting convergent substitutions in phylogenetic trees.

## Installation

```r
# install.packages("devtools")
devtools::install_github("SamT123/convergence")
```

## Usage

```r
library(convergence)

# Create tree and sequences object
tree_and_sequences <- makeTreeAndSequences(tree, sequences)

# Add ancestral state reconstruction
tree_and_sequences <- addASRusher(tree_and_sequences, aa_ref, nuc_ref)

# Calculate tree size and nucleotide rates
tree_info <- getTreeSizeAndNucRates(
  tree_and_sequences,
  noise_model = maxlike_models$normal_model
)

# Get substitution ratios for positions of interest
ratios <- getSubstitutionRatios(
  tree_and_sequences,
  tree_info$nuc_rates,
  tree_info$tree_size_fn,
  tree_info$sampling_function,
  positions = 1:550
)
```

## External Dependencies

Requires external tools:

- [UShER](https://github.com/yatisht/usher) (matUtils, faToVcf)
- [Chronumental](https://github.com/theosanderson/chronumental) (optional, for node dating)
