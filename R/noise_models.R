normal_model_specification = "data {
      array[N] int O;
      vector[N] expected_n;
      array[N] int mutation;
      array[N] int mutclass;
  }
  parameters{
       vector[N] bias;
       vector<lower=0>[12] std;
       //vector<lower=0>[12] nu;
       //vector<lower=0, upper = 1000>[12] linear_mean;
       vector[12] log_mean;
  }
  transformed parameters {
       //vector[12] log_mean = log2(linear_mean);

  }
  model{
      vector[N] lambda;

      //mutation_rate ~ student_t(nu[mutclass], log_mean[mutclass], std[mutclass]);

      bias ~ normal(log_mean[mutclass], std[mutclass]);

      for ( i in 1:N ) {
          lambda[i] = 2^(log2(expected_n[i]) + bias[mutation[i]]);
      }
      O ~ poisson( lambda );
  }
  generated quantities{
      vector[N] log_lik;
       vector[N] lambda;
      for ( i in 1:N ) {
          lambda[i] = 2^(log2(expected_n[i]) + bias[mutation[i]]);
      }
      for ( i in 1:N ) log_lik[i] = poisson_lpmf( O[i] | lambda[i] );
  }"

extractParametersFromNormalModel = function(model_fit) {
  mutation_order = model_fit$mutation_order
  parameters = list()
  for (i in seq_along(mutation_order)) {
    mutation = mutation_order[[i]]
    parameters[[mutation]] = list(
      centrality = bayestestR::map_estimate(
        model_fit$model_samples$log_mean[, i]
      )[[2]],
      scale = bayestestR::map_estimate(
        model_fit$model_samples$std[, i]
      )[[2]]
    )
  }
  parameters
}

sampleFromNormalModel = function(n, parameters, mutations_per_site, pois = T) {
  means = 2**(log2(mutations_per_site) +
    rnorm(
      n,
      parameters$centrality,
      parameters$scale
    ))

  if (pois) {
    return(rpois(n, means))
  } else {
    return(means)
  }
}


t_model_specification = "data {
      array[N] int O;
      vector[N] expected_n;
      array[N] int mutation;
      array[N] int mutclass;
  }
  parameters{
       vector[N] bias;
       vector<lower=0>[12] std;
       vector<lower=0>[12] nu;
       //vector<lower=0, upper = 1000>[12] linear_mean;
       vector[12] log_mean;
  }
  transformed parameters {
       //vector[12] log_mean = log2(linear_mean);

  }
  model{
      vector[N] lambda;
      nu ~ exponential( 0.5 );
      bias ~ student_t(nu[mutclass], log_mean[mutclass], std[mutclass]);

      for ( i in 1:N ) {
          lambda[i] = 2^(log2(expected_n[i]) + bias[mutation[i]]);
      }
      O ~ poisson( lambda );
  }
  generated quantities{
      vector[N] log_lik;
       vector[N] lambda;
      for ( i in 1:N ) {
          lambda[i] = 2^(log2(expected_n[i]) + bias[mutation[i]]);
      }
      for ( i in 1:N ) log_lik[i] = poisson_lpmf( O[i] | lambda[i] );
  }"

extractParametersFromTModel = function(model_fit) {
  mutation_order = model_fit$mutation_order
  parameters = list()
  for (i in seq_along(mutation_order)) {
    mutation = mutation_order[[i]]
    parameters[[mutation]] = list(
      centrality = bayestestR::map_estimate(
        model_fit$model_samples$log_mean[, i]
      )[[2]],
      scale = bayestestR::map_estimate(
        model_fit$model_samples$std[, i]
      )[[2]],
      df = bayestestR::map_estimate(
        model_fit$model_samples$nu[, i]
      )[[2]]
    )
  }
  parameters
}

sampleFromTModel = function(n, parameters, mutations_per_site, pois = T) {
  means = 2**(log2(mutations_per_site) +
    ggdist::rstudent_t(
      n,
      parameters$df,
      parameters$centrality,
      parameters$scale
    ))

  if (pois) {
    return(rpois(n, means))
  } else {
    return(means)
  }
}

#' Synonymous rate variation models
models = list(
  normal_model = list(
    specification = normal_model_specification,
    extract_parameters_function = extractParametersFromNormalModel,
    sampling_function = sampleFromNormalModel
  ),

  t_model = list(
    specification = t_model_specification,
    extract_parameters_function = extractParametersFromTModel,
    sampling_function = sampleFromTModel
  )
)


fitNoiseModel = function(
  synonymous_nuc_table,
  model_template,
  n_chains = 3,
  n_iter = 2000,
  n_noise_samples = 10000,
  seed = NULL,
  model_name = NULL
) {
  if (is.null(seed)) {
    seed = sample.int(.Machine$integer.max, 1)
  }

  mutation_order = unique(
    paste0(synonymous_nuc_table$from, synonymous_nuc_table$to)
  )

  synonymous_nuc_table_formatted = synonymous_nuc_table %>%
    mutate(
      mutclass = as.integer(factor(mutation_class, mutation_order)),
      mutation = as.integer(factor(mutation, unique(mutation))),
      expected_n = expected_n,
      O = n
    ) %>%
    select(mutation, mutclass, expected_n, O)

  model = stringr::str_replace_all(
    model_template,
    "N",
    as.character(nrow(synonymous_nuc_table_formatted))
  )

  model_fit <- rstan::stan(
    model_code = model,
    data = synonymous_nuc_table_formatted,
    verbose = FALSE,
    chains = n_chains,
    iter = n_iter,
    seed = seed
  )

  model_samples = rstan::extract(model_fit)

  list(
    model_fit = model_fit,
    model_samples = model_samples,
    mutation_order = mutation_order
  )
}
#
# getNoiseModelFromNucCounts = function(
#     nuc_counts,
#     model,
#     model_name,
#     n_chains = 3,
#     n_iter = 2000
# ) {
#   nuc_counts_long = nuc_counts %>%
#     mutate(
#       n = n_mutations_individual_positions,
#       mutation = purrr::pmap(
#         list(from, to, n),
#         \(f,t,x) paste0(f, names(x), t)
#       ),
#       mutation_class = paste0(from, to),
#       expected_n = mutations_per_site
#     ) %>%
#     tidyr::unnest(c(mutation, n)) %>%
#     select(
#       from,
#       to,
#       mutation_class,
#       mutation,
#       expected_n,
#       n
#     )
#
#   fit = fitNoiseModel(
#     nuc_counts_long,
#     model,
#     n_chains = n_chains,
#     n_iter = n_iter,
#     model_name = model_name
#   )
#
#   fit
# }

addPValuesToMutationTable = function(
  mutation_table,
  nuc_counts,
  sampling_function,
  alternative = c("greater", "two.sided", "less")
) {
  alternative = match.arg(alternative)

  whole_tree_parameters = setNames(
    nuc_counts$parameters,
    nuc_counts$from_to
  )

  mutation_table$n_samples = purrr::map(
    mutation_table$components,
    function(components) {
      component_samples = purrr::map(
        components,
        function(component) {
          params = whole_tree_parameters[[component[["nt_mutation_class"]]]]
          params[["centrality"]] = component[["expected_n"]]
          sampling_function(
            n = 1e4,
            parameters = params
          )
        }
      )

      purrr::reduce(
        component_samples,
        \(x, y) x + y
      )
    }
  )

  mutation_table$p = purrr::map2_dbl(
    mutation_table$n,
    mutation_table$n_samples,
    ~ hypothesisTest(
      samples = .y,
      observed = .x,
      alternative = alternative
    )
  )
  select(mutation_table, -n_samples)
}

hypothesisTest = function(
  samples,
  observed,
  alternative = c("greater", "two.sided", "less")
) {
  alternative = match.arg(alternative)
  p = dplyr::case_match(
    alternative,
    "greater" ~ mean(samples >= observed),
    "less" ~ mean(samples <= observed),
    "two.sided" ~ 2 * min(mean(samples <= observed), mean(samples >= observed))
  )

  min(p, 1)
}

#' Add CIs to mutation table
#'
#' inverts the hypothesis test in addPValuesToMutationTable to compute
#' confidence intervals (or compatibility set) for convergence ratios. CIs are
#' constructed by inversion to CI overlap with ratio=1.0 corresponds to p < alpha.
#'
#' exported because it is slow, so is not added to mutation_table by default
#'
#' @param mutation_table a mutation table
#' @param nuc_counts from getTreeInfo
#' @param sampling_function from getTreeInfo
#' @param confidence_level default 0.95
#' @param alternative "two.sided", "greater" (matching addPValuesToMutationTable), or "less"
#'
#' @return mutation_table with ratio_lower and ratio_upper cols
#'
#' @export
addCIsToMutationTable = function(
  mutation_table,
  nuc_counts,
  sampling_function,
  confidence_level = 0.95,
  alternative = c("two.sided", "greater", "less")
) {
  alternative = match.arg(alternative)

  whole_tree_parameters = setNames(
    nuc_counts$parameters,
    nuc_counts$from_to
  )

  alpha = 1 - confidence_level

  # helper function
  computePValueForRatio = function(
    ratio,
    observed,
    components,
    n_samples = 1e4
  ) {
    component_samples = purrr::map(
      components,
      function(component) {
        params = whole_tree_parameters[[component[["nt_mutation_class"]]]]
        params[["centrality"]] = ratio * component[["expected_n"]] # scale centrality by ratio
        sampling_function(
          n = n_samples,
          parameters = params
        )
      }
    )

    null_samples = purrr::reduce(
      component_samples,
      \(x, y) x + y
    )

    hypothesisTest(
      samples = null_samples,
      observed = observed,
      alternative = alternative
    )
  }

  findCIBounds = function(observed, expected_n, components) {
    observed_ratio = observed / expected_n

    if (alternative == "two.sided") {
      # 2 sided CI: p-value has n-shape as function of ratio
      ## (low at extremes, high near observed_ratio)
      # need to find where p = alpha on each sides of observed_ratio

      if (observed == 0) {
        lower_bound = 0
      } else {
        # lower bound is below observed ratio
        tryCatch(
          {
            lower_result = uniroot(
              function(ratio) {
                computePValueForRatio(ratio, observed, components) - alpha
              },
              interval = c(2**-12, observed_ratio),
              tol = 0.01
            )
            lower_bound = lower_result$root
          },
          error = function(e) {
            lower_bound = 0
          }
        )
      }

      # upper bound is above observed ratio
      tryCatch(
        {
          upper_result = uniroot(
            function(ratio) {
              computePValueForRatio(ratio, observed, components) - alpha
            },
            interval = c(observed_ratio, 2**12),
            tol = 0.01
          )
          upper_bound = upper_result$root
        },
        error = function(e) {
          upper_bound = Inf
        }
      )
    } else {
      # 1 sided CI (where is p = alpha?)

      if (observed == 0 && alternative == "greater") {
        bound = 0
      } else {
        tryCatch(
          {
            result = uniroot(
              function(ratio) {
                computePValueForRatio(ratio, observed, components) - alpha
              },
              interval = c(2**-12, 2**12),
              tol = 0.01
            )
            bound = result$root
          },
          error = function(e) {
            bound = if (alternative == "greater") 0 else Inf
          }
        )
      }

      if (alternative == "greater") {
        lower_bound = bound
        upper_bound = Inf
      } else {
        lower_bound = 0
        upper_bound = bound
      }
    }

    list(lower = lower_bound, upper = upper_bound)
  }

  cis = purrr::pmap(
    list(
      mutation_table$n,
      mutation_table$expected_n,
      mutation_table$components
    ),
    findCIBounds
  )

  mutation_table$ratio_lower = purrr::map_dbl(cis, "lower")
  mutation_table$ratio_upper = purrr::map_dbl(cis, "upper")

  mutation_table
}
