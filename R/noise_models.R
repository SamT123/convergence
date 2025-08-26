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
  sampling_function
) {
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
    ~ mean(.y >= .x)
  )
  select(mutation_table, -n_samples)
}
