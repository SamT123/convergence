poisson_lognormal_nll <- function(params, y) {
  mu <- params[1]
  sigma <- exp(params[2]) # enforce sigma > 0

  # Marginal probability for each y_i
  logliks <- sapply(y, function(yi) {
    f <- function(lambda) {
      dpois(yi, lambda) * dlnorm(lambda, meanlog = mu, sdlog = sigma)
      # dpois(yi, lambda) * 2**dnorm(lambda, mean = mu, sd = sigma)
    }
    val <- integrate(
      f,
      lower = 0,
      upper = Inf,
      rel.tol = 1e-10,
      abs.tol = 0
    )$value
    if (val <= 0) {
      return(-Inf)
    }
    log(val)
  })

  return(-sum(logliks)) # negative log-likelihood
}

fit_poisson_lognormal <- function(
  y,
  init_mu = log(mean(y) + 0.1),
  init_sigma = sd(log(y + 0.1)),
  attempts = 1000
) {
  fits = list()
  for (i in seq_len(attempts)) {
    if (i > 1 & i %% 100 == 1) {
      message("attempt ", i, " / ", attempts)
    }

    mu_mask = c(1, 1, -1, -1)
    sd_mask = c(1, -1, 1, -1)

    # fmt: skip
    par_jitter = c(init_mu, log(init_sigma)) *
      2**(
          c(
            mu_mask[((i - 1) %% 4 + 1)],
            sd_mask[((i - 1) %% 4 + 1)]
          ) *
          (i - 1) / attempts
         )

    opt <- try(
      optim(
        par = par_jitter,
        fn = poisson_lognormal_nll,
        y = y,
        method = "BFGS",
        control = list(maxit = 100000)
      )
    )

    if (class(opt) != "try-error") {
      fits = append(
        fits,
        list(opt)
      )
    }
  }

  min_nll = which.min(purrr::map_dbl(fits, "value"))
  opt = fits[[min_nll]]

  list(
    mu_hat = opt$par[1],
    sigma_hat = exp(opt$par[2]),

    mu_hat_linear = exp(opt$par[1]),
    sigma_hat_log2 = log2(exp(exp(opt$par[2]))),

    nll = opt$value,
    convergence = opt$convergence
  )
}


extractParametersFromMaxLikeNormalModel = function(model_fit) {
  list(
    centrality = model_fit[["mu_hat_linear"]],
    scale = model_fit[["sigma_hat_log2"]]
  )
}

sampleFromMaxLikeNormalModel = function(n, parameters, pois = T, ...) {
  means = 2**(log2(parameters$centrality) + rnorm(n, 0, parameters$scale))

  if (pois) {
    return(rpois(n, means))
  } else {
    return(means)
  }
}

#' Synonymous rate variation models
#'@export
maxlike_models = list(
  normal_model = list(
    fit = fit_poisson_lognormal,
    extract_parameters_function = extractParametersFromMaxLikeNormalModel,
    sampling_function = sampleFromMaxLikeNormalModel
  )
)
