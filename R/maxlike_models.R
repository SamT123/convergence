# ---- utilities ----
log_sum_exp <- function(a) {
  # Stable log(sum(exp(a))) for a vector 'a'
  m <- max(a)
  m + log(sum(exp(a - m)))
}

# Core log-likelihood for Poisson-lognormal via Gauss-Hermite (probabilist form)
# y      : integer vector of counts
# mu     : mean of log lambda
# sigma  : sd of log lambda ( > 0 )
# nodes, weights: from statmod::gauss.quad.prob(n, "normal")
pln_loglik_gh <- function(y, mu, sigma, nodes, weights) {
  # For each y_i, integral wrt Z ~ N(0,1): E[ Poi(y_i | exp(mu + sigma Z)) ]
  # Because weights integrate against the standard normal, we just need dpois.
  # Do everything in log-space with log-sum-exp.
  ll_i <- vapply(
    y,
    function(yi) {
      # log pmf at all nodes
      lpmf_nodes <- dpois(yi, lambda = exp(mu + sigma * nodes), log = TRUE)
      # log-sum over quadrature nodes: log( sum_k w_k * exp(lpmf_k) )
      log_sum_exp(log(weights) + lpmf_nodes)
    },
    numeric(1)
  )
  sum(ll_i)
}

# Wrapper that 'optim' can use on unconstrained scale for sigma via log-sigma
# par = c(mu, log_sigma)
neg_loglik_par <- function(par, y, nodes, weights) {
  mu <- par[1]
  sigma <- exp(par[2]) # enforce sigma > 0
  -pln_loglik_gh(y, mu, sigma, nodes, weights)
}

# ---- main MLE function ----
pln_mle <- function(
  y,
  n_quad = 100,
  start = c(mean(log(y + 0.2)), sd(log(y + 0.2))),
  sigma_lower = 0.01,
  sigma_upper = 5
) {
  if (any(y < 0) || any(y != floor(y))) {
    stop("y must be nonnegative integers")
  }
  if (length(y) == 0) {
    stop("y is empty")
  }
  gh <- statmod::gauss.quad.prob(n_quad, "normal") # nodes/weights for standard normal
  nodes <- gh$nodes
  weights <- gh$weights

  start <- c(start[1], log(max(start[2], sigma_lower)))

  opt <- optim(
    par = start,
    fn = neg_loglik_par,
    y = y,
    nodes = nodes,
    weights = weights,
    method = "L-BFGS-B",
    lower = c(-Inf, log(sigma_lower)),
    upper = c(Inf, log(sigma_upper)),
    hessian = TRUE
  )

  # Extract estimates
  mu_hat <- opt$par[1]
  sigma_hat <- exp(opt$par[2])
  logLik_hat <- -opt$value

  # Standard errors via Hessian (on (mu, log sigma) scale), then delta method for sigma
  vcov_par <- tryCatch(solve(opt$hessian), error = function(e) {
    matrix(NA_real_, 2, 2)
  })
  se_mu <- sqrt(vcov_par[1, 1])
  se_logsigma <- sqrt(vcov_par[2, 2])
  se_sigma <- if (is.finite(se_logsigma)) sigma_hat * se_logsigma else NA_real_

  out <- list(
    par = c(mu = mu_hat, sigma = sigma_hat),
    par_scale = c(mu = exp(mu_hat), sigma = log2(exp(sigma_hat))),
    se = c(mu = se_mu, sigma = se_sigma),
    vcov_par = vcov_par, # vcov for (mu, log sigma)
    logLik = logLik_hat,
    convergence = opt$convergence,
    message = opt$message,
    n_quad = n_quad,
    call = match.call()
  )
  class(out) <- "pln_mle"
  out
}


fit_poisson_lognormal <- function(
  y,
  init_mu = log(mean(y) + 0.2),
  init_sigma = sd(log(y + 0.2))
) {
  opt = pln_mle(y)

  list(
    mu_hat = opt$par[1],
    sigma_hat = opt$par[2],

    mu_hat_linear = opt$par_scale[1],
    sigma_hat_log2 = opt$par_scale[2],

    nll = opt$logLik
  )
}

extractParametersFromMaxLikeNormalModel = function(model_fit) {
  list(
    centrality = model_fit[["mu_hat_linear"]],
    scale = model_fit[["sigma_hat_log2"]]
  )
}

sampleFromMaxLikeNormalModel = function(
  n,
  parameters,
  pois = T
) {
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
