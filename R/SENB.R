#' Fit Self-Exciting Negative Binomial Model with Prediction
#'
#' Fits a self-exciting negative binomial (SE-NB) model using JAGS, with an optional
#' design matrix of covariates and full inprod for mean structure, and
#' can generate posterior predictive counts including self-excitation.
#'

#'
#' @importFrom R2jags jags
#' @importFrom coda as.mcmc
#' @param cases Vector of observed counts (length N)
#' @param pop   Optional vector of population offsets (length N)
#' @param covariates Optional numeric matrix (N x P) of covariates for the count component.
#' @param covariatespred Optional numeric matrix (M x P) of new covariates for count prediction.
#' @param poppred          Optional vector of population offsets (length M) for prediction.
#' @param casespred        Optional vector of true counts (length M) for prediction performance.
#' @param casesoldold        Optional parameter of the cases of 1 timepoint previous than the start of timepoints fit.
#' @param casesoldpred        Optional parameter of the cases of 1 timepoint previous than the start of the prediction.
#' @param beta_init        Optional list of length n_chains for beta, count coefficients initial values.
#' @param r_init           Optional numeric vector of length n_chains for dispersion parameter.
#' @param beta_prior_mean  Mean for beta prior (default: 0)
#' @param beta_prior_sd    SD   for beta prior (default: 10)
#' @param r_prior_shape    Shape for r ~ dgamma (default: 1)
#' @param r_prior_rate     Rate  for r ~ dgamma (default: 1)
#' @param n_iter           Total MCMC iterations (default: 100000)
#' @param n_burnin         Burn-in iterations (default: 10000)
#' @param n_chains         Number of chains (default: 3)
#' @param n_thin           Thinning interval (default: 1)
#' @param save_params      Character vector of parameters to save (default c("beta","delta","r"))
#' @return A list with MCMC summary, samples, DIC, and if prediction data provided:
#'         prediction_matrix, prediction_mean, mae, rmse
#'
SENB <- function(
    cases,
    pop         = NULL,
    casesoldold=0,
    covariates  = NULL,
    covariatespred    = NULL,
    poppred     = NULL,
    casesoldpred = 0,
    casespred   = NULL,
    beta_init   = NULL,
    r_init      = NULL,
    beta_prior_mean = 0,
    beta_prior_sd   = 10,
    r_prior_shape   = 1,
    r_prior_rate    = 1,
    n_iter     = 100000,
    n_burnin   = 10000,
    n_chains   = 3,
    n_thin     = 1,
    save_params = c("beta","r","eta")
) {
  if (!requireNamespace("R2jags", quietly=TRUE)) stop("Package R2jags must be installed.")

  N <- length(cases)
  # Build design matrix
  if (!is.null(covariates)) {
    X1 <- as.matrix(covariates)
    if (nrow(X1)!=N) stop("covariates must match length of cases.")
  } else {
    X1 <- matrix(0, N, 0)
  }
  X <- cbind(Intercept=1, X1)
  K <- ncol(X)

  # Offsets
  if (is.null(pop)) {
    pop_vec <- rep(1, N)
    offset  <- ""
    offset1 <- ""
  } else {
    pop_vec <- as.numeric(pop)
    offset  <- "log(pop[t]) + "
    offset1 <- "log(pop[1]) + "
  }

  # JAGS model string (uses casesoldold in model)
  model_string <- paste(
    "model{",
    "  Y[1] ~ dnegbin(pr[1], r)",
    "  pr[1] <- r/(r + mu[1])",
    "  mu[1] <- mu0[1] + eta * casesoldold",
    "  mu0[1] <- exp(lambda0[1])",
    paste0("  lambda0[1] <- ", offset1, "inprod(X[1,1:K], beta[1:K])"),
    "  for(t in 2:N){",
    "    Y[t] ~ dnegbin(pr[t], r)",
    "    pr[t] <- r/(r + mu[t])",
    "    mu[t] <- mu0[t] + eta * Y[t-1]",
    "    mu0[t] <- exp(lambda0[t])",
    paste0("    lambda0[t] <- ", offset, "inprod(X[t,1:K], beta[1:K])"),
    "  }",
    "  # Priors",
    paste0("  for(k in 1:K){ beta[k] ~ dnorm(", beta_prior_mean, ", 1/", beta_prior_sd^2, ") }"),
    paste0("  r ~ dgamma(", r_prior_shape, ", ", r_prior_rate, ")"),
    "  eta ~ dbeta(1,1)",
    "}", sep="\n"
  )

  # Write model and run
  model_file <- tempfile(fileext=".bug")
  writeLines(model_string, model_file)
  on.exit(unlink(model_file))

  # Initial values
  if (is.null(beta_init)) beta_init <- lapply(1:n_chains, function(i) rep(0, K))
  if (is.null(r_init))    r_init    <- seq(0.5, 0.5 + 0.5 * (n_chains - 1), length.out = n_chains)
  inits <- lapply(1:n_chains, function(i) list(
    beta = beta_init[[i]],
    r    = r_init[i],
    eta  = 0.5
  ))

  data_jags <- list(
    Y       = cases,
    N       = N,
    X       = X,
    pop     = pop_vec,
    K       = K,
    casesoldold = casesoldold
  )

  jags_out <- R2jags::jags(
    data              = data_jags,
    inits             = inits,
    parameters.to.save= save_params,
    model.file        = model_file,
    n.iter            = n_iter,
    n.burnin          = n_burnin,
    n.chains          = n_chains,
    n.thin            = n_thin
  )

  sum_df <- as.data.frame(jags_out$BUGSoutput$summary)
  sum_df$dic <- jags_out$BUGSoutput$DIC
  res <- list(
    mcmc_summary = sum_df,
    mcmc_samples = coda::as.mcmc(jags_out),
    dic          = sum_df$dic[1]
  )

  # WAIC
  waic_s <- rjags::jags.samples(jags_out$model, c("WAIC","deviance"), type="mean", n.iter=1000)
  p_waic <- sum(waic_s$WAIC)
  dev    <- sum(waic_s$deviance)
  res$waic <- round(c(waic = dev + p_waic, p_waic = p_waic), 1)

  # Prediction block uses casesoldpred
  if (!is.null(covariatespred)) {
    Xp1 <- as.matrix(covariatespred)
    M   <- nrow(Xp1)
    if (ncol(Xp1) != ncol(X1)) stop("covariatespred must match covariates cols.")
    Xpred <- cbind(Intercept=1, Xp1)
    sims <- jags_out$BUGSoutput$sims.matrix
    beta_post <- sims[, grep("^beta\\[", colnames(sims)), drop=FALSE]
    r_post    <- sims[, "r"]
    eta_post  <- sims[, "eta"]
    npost <- nrow(beta_post)
    pred_mat <- matrix(NA, npost, M)
    for (i in 1:npost) {
      # first step uses casesoldpred
      mu0_1       <- exp((if (is.null(poppred)) 0 else log(poppred[1])) +
                           as.numeric(Xpred[1, ] %*% beta_post[i, ]))
      mu_1        <- mu0_1 + eta_post[i] * casesoldpred
      pr_1        <- r_post[i] / (r_post[i] + mu_1)
      pred_mat[i,1] <- rnbinom(1, size = r_post[i], prob = pr_1)
      # subsequent steps
      for (t in 2:M) {
        mu0_t <- exp((if (is.null(poppred)) 0 else log(poppred[t])) +
                       as.numeric(Xpred[t, ] %*% beta_post[i, ]))
        mu_t  <- mu0_t + eta_post[i] * pred_mat[i, t-1]
        pr_t  <- r_post[i] / (r_post[i] + mu_t)
        pred_mat[i,t] <- rnbinom(1, size = r_post[i], prob = pr_t)
      }
    }
    res$pred_matrix <- pred_mat
    res$pred_mean   <- colMeans(pred_mat)
    if (!is.null(casespred)) {
      if (length(casespred) != M) stop("casespred length must equal M")
      res$mae  <- mean(abs(res$pred_mean - casespred))
      res$rmse <- sqrt(mean((res$pred_mean - casespred)^2))
    }
  }

  return(res)
}
