#' Fit Negative Binomial Model with Arbitrary Covariates
#'
#' Fits a negative binomial (NB) model using JAGS, with an optional
#' design matrix of covariates and full inprod for mean structure, and
#' can generate posterior predictive counts for new covariate data.
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
NB <- function(
    cases,
    pop         = NULL,
    casespred=NULL,
    covariates  = NULL,
    covariatespred    = NULL,
    poppred     = NULL,
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
    save_params = c("beta","r")
) {
  if (!requireNamespace("R2jags", quietly=TRUE))
    stop("Package R2jags must be installed.")

  N <- length(cases)
  # Construct design matrix with intercept
  if (!is.null(covariates)) {
    X1 <- as.matrix(covariates)
    if (nrow(X1)!=N) stop("covariates must match length of cases.")
  } else {
    X1 <- matrix(0, N, 0)
  }
  X <- cbind(Intercept=1, X1)
  K <- ncol(X)

  # Offsets for estimation
  if (is.null(pop)) {
    pop_vec <- rep(1, N)
    offset  <- ""
  } else {
    pop_vec <- as.numeric(pop)
    offset  <- "log(pop[t]) + "
  }

  # JAGS model string
  model_string <- paste(
    "model{",
    "  for(t in 1:N){",
    "    Y[t] ~ dnegbin(pr[t], r)",
    "    pr[t] <- r / (r + mu[t])",
    paste0("    mu[t]  <- exp(", offset,
           "inprod(X[t,1:K], beta[1:K]))"),
    "  }",
    "  for(k in 1:K){ beta[k] ~ dnorm(", beta_prior_mean,
    ", 1/", beta_prior_sd^2, ") }",
    paste0("  r ~ dgamma(", r_prior_shape, ", ", r_prior_rate, ")"),
    "}", sep="\n"
  )

  # Write model file
  model_file <- tempfile(fileext=".bug")
  writeLines(model_string, model_file)
  on.exit(unlink(model_file))

  # Initial values
  if (is.null(beta_init)) beta_init <- lapply(1:n_chains, function(i) rep(0, K))
  if (is.null(r_init))    r_init    <- seq(0.5, 0.5+0.5*(n_chains-1), length.out=n_chains)

  inits <- lapply(1:n_chains, function(i) list(
    beta = beta_init[[i]],
    r    = r_init[i]
  ))

  data_jags <- list(
    Y   = cases,
    N   = N,
    X   = X,
    pop = pop_vec,
    K   = K
  )

  # Run JAGS
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

  # Summaries
  summary_df <- as.data.frame(jags_out$BUGSoutput$summary)
  summary_df$dic <- jags_out$BUGSoutput$DIC
  result <- list(
    mcmc_summary = summary_df,
    mcmc_samples = coda::as.mcmc(jags_out),
    dic          = summary_df$dic[1]
  )


  if (!is.null(covariatespred)) {
    Xp1 <- as.matrix(covariatespred)
    M   <- nrow(Xp1)
    if (ncol(Xp1)!=ncol(X1)) stop("covariatespred must have same columns as covariates.")
    Xpred <- cbind(Intercept=1, Xp1)
    sims <- jags_out$BUGSoutput$sims.matrix
    beta_post <- sims[, grep("^beta\\[", colnames(sims)), drop=FALSE]
    r_post    <- sims[,"r"]
    npost <- nrow(beta_post)
    pred_matrix <- matrix(NA, npost, M)
    for (i in 1:npost) {
      for (t in 1:M) {
        mu <- exp((if(is.null(poppred)) 0 else log(poppred[t])) +
                    as.numeric(Xpred[t,] %*% beta_post[i,]))
        pr <- r_post[i] / (r_post[i] + mu)
        pred_matrix[i,t] <- rnbinom(1, size=r_post[i], prob=pr)
      }
    }
    result$pred_matrix <- pred_matrix
    result$pred_mean   <- colMeans(pred_matrix)
    # Compute MAE/RMSE if true values provided
    if (!is.null(casespred)) {
      if (length(casespred)!=M) stop("ytrue must match number of prediction points.")
      result$mae  <- mean(abs(result$pred_mean - casespred))
      result$rmse <- sqrt(mean((result$pred_mean - casespred)^2))
    }
  }

  return(result)
}
