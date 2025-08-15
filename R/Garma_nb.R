#' Fit Negative Binomial GARMA Model with Prediction
#'
#' This function fits a generalized autoregressive moving average (GARMA-NB)
#' model for count data using a negative binomial distribution, and optionally
#' generates posterior predictive counts for future covariate inputs.
#'
#'
#' @importFrom R2jags jags
#' @importFrom coda as.mcmc
#' @param cases Vector of observed counts (length N)
#' @param pop   Optional vector of population offsets (length N)
#' @param covariates Optional numeric matrix (N x P) of covariates for the count component.
#' @param covariatespred Optional numeric matrix (M x P) of new covariates for count prediction.
#' @param p Integer, autoregressive order
#' @param q Integer, moving average order
#' @param c Constant added before log (default 1)
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
#'         pred_matrix, pred_mean, mae, rmse
#'         
GARMA_NB <- function(
    cases,
    pop = NULL,
   covariates,
    p = 2,
    q = 2,
    c = 1,
    beta_init = NULL,
    r_init = NULL,
    beta_prior_mean = 0,
    beta_prior_sd = 10,
    r_prior_shape = 1,
    r_prior_rate = 1,
    n_iter = 100000,
    n_burnin = 10000,
    n_chains = 3,
    n_thin = 1,
    save_params = c("r", "beta", "phi", "theta"),
    covariatespred = NULL,
    poppred = NULL,
    casespred = NULL
) {
  if (!requireNamespace("R2jags", quietly = TRUE)) stop("Package R2jags must be installed.")
  
  N <- length(cases)
  X1 <- as.matrix(covariates)
  if (nrow(X1) != N) stop("covariates rows must equal length(cases)")
  X <- cbind(Intercept = 1, X1)
  K1 <- ncol(X)
  
  if (is.null(pop)) {
    pop_vec <- rep(1, N)
    offset_term <- ""
  } else {
    pop_vec <- as.numeric(pop)
    offset_term <- "log(pop[t]) +"
  }
  
  # Build single-assignment AR code
  if (p > 0) {
    ar_terms <- sapply(1:p, function(i) paste0(
      "phi[", i, "] * (log(max(c, y[t-", i, "])) - inprod(covariates[t-", i, ",], beta[]))"
    ))
    ar_code <- paste0("    ZAR[t] <- ", paste(ar_terms, collapse = " + "))
  } else {
    ar_code <- "    ZAR[t] <- 0"
  }
  
  # Build single-assignment MA code
  if (q > 0) {
    ma_terms <- sapply(1:q, function(i) paste0(
      "theta[", i, "] * (log(max(c, y[t-", i, "])) - mu[t-", i, "])"
    ))
    ma_code <- paste0("    ZMA[t] <- ", paste(ma_terms, collapse = " + "))
  } else {
    ma_code <- "    ZMA[t] <- 0"
  }
  
  # Compose model string
  model_string <- paste(
    "model {",
    "  for (t in 1:N) {",
    "    y[t] ~ dnegbin(pr[t], r)",
    "    pr[t] <- r / (r + lambda[t])",
    "    lambda[t] <- exp(mu[t])",
    paste0("    mu[t] <-", offset_term, "inprod(covariates[t,], beta[]) + ZAR[t] + ZMA[t]"),
    "  }",
    "  for (t in 1:p) { ZAR[t] <- 0 }",
    "  for (t in (p+1):N) {", ar_code, " }",
    "  for (t in 1:q) { ZMA[t] <- 0 }",
    "  for (t in (q+1):N) {", ma_code, " }",
    "  for (k in 1:p) { phi[k] ~ dunif(-1,1) }",
    "  for (k in 1:q) { theta[k] ~ dunif(-1,1) }",
    paste0("  for (k in 1:K1) { beta[k] ~ dnorm(", beta_prior_mean, ", 1/", beta_prior_sd^2, ") }"),
    paste0("  r ~ dgamma(", r_prior_shape, ", ", r_prior_rate, ")"),
    "}", sep = "\n"
  )
  model_file <- tempfile(fileext = ".bug")
  writeLines(model_string, model_file)
  on.exit(unlink(model_file))
  
  # Initial values
  if (is.null(beta_init)) beta_init <- lapply(1:n_chains, function(i) rep(0.3 * (i - 1), K1))
  if (is.null(r_init)) r_init <- seq(0.5, 0.5 + 0.5 * (n_chains - 1), length.out = n_chains)
  inits <- lapply(1:n_chains, function(i) list(beta = beta_init[[i]], r = r_init[i]))
  
  data4Jags <- list(
    y = cases,
    covariates = X,
    N = N,
    K1 = K1,
    pop = pop_vec,
    p = p,
    q = q,
    c = c
  )
  jags.out <- R2jags::jags(
    data = data4Jags,
    inits = inits,
    parameters.to.save = save_params,
    model.file = model_file,
    n.iter = n_iter,
    n.burnin = n_burnin,
    n.chains = n_chains,
    n.thin = n_thin
  )
  
  summary_df <- as.data.frame(jags.out$BUGSoutput$summary)
  summary_df$dic <- jags.out$BUGSoutput$DIC
  res <- list(
    mcmc_summary = summary_df,
    mcmc_samples = coda::as.mcmc(jags.out),
    dic = summary_df$dic[1]
  )
  
  # Prediction if covariatespred given
  if (!is.null(covariatespred)) {
    M <- nrow(covariatespred)
    Xp <- cbind(Intercept = 1, as.matrix(covariatespred))
    pops <- if (is.null(poppred)) rep(1, M) else poppred
    sims <- as.matrix(coda::as.mcmc(jags.out))
    
    #
    # —————————————————————————————————————
    #  extract MCMC samples for phi, theta, delta
    # —————————————————————————————————————
    #
    get_post_mat <- function(sims, name, length_expected) {
      # look for name[1], name[2], … 
      idx <- grep(paste0("^", name, "\\["), colnames(sims))
      # if we expected exactly 1 and found none, look for a scalar "name"
      if (length_expected == 1 && length(idx) == 0 && name %in% colnames(sims)) {
        idx <- which(colnames(sims) == name)
      }
      # return as a matrix (npost × length_expected, or 0 columns if none found)
      if (length(idx) > 0) {
        return( as.matrix(sims[, idx, drop = FALSE]) )
      } else {
        return( matrix(nrow = nrow(sims), ncol = 0) )
      }
    }
    
    npost      <- nrow(sims)
    beta_post  <- get_post_mat(sims, "beta",  K1)
    phi_post   <- if (p > 0) get_post_mat(sims, "phi",   p) else matrix(nrow = npost, ncol = 0)
    theta_post <- if (q > 0) get_post_mat(sims, "theta", q) else matrix(nrow = npost, ncol = 0)
    r_post     <- sims[, "r"]
    npost <- nrow(beta_post)
    pred_mat <- matrix(NA, npost, M)
    
    # initialize full series
    y_full <- c(cases, rep(NA, M))
    lam_full <- numeric(N + M)
    mu_full <- numeric(N + M)
    for (i in 1:npost) {
      lam_full[1:N] <- exp(drop(X %*% beta_post[i, ]))
      mu_full[1:N] <- lam_full[1:N]
      for (h in 1:M) {
        t <- N + h
        lam0 <- exp(
          (if (is.null(poppred)) 0 else log(pops[h])) +
            as.numeric(Xp[h, ] %*% beta_post[i, ])
        )
        # AR term
        zar <- 0
        if (p > 0) {
          for (j in 1:p) {
            idx <- t - j
            Xrow <- if (idx <= N) X[idx, ] else Xp[idx - N, ]
            resid_ar <- log(pmax(c, y_full[idx])) - as.numeric(Xrow %*% beta_post[i, ])
            zar <- zar + phi_post[i, j] * resid_ar
          }
        }
        # MA term
        zma <- 0
        if (q > 0) {
          for (j in 1:q) {
            resid_ma <- log(pmax(c, y_full[t - j])) - mu_full[t - j]
            zma <- zma + theta_post[i, j] * resid_ma
          }
        }
        # update forecasts
        mu_f <- lam0 * exp(zar + zma)
        pr_f <- r_post[i] / (r_post[i] + mu_f)
        y_full[t] <- rnbinom(1, size = r_post[i], prob = pr_f)
        pred_mat[i, h] <- y_full[t]
        lam_full[t] <- lam0
        mu_full[t] <- mu_f
      }
    }
    res$pred_matrix <- pred_mat
    res$pred_mean <- colMeans(pred_mat)
    if (!is.null(casespred)) {
      if (length(casespred) != M) stop("casespred length must equal M")
      res$mae <- mean(abs(res$pred_mean - casespred))
      res$rmse <- sqrt(mean((res$pred_mean - casespred)^2))
    }
  }
  
  return(res)
}

