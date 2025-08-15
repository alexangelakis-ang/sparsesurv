##' Fit Zero-Inflated Negative Binomial GARMA Model with Prediction
##'
##' This function fits a generalized autoregressive moving average (GARMA-ZINB)
##' model for count data using a zero-inflated negative binomial distribution,
##' allowing separate covariates for the count and zero-inflation parts,
##' and optionally generates posterior predictive counts for future covariate inputs.
##'
#'
#' @importFrom R2jags jags
#' @importFrom coda as.mcmc
#' @param cases Vector of observed counts (length N)
#' @param pop   Optional vector of population offsets (length N)
#' @param covariates_count Optional numeric matrix (N x P) of covariates for the count component.
#' @param covariates_zero  Optional numeric matrix (N x Q) of covariates for the zero-inflation component.
#' @param covariatespred_count Optional numeric matrix (M x P) of new covariates for count prediction.
#' @param covariatespred_zero  Optional numeric matrix (M x Q) of new covariates for zero-inflation prediction.
#' @param poppred          Optional vector of population offsets (length M) for prediction.
#' @param casespred        Optional vector of true counts (length M) for prediction performance.
#' @param p Integer, autoregressive order
#' @param q Integer, moving average order
#' @param c Constant added before log (default 1)
#' @param beta_init        Optional list of length n_chains for beta, count coefficients initial values.
#' @param delta_init       Optional list of length n_chains for delta, zero-inflation coefficients.
#' @param r_init           Optional numeric vector of length n_chains for dispersion parameter.
#' @param beta_prior_mean  Mean for beta prior (default: 0)
#' @param beta_prior_sd    SD   for beta prior (default: 10)
#' @param delta_prior_mean Mean for delta prior (default: 0)
#' @param delta_prior_sd   SD   for delta prior (default: 10)
#' @param r_prior_shape    Shape for r ~ dgamma (default: 1)
#' @param r_prior_rate     Rate  for r ~ dgamma (default: 1)
#' @param n_iter           Total MCMC iterations (default: 100000)
#' @param n_burnin         Burn-in iterations (default: 10000)
#' @param n_chains         Number of chains (default: 3)
#' @param n_thin           Thinning interval (default: 1)
#' @param save_params      Character vector of parameters to save (default c("beta","delta","r"))
#' @return A list with MCMC summary, samples, DIC, and if prediction data provided:
#'         pred_matrix, pred_mean, mae, rmse
#' @export
#'
GARMA_ZINB <- function(
    cases,
    pop = NULL,
    covariates_count = NULL,
    covariates_zero = NULL,
    p = 2,
    q = 2,
    c = 1,
    beta_init = NULL,
    delta_init = NULL,
    r_init = NULL,
    beta_prior_mean = 0,
    beta_prior_sd = 10,
    delta_prior_mean = 0,
    delta_prior_sd = 10,
    r_prior_shape = 1,
    r_prior_rate = 1,
    n_iter = 100000,
    n_burnin = 10000,
    n_chains = 3,
    n_thin = 1,
    save_params = c("r", "beta", "phi", "theta", "delta"),
    covariatespred_count = NULL,
    covariatespred_zero = NULL,
    poppred = NULL,
    casespred = NULL
) {
  if (!requireNamespace("R2jags", quietly = TRUE)) stop("Package R2jags must be installed.")

  N <- length(cases)
  # count part design
  if (!is.null(covariates_count)) {
    Xc1 <- as.matrix(covariates_count)
    if (nrow(Xc1) != N) stop("covariates_count must have same rows as cases")
  } else {
    Xc1 <- matrix(0, N, 0)
  }
  Xc <- cbind(Intercept = 1, Xc1)
  Kc <- ncol(Xc)

  # zero inflation design
  if (!is.null(covariates_zero)) {
    Xz1 <- as.matrix(covariates_zero)
    if (nrow(Xz1) != N) stop("covariates_zero must have same rows as cases")
  } else {
    Xz1 <- matrix(0, N, 0)
  }
  Xz <- cbind(Intercept = 1, Xz1)
  Kz <- ncol(Xz)

  # offsets
  if (is.null(pop)) {
    pop_vec <- rep(1, N)
    offset_str <- ""
  } else {
    pop_vec <- as.numeric(pop)
    offset_str <- "log(pop[t]) + "
  }

  # build AR/MA code
  ar_code <- if (p > 0) {
    paste0("    ZAR[t] <- ", paste(
      sapply(1:p, function(i) sprintf(
        "phi[%d] * (log(max(c, y[t-%d])) - mu[t-%d])", i, i, i
      )), collapse = " + "))
  } else {
    "    ZAR[t] <- 0"
  }
  ma_code <- if (q > 0) {
    paste0("    ZMA[t] <- ", paste(
      sapply(1:q, function(i) sprintf(
        "theta[%d] * (log(max(c, y[t-%d])) - mu[t-%d])", i, i, i
      )), collapse = " + "))
  } else {
    "    ZMA[t] <- 0"
  }

  # compose BUGS model
  model_string <- paste(
    "model {",
    "  for (t in 1:N) {",
    "    y[t] ~ dnegbin(pr[t], r)",
    "    pr[t] <- r / (r + (1+ze[t])*lambda[t])",
    "    lambda[t] <- exp(mu[t])",
    sprintf("    mu[t] <- %sinprod(Xc[t,], beta[]) + ZAR[t] + ZMA[t]", offset_str),
    "    ze[t] ~ dbern(pi[t])",
    sprintf("    logit(pi[t]) <- %sinprod(Xz[t,], delta[])", offset_str),
    "  }",
    "  for (t in 1:p) { ZAR[t] <- 0 }",
    sprintf("  for (t in %d:N) {", p+1),
    ar_code,
    "  }",
    "  for (t in 1:q) { ZMA[t] <- 0 }",
    sprintf("  for (t in %d:N) {", q+1),
    ma_code,
    "  }",
    "  for (k in 1:p) { phi[k] ~ dunif(-1,1) }",
    "  for (k in 1:q) { theta[k] ~ dunif(-1,1) }",
    sprintf("  for (k in 1:Kc) { beta[k] ~ dnorm(%f, 1/%f) }", beta_prior_mean, beta_prior_sd^2),
    sprintf("  r ~ dgamma(%f, %f)", r_prior_shape, r_prior_rate),
    sprintf("  for (k in 1:Kz) { delta[k] ~ dnorm(%f, 1/%f) }", delta_prior_mean, delta_prior_sd^2),
    "}", sep="\n"
  )
  model_file <- tempfile(fileext = ".bug")
  writeLines(model_string, model_file)
  on.exit(unlink(model_file), add = TRUE)

  # initials
  if (is.null(beta_init))  beta_init  <- replicate(n_chains, rep(0, Kc), simplify = FALSE)
  if (is.null(delta_init)) delta_init <- replicate(n_chains, rep(0, Kz), simplify = FALSE)
  if (is.null(r_init))     r_init     <- rep(1, n_chains)
  inits <- lapply(seq_len(n_chains), function(i) list(
    beta = beta_init[[i]], delta = delta_init[[i]], r = r_init[i],
    phi = rep(0, p), theta = rep(0, q)
  ))

  data_jags <- list(N = N, y = cases, Xc = Xc, Xz = Xz,
                    p = p, q = q, c = c, pop = pop_vec,
                    Kc = Kc, Kz = Kz)
  jags_out <- R2jags::jags(
    data = data_jags, inits = inits, parameters.to.save = save_params,
    model.file = model_file, n.iter = n_iter,
    n.burnin = n_burnin, n.chains = n_chains, n.thin = n_thin
  )

  summary_df <- as.data.frame(jags_out$BUGSoutput$summary)
  summary_df$dic <- jags_out$BUGSoutput$DIC
  res <- list(mcmc_summary = summary_df,
              mcmc_samples = coda::as.mcmc(jags_out),
              dic = summary_df$dic[1])

  ## Prediction
  if (!is.null(covariatespred_count)) {
    M <- nrow(covariatespred_count)

    # --- count prediction design ---
    if (!is.null(covariatespred_count)) {
      Cc1 <- as.matrix(covariatespred_count)
      if (nrow(Cc1) != M)
        stop("covariatespred_count must have same rows as covariatespred_count")
      if (ncol(Cc1) != Kc - 1)
        stop("covariatespred_count must have same columns as covariates_count")
    } else {
      # no new count covariates ⇒ all zero columns
      Cc1 <- matrix(0, nrow = M, ncol = Kc - 1)
    }
    Xc_pred <- cbind(Intercept = 1, Cc1)

    # --- zero‐inflation prediction design ---
    if (!is.null(covariatespred_zero)) {
      Zp1 <- as.matrix(covariatespred_zero)
      if (nrow(Zp1) != M)
        stop("covariatespred_zero must have same rows as covariatespred_count")
      if (ncol(Zp1) != Kz - 1)
        stop("covariatespred_zero must have same columns as covariates_zero")
    } else {
      # no new covariates ⇒ all zero columns
      Zp1 <- matrix(0, nrow = M, ncol = Kz - 1)
    }
    Xz_pred <- cbind(Intercept = 1, Zp1)

    pops <- if (is.null(poppred)) rep(1, M) else poppred
    sims <- as.matrix(coda::as.mcmc(jags_out))


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
    beta_post  <- get_post_mat(sims, "beta",  Kc)
    phi_post   <- if (p > 0) get_post_mat(sims, "phi",   p) else matrix(nrow = npost, ncol = 0)
    theta_post <- if (q > 0) get_post_mat(sims, "theta", q) else matrix(nrow = npost, ncol = 0)
    delta_post <- get_post_mat(sims, "delta", Kz)
    r_post     <- sims[, "r"]


    npost <- nrow(sims)
    pred_mat <- matrix(NA, npost, M)

    # full history container
    y_full    <- c(cases, rep(NA, M))
    mu_full   <- c(rep(NA, N+M))
    lam_full  <- c(rep(NA, N+M))
    # initialize mu and lam for observed
    for (i in 1:npost) {
      lam_full[1:N] <- exp(drop(Xc %*% beta_post[i, ]))
      mu_full[1:N]  <- lam_full[1:N]
      for (h in 1:M) {
        t_idx <- N + h
        # baseline
        lam0 <- exp((if (is.null(poppred)) 0 else log(pops[h])) +
                      drop(Xc_pred[h, ] %*% beta_post[i, ]))
        # AR
        zar <- 0
        if (p>0) for (j in 1:p) {
          idx <- t_idx - j
          mu_lag <- mu_full[idx]
          resid_ar <- log(max(c, y_full[idx])) - mu_lag
          zar <- zar + phi_post[i, j] * resid_ar
        }
        # MA
        zma <- 0
        if (q>0) for (j in 1:q) {
          idx <- t_idx - j
          resid_ma <- log(max(c, y_full[idx])) - mu_full[idx]
          zma <- zma + theta_post[i, j] * resid_ma
        }
        mu_f <- lam0 * exp(zar + zma)
        # zero inflation probability
        logit_pi <- (if (is.null(poppred)) 0 else log(pops[h])) +
          drop(Xz_pred[h, ] %*% delta_post[i, ])

        logit_pi_f <- 1/(1+exp(-logit_pi))
        zet <- rbinom(1,1,logit_pi_f)
        pr_f <- r_post[i] / (r_post[i] + (1-zet)*mu_f)
        # simulate
        y_sim <-  rnbinom(1, size=r_post[i], prob=pr_f)
        # record
        y_full[t_idx] <- y_sim
        lam_full[t_idx] <- lam0
        mu_full[t_idx] <- mu_f
        pred_mat[i,h] <- y_sim
      }
    }
    res$pred_matrix <- pred_mat
    res$pred_mean   <- colMeans(pred_mat)
    if (!is.null(casespred)) {
      if (length(casespred)!=M) stop("casespred length must equal M")
      res$mae  <- mean(abs(res$pred_mean-casespred))
      res$rmse <- sqrt(mean((res$pred_mean-casespred)^2))
    }
  }

  return(res)
}
