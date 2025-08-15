#' Fit Self exciting Zero-Inflated Negative Binomial Model with Arbitrary Covariates and Prediction
#'
#'
#' Fits a Self exciting zero-inflated negative binomial (SE-ZINB) model using JAGS, with an optional
#' design matrix of covariates and full inprod for mean structure, and
#' can generate posterior predictive counts for new covariate data.
#'
#' @importFrom R2jags jags
#' @importFrom coda as.mcmc
#' @param cases Vector of observed counts (length N)
#' @param pop   Optional vector of population offsets (length N)
#' @param casesoldold        Optional parameter of the cases of 1 timepoint previous than the start of timepoints fit.
#' @param casesoldpred        Optional parameter of the cases of 1 timepoint previous than the start of the prediction.
#' @param covariates_count Optional numeric matrix (N x P) of covariates for the count component.
#' @param covariates_zero  Optional numeric matrix (N x Q) of covariates for the zero-inflation component.
#' @param covariatespred_count Optional numeric matrix (M x P) of new covariates for count prediction.
#' @param covariatespred_zero  Optional numeric matrix (M x Q) of new covariates for zero-inflation prediction.
#' @param poppred          Optional vector of population offsets (length M) for prediction.
#' @param casespred        Optional vector of true counts (length M) for prediction performance.
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


SEZINB <- function(
    cases,
    pop = NULL,
    casesoldold =0,
    covariates_count = NULL,
    covariates_zero  = NULL,
    covariatespred_count = NULL,
    covariatespred_zero  = NULL,
    poppred             = NULL,
    casesoldpred            = 0,
    casespred           = NULL,
    beta_init           = NULL,
    delta_init          = NULL,
    r_init              = NULL,
    beta_prior_mean     = 0,
    beta_prior_sd       = 10,
    delta_prior_mean    = 0,
    delta_prior_sd      = 10,
    r_prior_shape       = 1,
    r_prior_rate        = 1,
    n_iter              = 100000,
    n_burnin            = 10000,
    n_chains            = 3,
    n_thin              = 1,
    save_params         = c("beta", "delta", "r", "eta")
) {
  if (!requireNamespace("R2jags", quietly = TRUE)) stop("Package R2jags is required.")

  N <- length(cases)

  # Count covariate matrix with intercept
  if (!is.null(covariates_count)) {
    Xc1 <- as.matrix(covariates_count)
    if (nrow(Xc1) != N) stop("covariates_count must match length of cases.")
  } else {
    Xc1 <- matrix(0, nrow = N, ncol = 0)
  }
  Xc <- cbind(Intercept = 1, Xc1)
  Kc <- ncol(Xc)

  # Zero-inflation covariate matrix with intercept
  if (!is.null(covariates_zero)) {
    Xz1 <- as.matrix(covariates_zero)
    if (nrow(Xz1) != N) stop("covariates_zero must match length of cases.")
  } else {
    Xz1 <- matrix(0, nrow = N, ncol = 0)
  }
  Xz <- cbind(Intercept = 1, Xz1)
  Kz <- ncol(Xz)

  # Offsets
  if (is.null(pop)) {
    pop_vec <- rep(1, N)
    off_str <- ""
    off1    <- ""
  } else {
    pop_vec <- as.numeric(pop)
    off_str <- "log(pop[t]) + "
    off1    <- "log(pop[1]) + "
  }

  # Compose BUGS model string
  model_string <- paste(
    "model{",
    "  Y[1] ~ dnegbin(pr[1], r)",
    "  pr[1] <- r / (r + (1 - ze[1]) * mu[1]) - 1e-10 * ze[1]",
    "  mu[1] <- mu0[1] + eta * casesoldold",
    "  mu0[1] <- exp(lambda0[1])",
    paste0("  lambda0[1] <- ", off1, "inprod(Xc[1,1:Kc], beta[1:Kc])"),
    "  ze[1] ~ dbern(pi[1])",
    paste0("  pi[1] <- ilogit(", off1, "inprod(Xz[1,1:Kz], delta[1:Kz]))"),
    "  for(t in 2:N){",
    "    Y[t] ~ dnegbin(pr[t], r)",
    "    pr[t] <- r / (r + (1 - ze[t]) * mu[t]) - 1e-10 * ze[t]",
    "    mu[t] <- mu0[t] + eta * Y[t-1]",
    "    mu0[t] <- exp(lambda0[t])",
    paste0("    lambda0[t] <- ", off_str, "inprod(Xc[t,1:Kc], beta[1:Kc])"),
    "    ze[t] ~ dbern(pi[t])",
    paste0("    pi[t] <- ilogit(", off_str, "inprod(Xz[t,1:Kz], delta[1:Kz]))"),
    "  }",
    "  r   ~ dgamma(", r_prior_shape, ", ", r_prior_rate, ")",
    "  eta ~ dbeta(1,1)",
    paste0("  for(k in 1:Kc){ beta[k]  ~ dnorm(", beta_prior_mean, ", 1/", beta_prior_sd^2, ") }"),
    paste0("  for(k in 1:Kz){ delta[k] ~ dnorm(", delta_prior_mean, ", 1/", delta_prior_sd^2, ") }"),
    "}",
    sep = "\n"
  )

  model_file <- tempfile(fileext = ".bug")
  writeLines(model_string, model_file)
  on.exit(unlink(model_file))

  # Initial values
  if (is.null(beta_init))  beta_init  <- lapply(1:n_chains, function(i) rep(0, Kc))
  if (is.null(delta_init)) delta_init <- lapply(1:n_chains, function(i) rep(0, Kz))
  if (is.null(r_init))     r_init     <- seq(0.5, 0.5 + 0.5*(n_chains-1), length.out = n_chains)

  inits <- lapply(1:n_chains, function(i) list(
    beta  = beta_init[[i]],
    delta = delta_init[[i]],
    r     = r_init[i],
    eta   = 0.5
  ))

  data4Jags <- list(
    Y       = cases,
    N       = N,
    Xc      = Xc,
    Xz      = Xz,
    pop     = pop_vec,
    casesoldold = casesoldold,
    Kc      = Kc,
    Kz      = Kz
  )

  jags.out <- R2jags::jags(
    data               = data4Jags,
    inits              = inits,
    parameters.to.save = save_params,
    model.file         = model_file,
    n.iter             = n_iter,
    n.burnin           = n_burnin,
    n.chains           = n_chains,
    n.thin             = n_thin
  )

  # Summaries
  summary_df <- as.data.frame(jags.out$BUGSoutput$summary)
  summary_df$dic <- jags.out$BUGSoutput$DIC
  s <- rjags::jags.samples(jags.out$model,
                           c("WAIC","deviance"), type="mean", n.iter=1000)
  p_waic <- sum(s$WAIC); dev <- sum(s$deviance)
  waic_vals <- round(c(waic=dev+p_waic, p_waic=p_waic),1)

  # Base result
  res <- list(
    mcmc_summary = summary_df,
    mcmc_samples = coda::as.mcmc(jags.out),
    dic          = summary_df$dic[1],
    waic         = waic_vals
  )

  # Prediction block
  if (!is.null(covariatespred_count) && !is.null(covariatespred_zero)) {
    Xc_pred <- cbind(Intercept = 1, as.matrix(covariatespred_count))
    Xz_pred <- cbind(Intercept = 1, as.matrix(covariatespred_zero))
    M       <- nrow(Xc_pred)
    if (ncol(Xc_pred) != Kc) stop("covariatespred_count must have same columns as covariates_count + intercept.")
    if (ncol(Xz_pred) != Kz) stop("covariatespred_zero must have same columns as covariates_zero + intercept.")

    sims       <- jags.out$BUGSoutput$sims.matrix
    beta_post  <- sims[, grep("^beta\\[", colnames(sims)), drop = FALSE]
    delta_post <- sims[, grep("^delta\\[", colnames(sims)), drop = FALSE]
    r_post     <- sims[, "r"]
    eta_post   <- sims[, "eta"]
    npost      <- nrow(beta_post)
    pred_mat   <- matrix(NA, npost, M)

    for (i in 1:npost) {
      # t = 1
      mu0_1 <- exp((if (is.null(poppred)) 0 else log(poppred[1])) +
                     as.numeric(Xc_pred[1, ] %*% beta_post[i, ]))
      mu_1  <- mu0_1 + eta_post[i] * casesoldpred
      pi_1  <- plogis((if (is.null(poppred)) 0 else log(poppred[1])) +
                        as.numeric(Xz_pred[1, ] %*% delta_post[i, ]))
      ze_1  <- rbinom(1, 1, pi_1)
      if (ze_1 == 1) {
        pred_mat[i, 1] <- 0
      } else {
        pr_1 <- r_post[i] / (r_post[i] + mu_1)
        pred_mat[i, 1] <- rnbinom(1, size = r_post[i], prob = pr_1)
      }

      # t > 1
      for (t in 2:M) {
        mu0_t <- exp((if (is.null(poppred)) 0 else log(poppred[t])) +
                       as.numeric(Xc_pred[t, ] %*% beta_post[i, ]))
        mu_t  <- mu0_t + eta_post[i] * pred_mat[i, t - 1]
        pi_t  <- plogis((if (is.null(poppred)) 0 else log(poppred[t])) +
                          as.numeric(Xz_pred[t, ] %*% delta_post[i, ]))
        ze_t  <- rbinom(1, 1, pi_t)
        if (ze_t == 1) {
          pred_mat[i, t] <- 0
        } else {
          pr_t <- r_post[i] / (r_post[i] + mu_t)
          pred_mat[i, t] <- rnbinom(1, size = r_post[i], prob = pr_t)
        }
      }
    }

    res$pred_matrix <- pred_mat
    res$pred_mean   <- colMeans(pred_mat)
    if (!is.null(casespred)) {
      if (length(casespred) != M) stop("casespred must match number of prediction rows.")
      res$mae  <- mean(abs(res$pred_mean - casespred))
      res$rmse <- sqrt(mean((res$pred_mean - casespred)^2))
    }
  }

  return(res)
}
