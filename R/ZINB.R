#' Fit Zero-Inflated Negative Binomial Model with Arbitrary Covariates and Prediction
#'
#'
#' Fits a zero-inflated negative binomial (ZINB) model using JAGS, with an optional
#' design matrix of covariates and full inprod for mean structure, and
#' can generate posterior predictive counts for new covariate data.
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
ZINB <- function(
    cases,
    pop              = NULL,
    covariates_count = NULL,
    covariates_zero  = NULL,
    covariatespred_count = NULL,
    covariatespred_zero  = NULL,
    poppred          = NULL,
    casespred        = NULL,
    beta_init        = NULL,
    delta_init       = NULL,
    r_init           = NULL,
    beta_prior_mean  = 0,
    beta_prior_sd    = 10,
    delta_prior_mean = 0,
    delta_prior_sd   = 10,
    r_prior_shape    = 1,
    r_prior_rate     = 1,
    n_iter           = 100000,
    n_burnin         = 10000,
    n_chains         = 3,
    n_thin           = 1,
    save_params      = c("beta","delta","r")
) {
  if (!requireNamespace("R2jags", quietly = TRUE))
    stop("Package R2jags is required.")
  
  N <- length(cases)
  
  # --- Build design matrices for estimation ---
  Xc <- if (!is.null(covariates_count)) {
    mat <- as.matrix(covariates_count)
    if (nrow(mat) != N) stop("covariates_count must match length of cases.")
    cbind(Intercept=1, mat)
  } else  matrix(1, nrow = N, ncol = 1)
  Kc <- ncol(Xc)
  Xz <- if (!is.null(covariates_zero)) {
    mat <- as.matrix(covariates_zero)
    if (nrow(mat) != N) stop("covariates_zero must match length of cases.")
    cbind(Intercept=1, mat)
  } else matrix(1, nrow = N, ncol = 1)
  Kz <- ncol(Xz)
  
  # --- Population offsets ---
  pop_vec <- if (is.null(pop)) rep(1, N) else as.numeric(pop)
  off_str <- if (is.null(pop)) "" else "log(pop[t]) + "
  
  # --- JAGS model string ---
  model_string <- paste(
    "model{",
    " for(t in 1:N){",
    "   ze[t] ~ dbern(pi[t])",
    "   y[t]  ~ dnegbin(pr[t], r)",
    "   pr[t] <- r / (r + (1-ze[t])*mu[t])",
    paste0("   mu[t]  <- exp(", off_str,
           "inprod(Xc[t,1:Kc], beta[1:Kc])) * (1 - ze[t])"),
    paste0("   logit(pi[t]) <- inprod(Xz[t,1:Kz], delta[1:Kz])"),
    " }",
    " for(k in 1:Kc){ beta[k]  ~ dnorm(", beta_prior_mean,
    ", 1/", beta_prior_sd^2, ") }",
    " for(k in 1:Kz){ delta[k] ~ dnorm(", delta_prior_mean,
    ", 1/", delta_prior_sd^2, ") }",
    paste0(" r ~ dgamma(", r_prior_shape, ", ", r_prior_rate, ")"),
    "}", sep = "\n"
  )
  
  # --- Write and run JAGS ---
  model_file <- tempfile(fileext = ".bug")
  writeLines(model_string, model_file)
  on.exit(unlink(model_file), add = TRUE)
  
  # --- Initial values ---
  if (is.null(beta_init))  beta_init  <- replicate(n_chains, rep(0, Kc), simplify=FALSE)
  if (is.null(delta_init)) delta_init <- replicate(n_chains, rep(0, Kz), simplify=FALSE)
  if (is.null(r_init))     r_init     <- rep(1, n_chains)
  inits <- lapply(seq_len(n_chains), function(i) list(
    beta  = beta_init[[i]],
    delta = delta_init[[i]],
    r     = r_init[i]
  ))
  
  data4j <- list(y = cases, N = N, Xc = Xc, Xz = Xz, pop = pop_vec, Kc = Kc, Kz = Kz)
  
  jags_out <- R2jags::jags(
    data               = data4j,
    inits              = inits,
    parameters.to.save = save_params,
    model.file         = model_file,
    n.iter             = n_iter,
    n.burnin           = n_burnin,
    n.chains           = n_chains,
    n.thin             = n_thin
  )
  
  summary_df <- as.data.frame(jags_out$BUGSoutput$summary)
  summary_df$dic <- jags_out$BUGSoutput$DIC
  result <- list(
    mcmc_summary = summary_df,
    mcmc_samples = coda::as.mcmc(jags_out),
    dic          = summary_df$dic[1]
  )
  
  # --- Prediction if requested ---
  if (!is.null(covariatespred_count)) {
    Xc1p <- as.matrix(covariatespred_count)
    M    <- nrow(Xc1p)
    if (ncol(Xc1p) != Kc-1) stop("covariatespred_count must match count covariates.")
    Xcp  <- cbind(Intercept=1, Xc1p)
    Xz1p <- if (!is.null(covariatespred_zero)) as.matrix(covariatespred_zero) else matrix(0, M, Kz-1)
    if (!is.null(covariatespred_zero) && ncol(Xz1p) != Kz-1) stop("covariatespred_zero must match zero covariates.")
    Xzp  <- cbind(Intercept=1, Xz1p)
    popp <- if (is.null(poppred)) rep(1, M) else as.numeric(poppred)
    
    sims   <- jags_out$BUGSoutput$sims.matrix
    beta_p <- sims[, grep("^beta", colnames(sims)), drop=FALSE]
    delta_p<- sims[, grep("^delta", colnames(sims)), drop=FALSE]
    r_p    <- sims[, "r"]
    npost  <- nrow(beta_p)
    
    pred_matrix <- matrix(NA, npost, M)
    for (i in seq_len(npost)) {
      for (t in seq_len(M)) {
        pi_t <- plogis(as.numeric(Xzp[t,] %*% delta_p[i,]))
        ze   <- rbinom(1, 1, pi_t)
        mu_t <- exp(log(popp[t]) + as.numeric(Xcp[t,] %*% beta_p[i,])) * (1-ze)
        pr   <- r_p[i] / (r_p[i] + (1-ze)*mu_t)
        pred_matrix[i,t] <- rnbinom(1, size=r_p[i], prob=pr)
      }
    }
    result$pred_matrix <- pred_matrix
    result$pred_mean   <- colMeans(pred_matrix)
    if (!is.null(casespred)) {
      if (length(casespred)!=M) stop("casespred must match number of prediction points.")
      result$mae  <- mean(abs(result$pred_mean - casespred))
      result$rmse <- sqrt(mean((result$pred_mean - casespred)^2))
    }
  }
  
  return(result)
}

