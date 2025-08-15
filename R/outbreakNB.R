
#' Fit an outbreak detection negative binomial outbreak model
#'
#' @param cases            Integer or numeric vector of observed case counts (length N).
#' @param pop              (Optional) Numeric vector of population offsets (length N). If NULL, offset = 1.
#' @param covariates   (Optional) Data.frame or matrix of covariates for the count model (N x p_c).
#' @param beta_init        (Optional) List of length `n_chains` giving initial values for beta (each a vector of length p_c+1).
#' @param r_init           (Optional) Numeric vector of length `n_chains` giving initial values for the NB dispersion parameter.
#' @param beta_prior_mean  Prior mean for beta coefficients of the Negative binomial part (default = 0).
#' @param beta_prior_sd    Prior SD   for beta coefficients of the Negative binomial part (default = 10).
#' @param r_prior_shape    Shape parameter of a prior on r (default = 1).
#' @param r_prior_rate     Rate  parameter of b prior on r (default = 1).
#' @param p_priors         Alpha parameters for the binomial priors on p00 and p11 (default = 1).
#' @param n_iter           Total number of MCMC iterations per chain (default = 100000).
#' @param n_burnin         Number of burn-in iterations (default = 10000).
#' @param n_chains         Number of MCMC chains (default = 3).
#' @param n_thin           Thinning interval for MCMC samples (default = 1).
#' @param save_params      Character vector of parameter names to save (must include "Z").
#' @param dates            (Optional) Vector of Date or POSIX dates for plotting Z; if NULL, uses index 1:N.
#' @param plot_Z           Logical; if TRUE, returns a ggplot2 object of the posterior mean Z over time.

#' @return A list with MCMC summary, samples, DIC, WAIC, and plot of the probability of being in an epidemic state.



outbreakNB <- function(
    cases,
    pop = NULL,
    covariates = NULL,
    beta_init = NULL,
    r_init = NULL,
    beta_prior_mean = 0,
    beta_prior_sd = 10,
    r_prior_shape = 1,
    r_prior_rate = 1,
    p_priors = 1,
    n_iter = 100000,
    n_burnin = 10000,
    n_chains = 3,
    n_thin = 1,
    save_params = c("beta", "r", "Z"),
    dates = NULL,      # optional date vector aligned with cases
    plot_Z = FALSE     # whether to build the Z plot
) {
  if (!requireNamespace("R2jags", quietly = TRUE)) stop("Package R2jags is required.")
  if (!requireNamespace("coda", quietly = TRUE)) stop("Package coda is required.")
  if (!requireNamespace("rjags", quietly = TRUE)) stop("Package rjags is required for WAIC/dev calculations.")

  N <- length(cases)

  # if save_params was explicitly passed, ensure Z is included
  if (!missing(save_params)) {
    if (!("Z" %in% save_params)) {
      save_params <- c(save_params, "Z")
    }
  }

  # Count covariate matrix with intercept
  if (!is.null(covariates)) {
    Xc1 <- as.matrix(covariates)
    if (nrow(Xc1) != N) stop("covariates must have length equal to cases.")
  } else {
    Xc1 <- matrix(0, nrow = N, ncol = 0)
  }
  Xc <- cbind(Intercept = 1, Xc1)
  Kc <- ncol(Xc)



  # Offsets (your original names)
  if (is.null(pop)) {
    pop_vec <- rep(1, N)
    off_str <- ""
    off1    <- ""
  } else {
    pop_vec <- as.numeric(pop)
    off_str <- "log(pop[t]) + "
    off1    <- "log(pop[1]) + "
  }

  # Build model string using paste / paste0 with inline precision as you had it
  model_lines <- c(
    "model{",
    "  # First time point",
    "  Y[1] ~ dnegbin(pr[1], r)",
    "  pr[1] <- r / (r + mu[1])",
    "  mu[1] <- mu0[1] + Z[1] * mu1[1]",
    "  mu0[1] <- exp(lambda0[1])",
    "  mu1[1] <- 0",
    paste0("  lambda0[1] <- ", off1, "inprod(Xc[1,1:Kc], beta[1:Kc])"),
    "",
    "  # Subsequent time points",
    "  for(t in 2:N){",
    "    Y[t] ~ dnegbin(pr[t], r)",
    "   pr[t] <- r / (r + mu[t])",
    "    mu[t] <- mu0[t] + Z[t] * mu1[t]",
    "    mu0[t] <- exp(lambda0[t])",
    "    mu1[t] <- mu[t-1]",
    paste0("    lambda0[t] <- ", off_str, "inprod(Xc[t,1:Kc], beta[1:Kc])"),
    "  }",
    "",
    "  # Latent self-excitation indicator Z with Markov-like structure",
    "  Z[1] ~ dbern(p[1])",
    "  p[1] ~ dunif(0,1)",
    "  for(t in 2:N){",
    "    Z[t] ~ dbern(p[t])",
    "    p[t] <- p[1] * (p11 + p00 - 1)^(t-1) + (1 - p00) * ((1 - (p11 + p00 - 1)^(t-1)) / (2 - p11 - p00))",
    "  }",
    "",
    "  # Priors",
    paste0("  r ~ dgamma(", r_prior_shape, ", ", r_prior_rate, ")"),
    "  eta ~ dbeta(1,1)",
    paste0("  p00 ~ dbeta(", p_priors, ", ", p_priors, ")"),
    paste0("  p11 ~ dbeta(", p_priors, ", ", p_priors, ")"),
    paste0("  for(k in 1:Kc){ beta[k]  ~ dnorm(", beta_prior_mean, ", 1/(", beta_prior_sd, "^2)) }"),
    "}",
    ""
  )
  model_string <- paste(model_lines, collapse = "\n")

  model_file <- tempfile(fileext = ".bug")
  writeLines(model_string, model_file)
  on.exit(unlink(model_file), add = TRUE)

  # Initial values
  if (is.null(beta_init))  beta_init  <- lapply(1:n_chains, function(i) rep(0, Kc))
  if (is.null(r_init))     r_init     <- seq(0.5, 0.5 + 0.5*(n_chains - 1), length.out = n_chains)

  inits <- lapply(1:n_chains, function(i) list(
    beta  = beta_init[[i]],
    r     = r_init[i]
  ))

  data4Jags <- list(
    Y  = cases,
    N  = N,
    Xc = Xc,
    pop = pop_vec,
    Kc = Kc
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


  full_summary <- as.data.frame(jags.out$BUGSoutput$summary)
  full_summary$dic <- jags.out$BUGSoutput$DIC

  # Filter out Z[...] rows for the main summary view
  summary_df <- full_summary[!grepl("^Z\\[", rownames(full_summary)), , drop = FALSE]


  # WAIC attempt (fallback-safe)
  s <- tryCatch({
    rjags::jags.samples(jags.out$model, c("WAIC", "deviance"), type = "mean", n.iter = 1000)
  }, error = function(e) NULL)

  if (!is.null(s) && !is.null(s$WAIC) && !is.null(s$deviance)) {
    p_waic <- sum(s$WAIC)
    dev    <- sum(s$deviance)
    waic_vals <- round(c(waic = dev + p_waic, p_waic = p_waic), 1)
  } else {
    waic_vals <- c(waic = NA, p_waic = NA)
  }




  ret <- list(
    mcmc_summary      = summary_df,
    mcmc_summary_full = full_summary,
    dic               = summary_df$dic[1],
    waic              = waic_vals,
    raw_output        = jags.out
  )

  # Optional Z plot
  if (plot_Z) {
    # Extract posterior mean Z
    Z_mean <- jags.out$BUGSoutput$mean$Z
    # Prepare date/index
    if (is.null(dates)) {
      dates_plot <- seq_len(length(Z_mean))
    } else {
      # try to coerce to Date, fallback to raw
      safe_date <- try(as.Date(dates), silent = TRUE)
      if (!inherits(safe_date, "try-error") && length(safe_date) == length(Z_mean)) {
        dates_plot <- safe_date
      } else if (length(dates) == length(Z_mean)) {
        dates_plot <- dates
      } else {
        stop("Provided 'dates' must match length of Z (number of timepoints).")
      }
    }

    # Build plotting data
    if (is.null(dim(Z_mean)) || length(dim(Z_mean)) == 1) {
      df_plot <- data.frame(date = dates_plot, value = as.numeric(Z_mean))
    } else {
      d <- reshape2::melt(Z_mean)
      # assume first dimension is time index
      if (is.numeric(d[[1]]) && length(dates_plot) >= max(as.integer(d[[1]]))) {
        d$date <- dates_plot[as.integer(d[[1]])]
      } else {
        d$date <- dates_plot
      }
      df_plot <- data.frame(date = d$date, value = d$value)
    }

    sp <- ggplot2::ggplot(df_plot, ggplot2::aes(x = date, y = value)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Posterior Mean of Latent Z over Time",
                    x = if (!is.null(dates)) "Date" else "Index",
                    y = expression(E[Z]))
    ret$plot_Z <- sp

    return(ret)
  }
}
