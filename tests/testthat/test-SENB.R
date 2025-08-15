test_that("SENB runs on toy data (local only)", {
  skip_on_cran()
  skip_if_not_installed("R2jags")

  # Skip if JAGS binary not found
  has_jags <- nzchar(Sys.which("jags")) ||
    nzchar(Sys.which("jags-terminal")) ||
    nzchar(Sys.which("jags.exe"))
  skip_if(!has_jags, "JAGS not installed")

  set.seed(1)
  cases <- rnbinom(72, size = 5, mu = 8)

  fit_nb <- NULL
  expect_error(
    fit_nb <- SENB(
      cases = cases,
      beta_prior_mean = 0,
      beta_prior_sd   = 5,
      r_prior_shape   = 2,
      r_prior_rate    = 0.5,
      n_iter  = 400,   # keep fast for CI
      n_burnin= 200,
      n_chains= 1,
      n_thin  = 2
    ),
    NA
  )

  expect_true(!is.null(fit_nb))
  # minimal structure check (adjust if you have a class)
  expect_true(is.list(fit_nb) || length(class(fit_nb)) > 0)
})
