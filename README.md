sparsesurv
================

# sparsesurv: Forecasting and Early Outbreak Detection for Sparse Count Data

**sparsesurv** provides tools for fitting, forecasting, and prospective
outbreak detection in sparse surveillance count time series. Models are
built around negative-binomial and zero-modified families (ZINB, hurdle)
with optional self-excitation and GARMA dynamics.

A key feature is **forecasting ahead using future climatic/environmental
covariates** you supply (e.g., weather or seasonality drivers). The
package also **separates endemic and epidemic states** and returns, for
each time point, the **probability of being in an epidemic state**,
enabling transparent flagging of emerging outbreaks.

**sparsesurv** provides tools for fitting, forecasting, and flagging
outbreaks in sparse surveillance count time series. It focuses on
negative-binomial families and zero-modified extensions common in
public-health surveillance, with optional self-excitation and
GARMA-style dynamics.

Potential users include biostatisticians, epidemiologists, and others
working with surveillance data (also useful wherever counts are
sparse/overdispersed).

## Features

- **Forecast with future covariates**: plug in known/forecasted
  climate/environmental inputs to get out-of-sample case forecasts with
  prediction intervals. (NB,SENB, GARMA_NB, ZINB, SEZINB, GARMA_ZINB)

- **Endemic vs. epidemic decomposition**: per-time-point posterior
  probability of the epidemic state, plus utilities to flag episodes by
  thresholding. (outbreakNB, outbreakZINB, outbreakHNB)

- **Model families**

  - Negative binomial (NB)
  - Self-exciting NB (SENB)
  - Generalised Autoregressive Moving Average NB (GARMA_NB)
  - NB Hidden Markov (outbreakNB)
  - Zero-inflated NB (ZINB)
  - Self-exciting ZINB (SEZINB)
  - Generalised Autoregressive Moving Average ZINB (GARMA_ZINB)
  - ZINB Hidden Markov (outbreakZINB)
  - Hurdle NB Hidden Markov (outbreakHNB)

- **Covariates**

  - Climatic/environmental covariates in the mean model and/or zero
    component

- **Forecasting**

  - Posterior predictive forecasts for near-term horizons (using
    NB,SENB, GARMA_NB, ZINB, SEZINB, GARMA_ZINB)

- **Outbreak detection**

  - Prospective detection rules for NB / ZINB / hurdle models (using
    outbreakNB, outbreakZINB, outbreakHNB)

- **Diagnostics** -summaries and residual checks to assess fit and
  calibration.

> System requirements: **JAGS (\>= 4.x)** for Bayesian fitting via
> *R2jags*.

## Installation

The stable release version of sparsesurv is hosted on the Comprehensive
R Archive Network (CRAN) at
<https://CRAN.R-project.org/package=surveillance> and can be installed
via

``` r

install.packages("sparsesurv")
```

Development version (GitHub):

``` r
# install.packages("remotes")
remotes::install_github("alexangelakis-ang/sparsesurv")
```

> Windows/macOS users: please install **JAGS** first so model fitting
> works.

## Quick start

Below is a lightweight example using the self-exciting NB fitter
`SENB()` on a toy series. (Kept short for README speed; increase
iterations/chains in real analyses.)

``` r
set.seed(1)
# Toy series: NB counts with mean ~8 and size ~5
cases <- rnbinom(72, size = 5, mu = 8)

# Fit a compact model (short MCMC for demonstration)
fit_nb <- SENB(
  cases = cases,
  beta_prior_mean = 0,
  beta_prior_sd   = 5,
  r_prior_shape   = 2,
  r_prior_rate    = 0.5,
  n_iter  = 400,
  n_burnin= 200,
  n_chains= 1,
  n_thin  = 2
)

print(fit_nb)
```

``` r
set.seed(1)
# Toy series: NB counts with mean ~8 and size ~5
cases <- rnbinom(72, size = 5, mu = 8)

# Fit a compact model (short MCMC for demonstration)

fit_outbreakzinb <- outbreakZINB(
  cases = cases,
  beta_prior_mean = 0,
  beta_prior_sd = 10,
  r_prior_shape = 2,
  r_prior_rate = 0.5,
  n_iter = 1000,
  n_burnin = 100,
  n_chains = 2,
  n_thin = 1,
  dates = dat1$date,
  plot_Z = TRUE
)
print(fit_outbreakzinb)
```

## Documentation

- Detailed help for each function can be found in the package
  documentation (e.g., ?SENB after installation) or the references.

## Background & references

- Early detection of outbreaks:

- Forecasting:

## Contributing

Issues and pull requests:
<https://github.com/alexangelakis-ang/sparsesurv/issues> or via e-mail
to maintainer(“sparsesurv”).

Please include a minimal reproducible example for bugs. For large
features, open an issue to discuss first.

## Citing

If you use **sparsesurv** in academic work, please cite the package.

## License

The sparsesurv package is free and open-source software, and you are
welcome to redistribute it under the terms of the GNU General Public
License, version 3. This program is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY.

## Funding

This work is supported by European Union’s Horizon 2020 research and
innovation programme under grant agreement No 101000365, project
PREPARE4VBD (A Cross- Disciplinary Alliance to Identify, PREdict and
prePARE for Emerging Vector-Borne Diseases).
