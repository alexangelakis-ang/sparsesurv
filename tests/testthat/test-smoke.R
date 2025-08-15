test_that("package loads and exports exist", {
  expect_s3_class(packageVersion("sparsesurv"), "package_version")
  expect_true(length(getNamespaceExports("sparsesurv")) >= 1)
})

test_that("skip heavy/JAGS tests on CRAN", {
  skip_on_cran()
  skip_if_not_installed("R2jags")
  # placeholder: put a tiny call that exercises a fast path if you have one
  # e.g. expect_no_error(fit_nb_mock(...))
  succeed()
})

