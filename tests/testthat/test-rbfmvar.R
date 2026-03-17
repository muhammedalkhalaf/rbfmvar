test_that("rbfmvar returns expected structure", {
  set.seed(1)
  n <- 80
  y1 <- cumsum(rnorm(n))
  y2 <- 0.5 * y1 + rnorm(n, sd = 0.5)
  Y  <- cbind(y1, y2)
  res <- rbfmvar(Y, lags = 2)
  expect_s3_class(res, "rbfmvar")
  expect_equal(res$n_vars, 2)
  expect_equal(res$nobs, n)
  expect_true(is.matrix(res$F_plus))
  expect_true(is.matrix(res$Sigma_e))
})

test_that("rbfmvar Bartlett kernel works", {
  set.seed(2)
  Y <- cbind(cumsum(rnorm(60)), cumsum(rnorm(60)))
  res <- rbfmvar(Y, lags = 2, kernel = "bartlett")
  expect_true(res$bandwidth > 0)
})

test_that("rbfmvar Parzen kernel works", {
  set.seed(3)
  Y <- cbind(cumsum(rnorm(60)), cumsum(rnorm(60)))
  res <- rbfmvar(Y, lags = 2, kernel = "parzen")
  expect_s3_class(res, "rbfmvar")
})

test_that("rbfmvar lag selection via AIC", {
  set.seed(4)
  Y <- cbind(cumsum(rnorm(80)), cumsum(rnorm(80)))
  res <- rbfmvar(Y, ic = "aic", max_lags = 4)
  expect_true(res$p_lags >= 2)
})

test_that("rbfmvar Granger test runs", {
  set.seed(5)
  n <- 80
  y1 <- cumsum(rnorm(n))
  y2 <- 0.5 * y1 + rnorm(n)
  Y  <- cbind(y1 = y1, y2 = y2)
  res <- rbfmvar(Y, lags = 2, granger = "y1 -> y2")
  expect_false(is.null(res$granger))
  expect_true(is.numeric(res$granger$wald_stat))
})

test_that("rbfmvar input validation", {
  expect_error(rbfmvar(matrix(rnorm(30), ncol = 1)))
  expect_error(rbfmvar(matrix(rnorm(30), ncol = 3)))
})

test_that("print.rbfmvar produces output", {
  set.seed(7)
  Y <- cbind(cumsum(rnorm(60)), cumsum(rnorm(60)))
  res <- rbfmvar(Y, lags = 2)
  out <- capture.output(print(res))
  expect_true(length(out) > 0)
})
