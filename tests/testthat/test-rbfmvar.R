# Tests for rbfmvar package

test_that("rbfmvar estimates a basic model", {
  set.seed(123)
  n <- 100
  e <- matrix(rnorm(n * 3), n, 3)
  y <- matrix(0, n, 3)
  colnames(y) <- c("y1", "y2", "y3")

  for (t in 3:n) {
    y[t, ] <- 0.3 * y[t - 1, ] + 0.2 * y[t - 2, ] + e[t, ]
  }

  fit <- rbfmvar(y, lags = 2)

  expect_s3_class(fit, "rbfmvar")
  expect_equal(fit$n_vars, 3)
  expect_equal(fit$p_lags, 2)
  expect_true(fit$T_eff > 0)
  expect_true(fit$bandwidth > 0)
})

test_that("rbfmvar validates inputs correctly", {
  # Less than 2 variables
  expect_error(rbfmvar(matrix(rnorm(100), 100, 1)), "At least 2 variables")

  # Invalid lags
  expect_error(rbfmvar(matrix(rnorm(200), 100, 2), lags = 0), "at least 1")

  # Invalid kernel
  expect_error(rbfmvar(matrix(rnorm(200), 100, 2), kernel = "invalid"),
               "bartlett.*parzen.*qs")

  # Invalid IC
  expect_error(rbfmvar(matrix(rnorm(200), 100, 2), ic = "invalid"),
               "aic.*bic.*hq.*none")

  # Insufficient observations
  expect_error(rbfmvar(matrix(rnorm(30), 15, 2)), "Insufficient observations")
})

test_that("different kernels work", {
  set.seed(456)
  y <- matrix(rnorm(200), 100, 2)
  colnames(y) <- c("x", "y")

  fit_bart <- rbfmvar(y, lags = 2, kernel = "bartlett")
  fit_parz <- rbfmvar(y, lags = 2, kernel = "parzen")
  fit_qs <- rbfmvar(y, lags = 2, kernel = "qs")

  expect_equal(fit_bart$kernel, "bartlett")
  expect_equal(fit_parz$kernel, "parzen")
  expect_equal(fit_qs$kernel, "qs")

  # Bandwidths should be positive
  expect_true(fit_bart$bandwidth > 0)
  expect_true(fit_parz$bandwidth > 0)
  expect_true(fit_qs$bandwidth > 0)
})

test_that("lag selection via IC works", {
  skip_on_cran()  # May be slow

  set.seed(789)
  n <- 150
  e <- matrix(rnorm(n * 2), n, 2)
  y <- matrix(0, n, 2)
  colnames(y) <- c("x", "y")

  for (t in 4:n) {
    y[t, ] <- 0.5 * y[t - 1, ] + 0.2 * y[t - 2, ] + e[t, ]
  }

  fit_aic <- suppressMessages(rbfmvar(y, max_lags = 4, ic = "aic"))
  fit_bic <- suppressMessages(rbfmvar(y, max_lags = 4, ic = "bic"))

  expect_true(fit_aic$p_lags >= 1)
  expect_true(fit_bic$p_lags >= 1)
  expect_equal(fit_aic$ic, "aic")
  expect_equal(fit_bic$ic, "bic")
})

test_that("summary method works", {
  set.seed(111)
  y <- matrix(rnorm(200), 100, 2)
  colnames(y) <- c("x", "y")

  fit <- rbfmvar(y, lags = 2)
  summ <- summary(fit)

  expect_s3_class(summ, "summary.rbfmvar")
  expect_true("coefficients" %in% names(summ))
  expect_equal(length(summ$coefficients), 2)
})

test_that("coef, residuals, fitted methods work", {
  set.seed(222)
  y <- matrix(rnorm(200), 100, 2)
  colnames(y) <- c("x", "y")

  fit <- rbfmvar(y, lags = 2)

  cf <- coef(fit)
  expect_true(is.matrix(cf))

  res <- residuals(fit)
  expect_true(is.matrix(res))
  expect_equal(ncol(res), 2)

  ft <- fitted(fit)
  expect_true(is.matrix(ft))
  expect_equal(ncol(ft), 2)
})

test_that("Granger causality test works", {
  set.seed(333)
  n <- 150
  e <- matrix(rnorm(n * 2), n, 2)
  y <- matrix(0, n, 2)
  colnames(y) <- c("x", "y")

  # x causes y
  for (t in 3:n) {
    y[t, "x"] <- 0.5 * y[t - 1, "x"] + e[t, 1]
    y[t, "y"] <- 0.3 * y[t - 1, "y"] + 0.4 * y[t - 1, "x"] + e[t, 2]
  }

  fit <- rbfmvar(y, lags = 2)
  test <- granger_test(fit, cause = "x", effect = "y")

  expect_s3_class(test, "rbfmvar_granger")
  expect_true(test$statistic >= 0)
  expect_true(test$p.value >= 0 && test$p.value <= 1)
  expect_true(test$df >= 1)
})

test_that("IRF computation works", {
  set.seed(444)
  y <- matrix(rnorm(200), 100, 2)
  colnames(y) <- c("x", "y")

  fit <- rbfmvar(y, lags = 2)
  ir <- irf(fit, horizon = 10)

  expect_s3_class(ir, "rbfmvar_irf")
  expect_equal(ir$horizon, 10)
  expect_true(is.array(ir$irf))
  expect_equal(dim(ir$irf)[1], 11)  # horizon + 1
})

test_that("FEVD computation works", {
  set.seed(555)
  y <- matrix(rnorm(200), 100, 2)
  colnames(y) <- c("x", "y")

  fit <- rbfmvar(y, lags = 2)
  fv <- fevd(fit, horizon = 10)

  expect_s3_class(fv, "rbfmvar_fevd")
  expect_equal(fv$horizon, 10)
  expect_true(is.array(fv$fevd))

  # FEVD should sum to 1 for each variable at each horizon
  for (h in 1:11) {
    for (i in 1:2) {
      total <- sum(fv$fevd[h, i, ])
      expect_equal(total, 1, tolerance = 1e-6)
    }
  }
})

test_that("forecast method works", {
  set.seed(666)
  y <- matrix(rnorm(200), 100, 2)
  colnames(y) <- c("x", "y")

  fit <- rbfmvar(y, lags = 2)
  fc <- forecast(fit, h = 5)

  expect_s3_class(fc, "rbfmvar_forecast")
  expect_equal(fc$horizon, 5)
  expect_true(is.matrix(fc$mean))
  expect_equal(ncol(fc$mean), 5)
  expect_equal(nrow(fc$mean), 2)

  # Lower bound should be less than upper
  expect_true(all(fc$lower <= fc$upper))
})

test_that("Granger matrix works", {
  set.seed(777)
  y <- matrix(rnorm(300), 100, 3)
  colnames(y) <- c("x", "y", "z")

  fit <- rbfmvar(y, lags = 2)
  gm <- granger_matrix(fit)

  expect_true(is.matrix(gm))
  expect_equal(nrow(gm), 3)
  expect_equal(ncol(gm), 3)

  # Diagonal should be NA
  expect_true(all(is.na(diag(gm))))

  # Off-diagonal should be valid p-values
  off_diag <- gm[!is.na(gm)]
  expect_true(all(off_diag >= 0 & off_diag <= 1))
})
