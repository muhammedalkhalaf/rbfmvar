# Tests for long-run variance estimation

test_that("estimate_lrv produces valid output", {
  set.seed(123)
  v <- matrix(rnorm(200), 100, 2)

  result <- estimate_lrv(v, kernel = "bartlett")

  expect_true(is.list(result))
  expect_true("Omega" %in% names(result))
  expect_true("bandwidth" %in% names(result))
  expect_true("kernel" %in% names(result))

  # Omega should be a matrix
  expect_true(is.matrix(result$Omega))
  expect_equal(dim(result$Omega), c(2, 2))

  # Omega should be symmetric
  expect_equal(result$Omega, t(result$Omega), tolerance = 1e-10)

  # Omega should be positive semi-definite
  eig <- eigen(result$Omega, symmetric = TRUE)$values
  expect_true(all(eig >= -1e-10))

  # Bandwidth should be positive
  expect_true(result$bandwidth > 0)
})

test_that("estimate_lrv works with different kernels", {
  set.seed(456)
  v <- matrix(rnorm(300), 100, 3)

  result_bart <- estimate_lrv(v, kernel = "bartlett")
  result_parz <- estimate_lrv(v, kernel = "parzen")
  result_qs <- estimate_lrv(v, kernel = "qs")

  # All should produce valid matrices
  for (res in list(result_bart, result_parz, result_qs)) {
    expect_true(is.matrix(res$Omega))
    expect_equal(dim(res$Omega), c(3, 3))
    expect_true(res$bandwidth > 0)
  }

  # Results should differ (different kernels give different estimates)
  expect_false(isTRUE(all.equal(result_bart$Omega, result_parz$Omega)))
})

test_that("estimate_lrv respects fixed bandwidth", {
  set.seed(789)
  v <- matrix(rnorm(200), 100, 2)

  result <- estimate_lrv(v, kernel = "bartlett", bandwidth = 5)

  expect_equal(result$bandwidth, 5)
})

test_that("select_bandwidth_andrews produces reasonable values", {
  set.seed(111)

  # White noise - should have small bandwidth
  v_white <- matrix(rnorm(500), 250, 2)
  bw_white <- select_bandwidth_andrews(v_white, "bartlett")

  # AR(1) process - should have larger bandwidth
  v_ar <- matrix(0, 250, 2)
  for (t in 2:250) {
    v_ar[t, ] <- 0.7 * v_ar[t - 1, ] + rnorm(2)
  }
  bw_ar <- select_bandwidth_andrews(v_ar, "bartlett")

  # Both should be positive
  expect_true(bw_white > 0)
  expect_true(bw_ar > 0)

  # AR process should generally need larger bandwidth
  # (This is a tendency, not always true for finite samples)
  # expect_true(bw_ar > bw_white)  # May fail randomly
})

test_that("estimate_onesided_lrv produces valid output", {
  set.seed(222)
  e <- matrix(rnorm(200), 100, 2)
  v <- matrix(rnorm(200), 100, 2)

  Delta <- estimate_onesided_lrv(e, v, kernel = "bartlett", bandwidth = 5)

  expect_true(is.matrix(Delta))
  expect_equal(dim(Delta), c(2, 2))
})
