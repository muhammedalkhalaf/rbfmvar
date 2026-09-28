## Checks of the RBFM-VAR estimator against the formulas of Chang (2000),
## Econometric Theory 16(6), 905-926, equations (10)-(13) and (20).

sim_chang <- function(TT, r1, r2, seed) {
  set.seed(seed)
  S <- matrix(c(1, .5, .5, 1), 2)
  e <- matrix(stats::rnorm(2 * (TT + 50)), ncol = 2) %*% chol(S)
  y <- matrix(0, TT + 50, 2)
  for (t in 3:(TT + 50)) {
    dy1 <- r1 * (y[t - 1, 1] - y[t - 2, 1]) +
      r2 * (y[t - 1, 1] - (y[t - 1, 2] - y[t - 2, 2])) + e[t, 1]
    y[t, 1] <- y[t - 1, 1] + dy1
    y[t, 2] <- 2 * y[t - 1, 2] - y[t - 2, 2] + e[t, 2]
  }
  y <- y[-(1:50), ]
  colnames(y) <- c("y1", "y2")
  y
}

test_that("F+ equals (Y'Z, Y+'W + T Delta+)(X'X)^{-1} computed directly", {
  y <- sim_chang(150, -0.3, -0.15, 11)
  f <- rbfmvar(y, lags = 3, bandwidth = 4)
  TT <- nrow(y)
  rows <- 4:TT
  d1 <- rbind(NA, diff(y)); d2 <- rbind(NA, diff(d1))
  L <- function(M, j) rbind(matrix(NA, j, ncol(M)), M[1:(nrow(M) - j), , drop = FALSE])
  Y <- d2[rows, ]
  Z <- L(d2, 1)[rows, ]
  W <- cbind(L(d1, 1)[rows, ], L(y, 1)[rows, ])
  X <- cbind(Z, W)
  e <- Y - X %*% solve(crossprod(X), crossprod(X, Y))
  dy1 <- L(d1, 1)[rows, ]; dy2 <- L(d1, 2)[rows, ]
  N <- t(solve(crossprod(dy2), crossprod(dy2, dy1)))
  v <- cbind(L(d2, 1)[rows, ], dy1 - dy2 %*% t(N))
  dW <- cbind(L(d2, 1)[rows, ], dy1)
  Tn <- length(rows); K <- 4
  k <- function(x) pmax(1 - abs(x), 0)
  G <- function(a, b, j) crossprod(a[(j + 1):Tn, , drop = FALSE], b[1:(Tn - j), , drop = FALSE]) / Tn
  Om <- function(a, b) Reduce(`+`, lapply(1:(Tn - 1), function(j) k(j / K) * (G(a, b, j) + t(G(b, a, j)))),
                              crossprod(a, b) / Tn)
  Oev <- Om(e, v); Ovv <- Om(v, v)
  Dvw <- Reduce(`+`, lapply(1:(Tn - 1), function(j) k(j / K) * G(v, dW, j)), crossprod(v, dW) / Tn)
  B <- Oev %*% solve((Ovv + t(Ovv)) / 2)
  Yp <- Y - v %*% t(B)
  Fp <- cbind(t(Y) %*% Z, t(Yp) %*% W + Tn * B %*% Dvw) %*% solve(crossprod(X))
  expect_equal(unname(t(f$F_plus)), unname(Fp), tolerance = 1e-10)
})

test_that("levels VAR implied by the estimates reproduces the fitted values", {
  y <- sim_chang(150, 0.5, 0, 3)
  f <- rbfmvar(y, lags = 3)
  A <- rbfmvar:::.rbfmvar_levels(f$F_plus, 2, 3)
  t <- nrow(y)
  lv <- A[[1]] %*% y[t - 1, ] + A[[2]] %*% y[t - 2, ] + A[[3]] %*% y[t - 3, ]
  expect_equal(as.numeric(lv - 2 * y[t - 1, ] + y[t - 2, ]),
               as.numeric(f$fitted[f$T_eff, ]), tolerance = 1e-8)
})

test_that("Granger test uses the full covariance of the restricted block", {
  y <- sim_chang(150, -0.3, -0.15, 5)
  f <- rbfmvar(y, lags = 2)
  g <- granger_test(f, "y2", "y1")
  i <- grep("\\.y2$", f$regnames)
  b <- f$F_plus[i, 1]
  V <- f$Sigma_e[1, 1] * f$XtX_inv[i, i]
  expect_equal(g$statistic, as.numeric(t(b) %*% solve(V, b)))
  expect_equal(g$df, 2L)
})

test_that("IRF bootstrap re-estimates the model and forecasts are in levels", {
  y <- sim_chang(120, 0.5, 0, 9)
  f <- rbfmvar(y, lags = 2)
  ir <- irf(f, horizon = 4, boot = 5, seed = 1)
  expect_equal(dim(ir$irf_lower), c(5, 2, 2))
  expect_true(all(ir$irf_lower <= ir$irf_upper + 1e-12))
  fc <- forecast(f, h = 2)
  A <- rbfmvar:::.rbfmvar_levels(f$F_plus, 2, 2)
  t <- nrow(y)
  expect_equal(as.numeric(fc$mean[, 1]),
               as.numeric(A[[1]] %*% y[t, ] + A[[2]] %*% y[t - 1, ]))
  expect_equal(as.numeric(fc$se[, 1]), sqrt(diag(f$Sigma_e)), ignore_attr = TRUE)
})
