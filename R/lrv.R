#' @title Long-Run Variance Estimation
#' @name lrv
#' @description Functions for estimating long-run variance (LRV) matrices using
#'   kernel-based methods with automatic bandwidth selection.
#'
#' @references
#' Andrews, D. W. K. (1991). Heteroskedasticity and Autocorrelation Consistent
#' Covariance Matrix Estimation. \emph{Econometrica}, 59(3), 817-858.
#' \doi{10.2307/2938229}
#'
#' Newey, W. K., & West, K. D. (1994). Automatic Lag Selection in Covariance
#' Matrix Estimation. \emph{Review of Economic Studies}, 61(4), 631-653.
#' \doi{10.2307/2297912}
#'
#' @keywords internal
NULL

#' Estimate Long-Run Variance Matrix
#'
#' @description Estimates the long-run variance matrix of a multivariate time
#'   series using kernel-weighted autocovariances.
#'
#' @param v Matrix of residuals (T x n).
#' @param kernel Character string specifying the kernel type:
#'   \code{"bartlett"}, \code{"parzen"}, or \code{"qs"}.
#' @param bandwidth Bandwidth parameter. If \code{-1} (default), automatic
#'   bandwidth selection via Andrews (1991) is used.
#' @param prewhiten Logical. Whether to prewhiten the series before LRV
#'   estimation. Default is \code{FALSE}.
#'
#' @return A list containing:
#'   \describe{
#'     \item{Omega}{Estimated long-run variance matrix (n x n).}
#'     \item{bandwidth}{Bandwidth used in estimation.}
#'     \item{kernel}{Kernel used.}
#'   }
#'
#' @references
#' Andrews, D. W. K. (1991). Heteroskedasticity and Autocorrelation Consistent
#' Covariance Matrix Estimation. \emph{Econometrica}, 59(3), 817-858.
#' \doi{10.2307/2938229}
#'
#' @keywords internal
estimate_lrv <- function(v, kernel = "bartlett", bandwidth = -1,
                         prewhiten = FALSE) {
  v <- as.matrix(v)
  TT <- nrow(v)
  n <- ncol(v)

  # Center the series
  v <- scale(v, center = TRUE, scale = FALSE)

  # Automatic bandwidth selection
  if (bandwidth < 0) {
    bandwidth <- select_bandwidth_andrews(v, kernel)
  }

  # Get kernel function
  kern_fn <- get_kernel_function(kernel)

  # Compute autocovariances and weighted sum
  Omega <- matrix(0, n, n)
  max_lag <- .kernel_max_lag(kernel, bandwidth, TT)

  for (j in 0:max_lag) {
    # Autocovariance at lag j
    if (j == 0) {
      Gamma_j <- crossprod(v) / TT
    } else {
      Gamma_j <- crossprod(v[(j + 1):TT, , drop = FALSE],
                           v[1:(TT - j), , drop = FALSE]) / TT
    }

    # Kernel weight
    w <- kern_fn(j / bandwidth)

    if (j == 0) {
      Omega <- Omega + w * Gamma_j
    } else {
      # Add both Gamma_j and Gamma_j' for symmetry
      Omega <- Omega + w * (Gamma_j + t(Gamma_j))
    }
  }

  # Ensure positive semi-definiteness
  Omega <- (Omega + t(Omega)) / 2
  eig <- eigen(Omega, symmetric = TRUE)
  eig$values <- pmax(eig$values, 0)
  Omega <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)

  list(
    Omega = Omega,
    bandwidth = bandwidth,
    kernel = kernel
  )
}

#' Andrews (1991) Automatic Bandwidth Selection
#'
#' @description Implements the data-dependent automatic bandwidth selection
#'   procedure of Andrews (1991).
#'
#' @param v Matrix of residuals (T x n).
#' @param kernel Character string specifying the kernel type.
#'
#' @return Optimal bandwidth.
#'
#' @references
#' Andrews, D. W. K. (1991). Heteroskedasticity and Autocorrelation Consistent
#' Covariance Matrix Estimation. \emph{Econometrica}, 59(3), 817-858.
#' \doi{10.2307/2938229}
#'
#' @keywords internal
select_bandwidth_andrews <- function(v, kernel = "bartlett") {
  v <- as.matrix(v)
  TT <- nrow(v)
  n <- ncol(v)

  q <- get_kernel_exponent(kernel)
  c_gamma <- get_kernel_constant(kernel)

  # Fit AR(1) to each series and estimate alpha
  alpha_num <- 0
  alpha_den <- 0

  for (i in 1:n) {
    vi <- v[, i]

    # AR(1) regression
    y <- vi[2:TT]
    x <- vi[1:(TT - 1)]

    rho <- sum(x * y) / sum(x^2)
    sigma2 <- sum((y - rho * x)^2) / (TT - 2)

    # Bound rho away from unit root
    rho <- min(max(rho, -0.99), 0.99)

    if (q == 1) {
      # Bartlett
      alpha_num <- alpha_num + 4 * rho^2 * sigma2^2 / ((1 - rho)^6 * (1 + rho)^2)
      alpha_den <- alpha_den + sigma2^2 / ((1 - rho)^4)
    } else {
      # Parzen/QS (q = 2)
      alpha_num <- alpha_num + 4 * rho^2 * sigma2^2 / ((1 - rho)^8)
      alpha_den <- alpha_den + sigma2^2 / ((1 - rho)^4)
    }
  }

  alpha <- alpha_num / alpha_den

  # Optimal bandwidth
  bw <- c_gamma * (alpha * TT)^(1 / (2 * q + 1))

  # Bound bandwidth
  bw <- max(1, min(bw, TT / 2))

  bw
}

#' Kernel Lag Truncation
#'
#' Bartlett and Parzen weights vanish beyond the bandwidth; the Quadratic
#' Spectral kernel has unbounded support, so all lags are used.
#' @keywords internal
.kernel_max_lag <- function(kernel, bandwidth, TT) {
  if (tolower(kernel) == "qs") TT - 1 else min(TT - 1, ceiling(bandwidth))
}

#' Estimate Two-Sided Long-Run Covariance
#'
#' @description \eqn{\hat\Omega_{ab} = \sum_{|j| < T} k(j/K)\hat\Gamma_{ab}(j)}
#'   with \eqn{\hat\Gamma_{ab}(j) = T^{-1}\sum_t a_{t+j} b_t'}.
#'
#' @param a Matrix (T x n1).
#' @param b Matrix (T x n2).
#' @param kernel Character string specifying the kernel type.
#' @param bandwidth Bandwidth parameter.
#'
#' @return Long-run covariance matrix (n1 x n2).
#' @keywords internal
estimate_cross_lrv <- function(a, b, kernel = "bartlett", bandwidth) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  TT <- nrow(a)
  kern_fn <- get_kernel_function(kernel)
  Omega <- crossprod(a, b) / TT
  for (j in seq_len(.kernel_max_lag(kernel, bandwidth, TT))) {
    w <- kern_fn(j / bandwidth)
    if (w == 0) next
    G_pos <- crossprod(a[(j + 1):TT, , drop = FALSE], b[1:(TT - j), , drop = FALSE]) / TT
    G_neg <- crossprod(a[1:(TT - j), , drop = FALSE], b[(j + 1):TT, , drop = FALSE]) / TT
    Omega <- Omega + w * (G_pos + G_neg)
  }
  Omega
}

#' Estimate One-Sided Long-Run Covariance
#'
#' @description \eqn{\hat\Delta_{ev} = \sum_{j=0}^{T-1} k(j/K) T^{-1}
#'   \sum_t e_{t+j} v_t'}, the kernel estimate of \eqn{\sum_{j \ge 0} E
#'   e_j v_0'} (the convention of Chang, 2000, p. 909).
#'
#' @param e Matrix of residuals e (T x n1).
#' @param v Matrix of residuals v (T x n2).
#' @param kernel Character string specifying the kernel type.
#' @param bandwidth Bandwidth parameter.
#'
#' @return One-sided long-run covariance matrix (n1 x n2).
#'
#' @keywords internal
estimate_onesided_lrv <- function(e, v, kernel = "bartlett", bandwidth) {
  e <- as.matrix(e)
  v <- as.matrix(v)
  TT <- nrow(e)
  kern_fn <- get_kernel_function(kernel)
  Delta <- crossprod(e, v) / TT
  for (j in seq_len(.kernel_max_lag(kernel, bandwidth, TT))) {
    w <- kern_fn(j / bandwidth)
    if (w == 0) next
    Delta <- Delta + w * crossprod(e[(j + 1):TT, , drop = FALSE],
                                   v[1:(TT - j), , drop = FALSE]) / TT
  }
  Delta
}
