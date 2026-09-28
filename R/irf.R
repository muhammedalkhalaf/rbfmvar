#' @title Impulse Response Functions
#'
#' @description Computes orthogonalized impulse response functions (IRF) from
#'   an RBFM-VAR model with optional bootstrap confidence intervals.
#'
#' @param object An \code{rbfmvar} object from \code{\link{rbfmvar}}.
#' @param horizon Integer. Number of periods for the IRF. Default is 20.
#' @param ortho Logical. If \code{TRUE} (default), compute orthogonalized IRFs
#'   using Cholesky decomposition of the error covariance matrix.
#' @param boot Integer. Number of bootstrap replications for confidence
#'   intervals. If 0 (default), no bootstrap is performed.
#' @param ci Numeric. Confidence level for bootstrap intervals (0-100).
#'   Default is 90.
#' @param seed Integer. Random seed for reproducibility. Default is \code{NULL}.
#'
#' @details
#' The IRF measures the response of each variable to a one-standard-deviation
#' shock in each of the structural innovations. When \code{ortho = TRUE},
#' the structural shocks are identified using the Cholesky decomposition of
#' the residual covariance matrix (recursive identification).
#'
#' The responses are those of the levels VAR implied by the RBFM-VAR
#' coefficients (see the internal function \code{.rbfmvar_levels}). With
#' \code{ortho = FALSE} the responses are to unit shocks in the reduced-form
#' errors. Bootstrap intervals are percentile intervals from a
#' recursive-design residual bootstrap: each replication rebuilds the sample
#' from the estimated levels VAR, the observed initial values and resampled
#' centred residuals, and re-estimates the RBFM-VAR with the same lag order,
#' kernel and bandwidth. No bias correction is applied.
#'
#' @return An object of class \code{"rbfmvar_irf"} containing:
#' \describe{
#'   \item{irf}{Array of IRF values (horizon x n x n). Element [h, i, j] is
#'     the response of variable i to a shock in variable j at horizon h.}
#'   \item{irf_lower}{Lower confidence bounds (if bootstrap was performed).}
#'   \item{irf_upper}{Upper confidence bounds (if bootstrap was performed).}
#'   \item{horizon}{IRF horizon.}
#'   \item{varnames}{Variable names.}
#'   \item{ortho}{Whether orthogonalized IRFs were computed.}
#'   \item{boot}{Number of bootstrap replications.}
#'   \item{ci}{Confidence level.}
#' }
#'
#' @references
#' Lutkepohl, H. (2005). \emph{New Introduction to Multiple Time Series Analysis}.
#' Springer-Verlag. \doi{10.1007/978-3-540-27752-1}
#'
#' @examples
#' # Simulate VAR data
#' set.seed(123)
#' n <- 200
#' e <- matrix(rnorm(n * 3), n, 3)
#' y <- matrix(0, n, 3)
#' colnames(y) <- c("y1", "y2", "y3")
#' for (t in 3:n) {
#'   y[t, ] <- 0.3 * y[t-1, ] + 0.2 * y[t-2, ] + e[t, ]
#' }
#'
#' fit <- rbfmvar(y, lags = 2)
#' ir <- irf(fit, horizon = 20)
#' plot(ir)
#'
#' # With bootstrap confidence intervals
#' ir_boot <- irf(fit, horizon = 20, boot = 500, ci = 95)
#' plot(ir_boot)
#'
#' @export
irf <- function(object, horizon = 20, ortho = TRUE, boot = 0, ci = 90,
                seed = NULL) {
  if (!inherits(object, "rbfmvar")) {
    stop("'object' must be of class 'rbfmvar'.")
  }
  if (horizon < 1) {
    stop("'horizon' must be at least 1.")
  }
  if (ci <= 0 || ci >= 100) {
    stop("'ci' must be between 0 and 100.")
  }

  n <- object$n_vars
  varnames <- object$varnames

  irf_array <- .rbfmvar_irf_array(.rbfmvar_levels(object$F_plus, n, object$p_lags),
                                  object$Sigma_e, horizon, ortho)
  dimnames(irf_array) <- list(paste0("h=", 0:horizon), varnames, varnames)

  result <- list(
    irf = irf_array,
    irf_lower = NULL,
    irf_upper = NULL,
    horizon = horizon,
    varnames = varnames,
    ortho = ortho,
    boot = boot,
    ci = ci
  )

  if (boot > 0) {
    if (is.null(object$data)) {
      stop("'object' does not contain the data; re-estimate with rbfmvar() >= 2.1.0.")
    }
    if (!is.null(seed)) set.seed(seed)

    Y <- object$data
    TT <- nrow(Y)
    A <- .rbfmvar_levels(object$F_plus, n, object$p_lags)
    P <- length(A)
    u <- scale(object$residuals, center = TRUE, scale = FALSE)
    T_eff <- nrow(u)
    t0 <- TT - T_eff + 1

    irf_boot <- array(NA_real_, dim = c(boot, horizon + 1, n, n))
    for (b in seq_len(boot)) {
      # Recursive-design residual bootstrap: rebuild the sample from the
      # estimated levels VAR, the observed initial values and resampled
      # centred residuals, then re-estimate the RBFM-VAR
      ystar <- Y
      e_star <- u[sample.int(T_eff, T_eff, replace = TRUE), , drop = FALSE]
      for (t in t0:TT) {
        yt <- e_star[t - t0 + 1, ]
        for (k in seq_len(P)) yt <- yt + A[[k]] %*% ystar[t - k, ]
        ystar[t, ] <- yt
      }
      est <- tryCatch(rbfmvar_estimate(ystar, object$p_lags, object$kernel,
                                       object$bandwidth), error = function(e) NULL)
      if (is.null(est)) next
      irf_boot[b, , , ] <- .rbfmvar_irf_array(.rbfmvar_levels(est$F_plus, n, object$p_lags),
                                             est$Sigma_e, horizon, ortho)
    }

    alpha <- (100 - ci) / 200
    q_lo <- function(x) stats::quantile(x, alpha, na.rm = TRUE, names = FALSE)
    q_hi <- function(x) stats::quantile(x, 1 - alpha, na.rm = TRUE, names = FALSE)
    irf_lower <- apply(irf_boot, c(2, 3, 4), q_lo)
    irf_upper <- apply(irf_boot, c(2, 3, 4), q_hi)
    dimnames(irf_lower) <- dimnames(irf_upper) <- dimnames(irf_array)
    result$irf_lower <- irf_lower
    result$irf_upper <- irf_upper
    result$boot_failed <- sum(is.na(irf_boot[, 1, 1, 1]))
  }

  class(result) <- "rbfmvar_irf"
  result
}

#' Levels VAR Coefficients Implied by the RBFM-VAR Estimates
#'
#' Solves \eqn{(1-L)^2 y_t = \sum_j \Gamma_j L^j (1-L)^2 y_t + \Pi_1 L (1-L)
#' y_t + \Pi_2 L y_t + e_t} for \eqn{y_t = \sum_k A_k y_{t-k} + e_t}.
#'
#' @param F Coefficient matrix (regressors x n) in the order
#'   \eqn{\Gamma_1, \dots, \Gamma_{p-2}, \Pi_1, \Pi_2}.
#' @param n Number of variables.
#' @param p Lag order.
#' @return List of n x n matrices \eqn{A_1, \dots, A_P}, \eqn{P = \max(p, 2)}.
#' @keywords internal
.rbfmvar_levels <- function(F, n, p) {
  Ft <- t(F)
  ng <- max(p - 2, 0)
  P <- max(p, 2)
  # C(L) = sum_k C_k L^k with y-coefficients; A_k = -C_k for k >= 1
  C <- lapply(0:P, function(k) matrix(0, n, n))
  addL <- function(C, k, M) { C[[k + 1]] <- C[[k + 1]] + M; C }
  I <- diag(n)
  C <- addL(C, 0, I); C <- addL(C, 1, -2 * I); C <- addL(C, 2, I)
  if (ng > 0) for (j in 1:ng) {
    G <- Ft[, (j - 1) * n + 1:n, drop = FALSE]
    C <- addL(C, j, -G); C <- addL(C, j + 1, 2 * G); C <- addL(C, j + 2, -G)
  }
  Pi1 <- Ft[, ng * n + 1:n, drop = FALSE]
  Pi2 <- Ft[, ng * n + n + 1:n, drop = FALSE]
  C <- addL(C, 1, -Pi1); C <- addL(C, 2, Pi1)
  C <- addL(C, 1, -Pi2)
  lapply(1:P, function(k) -C[[k + 1]])
}

#' Moving Average Coefficients of a Levels VAR
#' @keywords internal
.rbfmvar_irf_array <- function(A, Sigma_e, horizon, ortho) {
  n <- nrow(Sigma_e)
  P <- length(A)
  S <- if (ortho) t(chol(Sigma_e)) else diag(n)
  Psi <- vector("list", horizon + 1)
  Psi[[1]] <- diag(n)
  for (h in seq_len(horizon)) {
    M <- matrix(0, n, n)
    for (k in seq_len(min(h, P))) M <- M + A[[k]] %*% Psi[[h - k + 1]]
    Psi[[h + 1]] <- M
  }
  out <- array(0, dim = c(horizon + 1, n, n))
  for (h in 0:horizon) out[h + 1, , ] <- Psi[[h + 1]] %*% S
  out
}

#' @export
print.rbfmvar_irf <- function(x, ...) {
  cat("\nImpulse Response Functions (RBFM-VAR)\n")
  cat("=====================================\n\n")
  cat("Horizon:", x$horizon, "periods\n")
  cat("Variables:", paste(x$varnames, collapse = ", "), "\n")
  cat("Orthogonalized:", x$ortho, "\n")

  if (x$boot > 0) {
    cat("Bootstrap:", x$boot, "replications,", x$ci, "% CI\n")
  }

  cat("\nUse plot() to visualize the IRFs.\n")
  invisible(x)
}

#' @export
plot.rbfmvar_irf <- function(x, shock = NULL, response = NULL, ...) {
  n <- length(x$varnames)
  horizon <- x$horizon

  # Determine which plots to show
  if (is.null(shock)) shock <- x$varnames
  if (is.null(response)) response <- x$varnames

  shock_idx <- which(x$varnames %in% shock)
  response_idx <- which(x$varnames %in% response)

  n_plots <- length(shock_idx) * length(response_idx)

  # Set up plotting grid
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  n_row <- ceiling(sqrt(n_plots))
  n_col <- ceiling(n_plots / n_row)
  graphics::par(mfrow = c(n_row, n_col), mar = c(4, 4, 2, 1))

  h <- 0:horizon

  for (j in shock_idx) {
    for (i in response_idx) {
      irf_vals <- x$irf[, i, j]

      ylim <- range(irf_vals)
      if (!is.null(x$irf_lower)) {
        ylim <- range(c(ylim, x$irf_lower[, i, j], x$irf_upper[, i, j]))
      }

      plot(h, irf_vals, type = "l", lwd = 2,
           xlab = "Horizon", ylab = "Response",
           main = paste(x$varnames[i], "<-", x$varnames[j]),
           ylim = ylim, ...)

      graphics::abline(h = 0, lty = 2, col = "gray")

      if (!is.null(x$irf_lower)) {
        graphics::polygon(c(h, rev(h)),
                          c(x$irf_lower[, i, j], rev(x$irf_upper[, i, j])),
                          col = grDevices::rgb(0, 0, 1, 0.2), border = NA)
      }
    }
  }
}
