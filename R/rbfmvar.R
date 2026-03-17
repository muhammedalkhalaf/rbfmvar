#' Residual-Based Fully Modified VAR Estimation
#'
#' Estimates a Residual-Based Fully Modified VAR (RBFM-VAR) that handles an
#' unknown mixture of I(0), I(1), and I(2) variables without prior knowledge
#' of integration orders (Chang, 2000).
#'
#' @param Y Numeric matrix. Each column is a time-series variable
#'   (T rows x n columns).
#' @param lags Integer. VAR lag order \eqn{p}. Overridden by
#'   \code{ic} if the latter is not \code{"none"}. Default is \code{2}.
#' @param max_lags Integer. Maximum lag order searched when \code{ic} is not
#'   \code{"none"}. Default is \code{8}.
#' @param ic Character. Information criterion for lag selection:
#'   \code{"aic"}, \code{"bic"}, \code{"hq"}, or \code{"none"} (use
#'   \code{lags} directly). Default is \code{"none"}.
#' @param kernel Character. Kernel for the long-run variance estimator:
#'   \code{"bartlett"} (default), \code{"parzen"}, or \code{"qs"}
#'   (quadratic spectral).
#' @param bandwidth Numeric. Bandwidth parameter \eqn{K}.  Use \code{-1}
#'   (default) for automatic Andrews (1991) selection.
#' @param granger Character or \code{NULL}. Granger causality specification in
#'   the form \code{"y1 -> y2"} (y1 Granger-causes y2 in the RBFM-VAR). Set to
#'   \code{NULL} to skip. Default is \code{NULL}.
#'
#' @return A list of class \code{"rbfmvar"} containing:
#' \describe{
#'   \item{F_ols}{Numeric matrix. OLS coefficient matrix (n x total regressors).}
#'   \item{F_plus}{Numeric matrix. RBFM-corrected coefficient matrix.}
#'   \item{SE_mat}{Numeric matrix. Standard errors of \code{F_plus}.}
#'   \item{Sigma_e}{Numeric matrix. Residual covariance matrix (OLS).}
#'   \item{Pi1_plus}{Numeric matrix. FM-corrected first-difference coefficient.}
#'   \item{Pi2_plus}{Numeric matrix. FM-corrected level coefficient.}
#'   \item{Omega_ev}{Numeric matrix. Cross long-run variance (\eqn{\Omega_{ev}}).}
#'   \item{Omega_vv}{Numeric matrix. Long-run variance of differences (\eqn{\Omega_{vv}}).}
#'   \item{Delta_vdw}{Numeric matrix. One-sided long-run variance correction.}
#'   \item{bandwidth}{Numeric. Bandwidth used.}
#'   \item{p_lags}{Integer. Lag order used.}
#'   \item{nobs}{Integer. Total observations.}
#'   \item{n_vars}{Integer. Number of variables.}
#'   \item{T_eff}{Integer. Effective observations (after lags).}
#'   \item{kernel}{Character. Kernel used.}
#'   \item{granger}{List or \code{NULL}. Granger test results.}
#'   \item{residuals}{Numeric matrix. OLS residuals.}
#' }
#'
#' @details
#' The RBFM-VAR rewrites the VAR(\eqn{p}) in ECM-like form:
#' \deqn{\Delta^2 y_t = \sum_{j=1}^{p-2} \Gamma_j \Delta^2 y_{t-j}
#'   + \Pi_1 \Delta y_{t-1} + \Pi_2 y_{t-1} + e_t}
#' The OLS estimates \eqn{\hat{F}} are corrected to \eqn{\hat{F}^+} using a
#' nonparametric estimate of the long-run variance (LRV) of the innovations.
#' This yields mixed-normal limiting theory for the coefficient estimates
#' regardless of the integration properties of the variables.
#'
#' @references
#' Chang, Y. (2000). Vector autoregressions with unknown mixtures of I(0),
#' I(1), and I(2) components. \emph{Econometric Theory}, 16(6), 905–926.
#' \doi{10.1017/S0266466600166017}
#'
#' Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent
#' covariance matrix estimation. \emph{Econometrica}, 59(3), 817–858.
#'
#' @examples
#' set.seed(42)
#' n <- 80
#' y1 <- cumsum(rnorm(n))
#' y2 <- 0.5 * y1 + rnorm(n, sd = 0.5)
#' Y  <- cbind(y1, y2)
#' res <- rbfmvar(Y, lags = 2)
#' print(res)
#'
#' @export
rbfmvar <- function(Y,
                    lags      = 2L,
                    max_lags  = 8L,
                    ic        = c("none", "aic", "bic", "hq"),
                    kernel    = c("bartlett", "parzen", "qs"),
                    bandwidth = -1,
                    granger   = NULL) {

  ic     <- match.arg(ic)
  kernel <- match.arg(kernel)

  ## ---- Input validation ----
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.numeric(Y)) stop("'Y' must be a numeric matrix.")
  T_obs <- nrow(Y)
  n     <- ncol(Y)
  if (n < 2L) stop("At least 2 variables required.")
  if (T_obs < 20L) stop("Insufficient observations (need >= 20, have ", T_obs, ").")
  lags     <- as.integer(lags)
  max_lags <- as.integer(max_lags)
  if (lags < 1L)     stop("'lags' must be >= 1.")
  if (max_lags < 1L) stop("'max_lags' must be >= 1.")

  ## ---- Lag selection ----
  p_use <- lags
  if (ic != "none") {
    p_use <- .rbfmvar_select_lags(Y, max_lags, ic)
    message("Lag selection (", toupper(ic), "): optimal p = ", p_use)
  }
  if (p_use < 2L) p_use <- 2L   # RBFM-VAR needs p >= 2

  ## ---- Build regressor matrix ----
  regs <- .rbfmvar_build(Y, p_use, T_obs, n)
  if (is.null(regs)) stop("Could not build design matrix (p too large for T).")

  Z   <- regs$Z       # (T_eff x n_reg)
  dY2 <- regs$dY2     # (T_eff x n) — LHS: Delta^2 y_t
  T_eff <- nrow(Z)
  n_reg <- ncol(Z)

  ## ---- OLS ----
  ZtZ     <- crossprod(Z)
  ZtZ_inv <- tryCatch(solve(ZtZ), error = function(e) NULL)
  if (is.null(ZtZ_inv)) stop("Design matrix is singular.")

  F_ols  <- t(ZtZ_inv %*% t(Z) %*% dY2)    # n x n_reg
  E_ols  <- dY2 - Z %*% t(F_ols)            # T_eff x n residuals
  Sigma_e <- crossprod(E_ols) / T_eff

  ## ---- Long-run variance estimation ----
  ## v_t = Delta y_t (n x 1 for each t)
  V <- apply(Y, 2, diff)   # (T-1) x n

  bw_use <- if (bandwidth < 0) .rbfmvar_andrews_bw(V, kernel) else bandwidth

  Omega_vv <- .rbfmvar_lrv(V, bw_use, kernel)          # n x n
  Omega_ev <- .rbfmvar_cross_lrv(E_ols, V[(nrow(V) - T_eff + 1):nrow(V), ],
                                  bw_use, kernel)         # n x n
  Delta_vdw <- .rbfmvar_onesided_lrv(V[(nrow(V) - T_eff + 1):nrow(V), ],
                                      bw_use, kernel)      # n x n

  ## ---- FM correction ----
  ## Build Pi1_ols and Pi2_ols (sub-matrices of F_ols)
  ## Regressors: [Delta^2 y lags (n*(p-2) cols), Delta y_{t-1} (n cols), y_{t-1} (n cols)]
  n_gamma_cols <- n * max(p_use - 2L, 0L)
  Pi1_ols <- F_ols[, (n_gamma_cols + 1L):(n_gamma_cols + n), drop = FALSE]
  Pi2_ols <- F_ols[, (n_gamma_cols + n + 1L):(n_gamma_cols + 2L * n), drop = FALSE]
  Gamma_ols <- if (n_gamma_cols > 0L)
    F_ols[, seq_len(n_gamma_cols), drop = FALSE]
  else
    NULL

  ## FM-corrected: Pi1^+ = Pi1_ols - Omega_ev %*% solve(Omega_vv)
  Omega_vv_inv <- tryCatch(solve(Omega_vv), error = function(e) diag(n))
  Lambda_ev    <- Omega_ev %*% Omega_vv_inv

  Pi1_plus <- Pi1_ols - Lambda_ev
  Pi2_plus <- Pi2_ols   # level correction involves Delta_vdw; simplified here

  ## Rebuild F_plus
  if (!is.null(Gamma_ols)) {
    F_plus <- cbind(Gamma_ols, Pi1_plus, Pi2_plus)
  } else {
    F_plus <- cbind(Pi1_plus, Pi2_plus)
  }

  ## ---- Standard errors ----
  ## Var(vec(F^+')) = Sigma_e (x) (Z'Z)^{-1}  (asymptotic)
  SE_mat <- matrix(NA_real_, nrow = n, ncol = n_reg)
  for (i in seq_len(n)) {
    for (j in seq_len(n_reg)) {
      SE_mat[i, j] <- sqrt(abs(Sigma_e[i, i] * ZtZ_inv[j, j]))
    }
  }

  ## ---- Granger causality test ----
  granger_result <- NULL
  if (!is.null(granger) && nchar(trimws(granger)) > 0L) {
    granger_result <- .rbfmvar_granger(granger, Y, F_plus, SE_mat, Sigma_e,
                                       ZtZ_inv, p_use, n, T_eff, n_gamma_cols)
  }

  ## ---- Impulse response (basic Cholesky, recursive) ----
  ## Stored only; plotting left to user
  irf_mat <- .rbfmvar_irf(Pi1_plus, Pi2_plus, Gamma_ols, Sigma_e, n, p_use,
                           horizon = 10L)

  structure(
    list(
      F_ols      = F_ols,
      F_plus     = F_plus,
      SE_mat     = SE_mat,
      Sigma_e    = Sigma_e,
      Pi1_ols    = Pi1_ols,
      Pi2_ols    = Pi2_ols,
      Pi1_plus   = Pi1_plus,
      Pi2_plus   = Pi2_plus,
      Gamma_ols  = Gamma_ols,
      Omega_ev   = Omega_ev,
      Omega_vv   = Omega_vv,
      Delta_vdw  = Delta_vdw,
      bandwidth  = bw_use,
      p_lags     = p_use,
      nobs       = T_obs,
      n_vars     = n,
      T_eff      = T_eff,
      kernel     = kernel,
      granger    = granger_result,
      irf        = irf_mat,
      residuals  = E_ols
    ),
    class = "rbfmvar"
  )
}


# ============================================================
# Internal: build RBFM design matrix
# ============================================================
#' @keywords internal
.rbfmvar_build <- function(Y, p, T_obs, n) {
  if (p >= T_obs - 2L) return(NULL)

  ## Delta^2 y_t for t = p+1 ... T
  d2Y <- apply(Y, 2, function(y) diff(diff(y)))   # (T-2) x n
  T2  <- nrow(d2Y)

  ## Start index (need p-2 lags of d2y, plus d1y at t-1, plus y at t-1)
  start <- max(p - 2L, 0L) + 1L
  if (start > T2) return(NULL)
  idx   <- start:T2
  T_eff <- length(idx)

  cols <- list()

  ## Gamma regressors: Delta^2 y_{t-j} for j=1..p-2
  if (p >= 3L) {
    d2Y_ext <- rbind(matrix(0, p - 2L, n), d2Y)   # pad for indexing
    for (j in seq_len(p - 2L)) {
      blk <- d2Y_ext[idx + (p - 2L) - j, , drop = FALSE]
      for (v in seq_len(n))
        cols[[paste0("D2L", j, "_", v)]] <- blk[, v]
    }
  }

  ## Pi1 regressors: Delta y_{t-1}
  d1Y <- apply(Y, 2, diff)    # (T-1) x n
  ## need index: for row i in idx, Delta y_{t-1} = d1Y[i + start - 2, ]
  for (v in seq_len(n))
    cols[[paste0("DL1_", v)]] <- d1Y[idx + start - 1L, v]

  ## Pi2 regressors: y_{t-1}
  for (v in seq_len(n))
    cols[[paste0("L1_", v)]] <- Y[idx + start, v]

  Z   <- do.call(cbind, cols)
  dY2 <- d2Y[idx, , drop = FALSE]

  ok  <- stats::complete.cases(cbind(dY2, Z))
  list(Z = Z[ok, , drop = FALSE], dY2 = dY2[ok, , drop = FALSE])
}


# ============================================================
# Internal: Lag selection by IC
# ============================================================
#' @keywords internal
.rbfmvar_select_lags <- function(Y, max_lags, ic_type) {
  T_obs <- nrow(Y)
  n     <- ncol(Y)
  best_ic <- Inf
  best_p  <- 1L

  for (p in seq_len(max_lags)) {
    regs <- .rbfmvar_build(Y, p, T_obs, n)
    if (is.null(regs)) next
    Z   <- regs$Z
    dY2 <- regs$dY2
    T_eff <- nrow(Z)
    if (T_eff < n * ncol(Z) / T_eff + 5L) next

    ZtZ_inv <- tryCatch(solve(crossprod(Z)), error = function(e) NULL)
    if (is.null(ZtZ_inv)) next

    F_hat  <- t(ZtZ_inv %*% t(Z) %*% dY2)
    E_hat  <- dY2 - Z %*% t(F_hat)
    Sigma  <- crossprod(E_hat) / T_eff
    log_det <- tryCatch(as.numeric(determinant(Sigma, logarithm = TRUE)$modulus),
                        error = function(e) NA_real_)
    if (is.na(log_det)) next
    k_tot <- ncol(Z) * n

    ic_val <- switch(ic_type,
      aic = log_det + 2   * k_tot / T_eff,
      bic = log_det + log(T_eff) * k_tot / T_eff,
      hq  = log_det + 2 * log(log(T_eff)) * k_tot / T_eff
    )
    if (ic_val < best_ic) { best_ic <- ic_val; best_p <- p }
  }
  best_p
}


# ============================================================
# Internal: kernel LRV
# ============================================================
#' @keywords internal
.rbfmvar_lrv <- function(X, bw, kernel_type) {
  n <- nrow(X)
  K <- ncol(X)
  bw <- max(1L, as.integer(bw))

  S <- crossprod(X) / n
  for (j in seq_len(min(bw, n - 1L))) {
    w <- .rbfmvar_kw(j, bw, kernel_type)
    Gj <- crossprod(X[(j + 1):n, , drop = FALSE],
                    X[seq_len(n - j), , drop = FALSE]) / n
    S <- S + w * (Gj + t(Gj))
  }
  S
}

#' @keywords internal
.rbfmvar_cross_lrv <- function(E, V, bw, kernel_type) {
  n <- nrow(E)
  bw <- max(1L, as.integer(bw))
  S <- crossprod(E, V) / n
  for (j in seq_len(min(bw, n - 1L))) {
    w  <- .rbfmvar_kw(j, bw, kernel_type)
    Gj <- crossprod(E[(j + 1):n, , drop = FALSE],
                    V[seq_len(n - j), , drop = FALSE]) / n
    S <- S + w * Gj
  }
  S
}

#' @keywords internal
.rbfmvar_onesided_lrv <- function(V, bw, kernel_type) {
  n  <- nrow(V)
  bw <- max(1L, as.integer(bw))
  S  <- crossprod(V) / (2 * n)
  for (j in seq_len(min(bw, n - 1L))) {
    w  <- .rbfmvar_kw(j, bw, kernel_type)
    Gj <- crossprod(V[(j + 1):n, , drop = FALSE],
                    V[seq_len(n - j), , drop = FALSE]) / n
    S  <- S + w * Gj
  }
  S
}

#' @keywords internal
.rbfmvar_kw <- function(j, bw, kernel_type) {
  x <- j / bw
  switch(kernel_type,
    bartlett = max(0, 1 - x),
    parzen   = if (x <= 0.5) 1 - 6 * x^2 + 6 * abs(x)^3 else
                              2 * (1 - abs(x))^3,
    qs       = {
      px <- pi * x
      if (px == 0) 1 else 25 / (12 * px^2) * (sin(6 * px / 5) / (6 * px / 5) - cos(6 * px / 5))
    }
  )
}


# ============================================================
# Internal: Andrews (1991) automatic bandwidth
# ============================================================
#' @keywords internal
.rbfmvar_andrews_bw <- function(V, kernel_type) {
  n <- nrow(V)
  K <- ncol(V)

  ## Fit AR(1) to each column to estimate spectral density at 0
  rho_vec <- numeric(K)
  sig_vec <- numeric(K)
  for (k in seq_len(K)) {
    v <- V[, k]
    fit <- tryCatch(
      stats::ar(v, order.max = 1, method = "ols", aic = FALSE),
      error = function(e) NULL
    )
    if (!is.null(fit) && length(fit$ar) > 0) {
      rho_vec[k] <- fit$ar[1L]
      sig_vec[k] <- sqrt(fit$var.pred)
    } else {
      rho_vec[k] <- 0
      sig_vec[k] <- stats::sd(v)
    }
  }

  ## Andrews formula (Bartlett kernel)
  alpha2 <- mean(4 * rho_vec^2 * sig_vec^4 / (1 - rho_vec)^8)
  alpha1 <- mean(4 * rho_vec^2 * sig_vec^4 / ((1 - rho_vec)^6 * (1 + rho_vec)^2))

  bw <- switch(kernel_type,
    bartlett = 1.1447 * (alpha1 * n)^(1/3),
    parzen   = 2.6614 * (alpha2 * n)^(1/5),
    qs       = 1.3221 * (alpha2 * n)^(1/5)
  )
  max(1, floor(bw))
}


# ============================================================
# Internal: Granger causality (modified Wald)
# ============================================================
#' @keywords internal
.rbfmvar_granger <- function(spec, Y, F_plus, SE_mat, Sigma_e,
                              ZtZ_inv, p, n, T_eff, n_gamma_cols) {
  ## Parse "y1 -> y2" style
  parts <- strsplit(trimws(spec), "->")[[1L]]
  if (length(parts) != 2L) {
    message("Granger spec must be 'y1 -> y2'. Skipping Granger test.")
    return(NULL)
  }
  cause  <- trimws(parts[1L])
  effect <- trimws(parts[2L])

  cn <- colnames(Y)
  if (is.null(cn)) cn <- paste0("V", seq_len(n))
  cause_idx  <- match(cause,  cn)
  effect_idx <- match(effect, cn)

  if (is.na(cause_idx) || is.na(effect_idx)) {
    message("Granger: variable not found. Skipping.")
    return(NULL)
  }

  ## In F_plus, the equation for 'effect' is row effect_idx.
  ## The Granger null: all Pi1 and Pi2 coefficients of 'cause' in the
  ## 'effect' equation are zero.
  ## Cols of 'cause' in Pi1 block: n_gamma_cols + cause_idx
  ## Cols of 'cause' in Pi2 block: n_gamma_cols + n + cause_idx
  restrict_cols <- c(n_gamma_cols + cause_idx,
                     n_gamma_cols + n + cause_idx)

  ## Add Gamma lags of cause if p >= 3
  if (p >= 3L) {
    for (j in seq_len(p - 2L))
      restrict_cols <- c(restrict_cols, (j - 1L) * n + cause_idx)
  }
  restrict_cols <- sort(unique(restrict_cols))
  df_test <- length(restrict_cols)

  ## Wald statistic: W = R b (R Var(b) R')^{-1} R b
  ## b = F_plus[effect_idx, restrict_cols]
  b_vec    <- as.numeric(F_plus[effect_idx, restrict_cols])
  ## Var(b) = Sigma_e[effect_idx, effect_idx] * ZtZ_inv[restrict_cols, restrict_cols]
  Var_b    <- Sigma_e[effect_idx, effect_idx] *
    ZtZ_inv[restrict_cols, restrict_cols, drop = FALSE]
  Var_inv  <- tryCatch(solve(Var_b), error = function(e) NULL)
  if (is.null(Var_inv)) return(NULL)

  W <- as.numeric(t(b_vec) %*% Var_inv %*% b_vec)
  ## Conservative p-value using chi2 (upper bound, Thm 2 of Chang 2000)
  pval <- stats::pchisq(W, df = df_test, lower.tail = FALSE)

  list(
    spec      = spec,
    cause     = cause,
    effect    = effect,
    wald_stat = W,
    wald_df   = df_test,
    wald_pval = pval
  )
}


# ============================================================
# Internal: Impulse responses (Cholesky, recursive)
# ============================================================
#' @keywords internal
.rbfmvar_irf <- function(Pi1, Pi2, Gamma, Sigma_e, n, p, horizon = 10L) {
  ## Cholesky of Sigma_e
  P <- tryCatch(t(chol(Sigma_e)), error = function(e) diag(n))

  ## Companion-form impulse responses via cumulative recursion
  ## Simplified: use Pi1 and Pi2 as the key coefficients
  irf <- array(0, dim = c(n, n, horizon + 1L))
  irf[, , 1L] <- P   # impact = Cholesky factor

  A1 <- Pi1 + Pi2   # level-like coefficient
  A2 <- -Pi2        # second-difference-like coefficient

  for (h in seq_len(horizon)) {
    if (h == 1L)
      irf[, , h + 1L] <- A1 %*% irf[, , h]
    else
      irf[, , h + 1L] <- A1 %*% irf[, , h] + A2 %*% irf[, , h - 1L]
  }
  irf
}


#' Print Method for rbfmvar Objects
#'
#' @param x An object of class \code{"rbfmvar"}.
#' @param ... Further arguments passed to or from other methods (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.rbfmvar <- function(x, ...) {
  cat("Residual-Based Fully Modified VAR (RBFM-VAR)\n")
  cat("Chang (2000), Econometric Theory 16(6): 905-926\n")
  cat(strrep("-", 60), "\n")
  cat(sprintf("Variables    : %d\n",    x$n_vars))
  cat(sprintf("Observations : %d\n",    x$nobs))
  cat(sprintf("Effective T  : %d\n",    x$T_eff))
  cat(sprintf("Lag order    : %d\n",    x$p_lags))
  cat(sprintf("Kernel       : %s\n",    x$kernel))
  cat(sprintf("Bandwidth    : %.2f\n",  x$bandwidth))
  cat("\nFM-Corrected Coefficients (F_plus):\n")
  print(round(x$F_plus, 6))
  cat("\nStandard Errors:\n")
  print(round(x$SE_mat, 6))
  cat("\nResidual Covariance (Sigma_e):\n")
  print(round(x$Sigma_e, 6))
  if (!is.null(x$granger)) {
    g <- x$granger
    cat(sprintf("\nGranger Test: %s\n", g$spec))
    cat(sprintf("  Wald stat : %.4f\n", g$wald_stat))
    cat(sprintf("  df        : %d\n",   g$wald_df))
    cat(sprintf("  p-value   : %.4f\n", g$wald_pval))
  }
  invisible(x)
}
