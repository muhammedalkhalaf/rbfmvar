#' @title Residual-Based Fully Modified VAR Estimation
#'
#' @description Estimates a Residual-Based Fully Modified Vector Autoregression
#'   (RBFM-VAR) model following Chang (2000). The RBFM-VAR procedure extends
#'   Phillips (1995) FM-VAR to handle any unknown mixture of I(0), I(1), and
#'   I(2) components without prior knowledge of the number or location of
#'   unit roots.
#'
#' @param data A numeric matrix or data frame containing the time series
#'   variables. Must have at least 2 columns.
#' @param lags Integer. The VAR lag order p. Must be at least 1. Default is 2.
#' @param max_lags Integer. Maximum number of lags to consider for information
#'   criterion selection. Default is 8.
#' @param ic Character string specifying the information criterion for lag
#'   selection: \code{"aic"}, \code{"bic"}, \code{"hq"}, or \code{"none"}
#'   (use \code{lags} directly). Default is \code{"none"}.
#' @param kernel Character string specifying the kernel for long-run variance
#'   estimation: \code{"bartlett"}, \code{"parzen"}, or \code{"qs"} (Quadratic
#'   Spectral). Default is \code{"bartlett"}.
#' @param bandwidth Numeric. Bandwidth for kernel estimation. If \code{-1}
#'   (default), automatic bandwidth selection via Andrews (1991) is used.
#' @param level Numeric. Confidence level for coefficient intervals (0-100).
#'   Default is 95.
#'
#' @details
#' The RBFM-VAR model is specified as:
#' \deqn{\Delta^2 y_t = \sum_{j=1}^{p-2} \Gamma_j \Delta^2 y_{t-j} + \Pi_1 \Delta y_{t-1} + \Pi_2 y_{t-1} + e_t}
#'
#' where \eqn{\Delta} is the difference operator and \eqn{\Delta^2 = \Delta \circ \Delta}.
#'
#' The coefficients are estimated by the RBFM-VAR procedure of Chang (2000,
#' equations 10 to 13); see the internal function \code{rbfmvar_estimate}
#' for the formulas. The automatic bandwidth applies Andrews (1991) to the
#' series \eqn{\hat v_t}; the correction depends on the bandwidth in finite
#' samples, and Chang (2000) does not state the kernel or bandwidth used in
#' its simulations.
#'
#' The FM+ correction eliminates the second-order asymptotic bias that arises

#' from the correlation between the regression errors and the innovations in
#' integrated regressors. The estimator achieves:
#' \itemize{
#'   \item Zero mean mixed normal limiting distribution
#'   \item Chi-square Wald statistics for hypothesis testing
#'   \item Robustness to unknown integration orders
#' }
#'
#' @return An object of class \code{"rbfmvar"} containing:
#' \describe{
#'   \item{F_ols}{OLS coefficient matrix.}
#'   \item{F_plus}{FM+ corrected coefficient matrix.}
#'   \item{SE_mat}{Standard errors for FM+ coefficients.}
#'   \item{Pi1_ols, Pi1_plus}{Coefficient matrices for \eqn{\Delta y_{t-1}}.}
#'   \item{Pi2_ols, Pi2_plus}{Coefficient matrices for \eqn{y_{t-1}}.}
#'   \item{Gamma_ols, Gamma_plus}{Coefficient matrices for \eqn{\Delta^2 y_{t-j}} (if p >= 3).}
#'   \item{Sigma_e}{Residual covariance matrix.}
#'   \item{Omega_ev, Omega_vv}{Long-run variance components.}
#'   \item{Delta_vdw}{One-sided long-run covariance for FM correction.}
#'   \item{N_hat}{OLS coefficient of \eqn{\Delta y_{t-1}} on \eqn{\Delta y_{t-2}} used in \eqn{\hat v_t}.}
#'   \item{XtX_inv}{\eqn{(X'X)^{-1}} of the regressors, used by \code{granger_test}.}
#'   \item{data}{The data matrix, used by \code{irf}, \code{forecast} and \code{ic_table}.}
#'   \item{residuals}{Matrix of residuals from FM+ estimation.}
#'   \item{fitted}{Matrix of fitted values.}
#'   \item{nobs}{Number of observations in original data.}
#'   \item{T_eff}{Effective sample size after differencing.}
#'   \item{n_vars}{Number of variables.}
#'   \item{p_lags}{VAR lag order used.}
#'   \item{bandwidth}{Bandwidth used for LRV estimation.}
#'   \item{kernel}{Kernel used for LRV estimation.}
#'   \item{ic}{Information criterion used (if any).}
#'   \item{varnames}{Variable names.}
#'   \item{call}{The matched call.}
#' }
#'
#' @references
#' Chang, Y. (2000). Vector Autoregressions with Unknown Mixtures of I(0), I(1),
#' and I(2) Components. \emph{Econometric Theory}, 16(6), 905-926.
#' \doi{10.1017/S0266466600166058}
#'
#' Phillips, P. C. B. (1995). Fully Modified Least Squares and Vector
#' Autoregression. \emph{Econometrica}, 63(5), 1023-1078.
#' \doi{10.2307/2171721}
#'
#' Andrews, D. W. K. (1991). Heteroskedasticity and Autocorrelation Consistent
#' Covariance Matrix Estimation. \emph{Econometrica}, 59(3), 817-858.
#' \doi{10.2307/2938229}
#'
#' @examples
#' # Simulate a simple VAR(2) process
#' set.seed(123)
#' n <- 200
#' e <- matrix(rnorm(n * 3), n, 3)
#' y <- matrix(0, n, 3)
#' for (t in 3:n) {
#'   y[t, ] <- 0.3 * y[t-1, ] + 0.2 * y[t-2, ] + e[t, ]
#' }
#' colnames(y) <- c("y1", "y2", "y3")
#'
#' # Estimate RBFM-VAR
#' fit <- rbfmvar(y, lags = 2)
#' summary(fit)
#'
#' # With automatic lag selection
#' fit_aic <- rbfmvar(y, max_lags = 6, ic = "aic")
#' summary(fit_aic)
#'
#' @export
rbfmvar <- function(data, lags = 2, max_lags = 8, ic = "none",
                    kernel = "bartlett", bandwidth = -1, level = 95) {

  # Capture call

  cl <- match.call()

  # Convert to matrix
  data <- as.matrix(data)
  if (!is.numeric(data)) {
    stop("'data' must be numeric.")
  }

  nobs <- nrow(data)
  n_vars <- ncol(data)

  # Get variable names
  varnames <- colnames(data)
  if (is.null(varnames)) {
    varnames <- paste0("V", seq_len(n_vars))
    colnames(data) <- varnames
  }

  # Validate inputs
  if (n_vars < 2) {
    stop("At least 2 variables required for VAR.")
  }

  if (lags < 1) {
    stop("'lags' must be at least 1.")
  }

  if (max_lags < 1) {
    stop("'max_lags' must be at least 1.")
  }

  kernel <- tolower(kernel)
  if (!kernel %in% c("bartlett", "parzen", "qs")) {
    stop("'kernel' must be 'bartlett', 'parzen', or 'qs'.")
  }

  ic <- tolower(ic)
  if (!ic %in% c("aic", "bic", "hq", "none")) {
    stop("'ic' must be 'aic', 'bic', 'hq', or 'none'.")
  }

  if (level <= 0 || level >= 100) {
    stop("'level' must be between 0 and 100.")
  }

  if (nobs < 20) {
    stop("Insufficient observations (need at least 20, have ", nobs, ").")
  }

  # Lag selection via IC
  p_use <- lags
  if (ic != "none") {
    ic_result <- select_lags_ic(data, max_lags, ic)
    p_use <- ic_result$best_p
    message("Lag selection (", toupper(ic), "): optimal p = ", p_use)
  }

  # Check sufficient observations for chosen lag
  min_obs <- p_use + 3  # Need at least p+3 for second differences
  if (nobs < min_obs) {
    stop("Insufficient observations for p = ", p_use,
         " (need at least ", min_obs, ", have ", nobs, ").")
  }

  # Run RBFM-VAR estimation
  result <- rbfmvar_estimate(data, p_use, kernel, bandwidth)

  # Add metadata
  result$nobs <- nobs
  result$n_vars <- n_vars
  result$p_lags <- p_use
  result$kernel <- kernel
  result$ic <- ic
  result$level <- level
  result$varnames <- varnames
  result$call <- cl
  result$data <- data

  class(result) <- "rbfmvar"
  result
}

#' Core RBFM-VAR Estimation
#'
#' Implements the RBFM-VAR estimator of Chang (2000), equations (10) to
#' (13). The regression is written in the second-difference form of
#' equation (2); since \eqn{y_t = \Delta^2 y_t + \Delta y_{t-1} + y_{t-1}},
#' the dependent variables of (2) and (3) differ by a linear function of
#' \eqn{w_t}, and the corrected estimates of the two forms differ only by
#' the identity added to \eqn{\Pi_1} and \eqn{\Pi_2}.
#'
#' With \eqn{z_t = (\Delta^2 y_{t-1}', \dots, \Delta^2 y_{t-p+2}')'} and
#' \eqn{w_t = (\Delta y_{t-1}', y_{t-1}')'},
#' \deqn{\hat F^+ = (Y'Z, Y^{+\prime}W + T\hat\Delta^+)(X'X)^{-1},}
#' \eqn{Y^{+\prime} = Y' - \hat\Omega_{\hat\varepsilon\hat v}
#' \hat\Omega_{\hat v\hat v}^{-1}\hat V'} and \eqn{\hat\Delta^+ =
#' \hat\Omega_{\hat\varepsilon\hat v}\hat\Omega_{\hat v\hat v}^{-1}
#' \hat\Delta_{\hat v\Delta w}}, where \eqn{\hat\varepsilon_t} are the OLS
#' residuals, \eqn{\hat v_t = (\Delta^2 y_{t-1}', \Delta y_{t-1}' - \hat
#' N\Delta y_{t-2}')'} (equation 11) with \eqn{\hat N} the OLS coefficient of
#' \eqn{\Delta y_{t-1}} on \eqn{\Delta y_{t-2}}, \eqn{\Delta w_t =
#' (\Delta^2 y_{t-1}', \Delta y_{t-1}')'}, \eqn{\hat\Omega} are two-sided
#' kernel long-run covariance estimates and \eqn{\hat\Delta_{\hat v\Delta
#' w} = \sum_{j \ge 0} k(j/K) T^{-1}\sum_t \hat v_{t+j}\Delta w_t'} is the
#' one-sided long-run covariance (Chang, 2000, p. 909:
#' \eqn{\Delta = \sum_{j\ge0} E u_j u_0'}).
#'
#' @param Y Data matrix (T x n).
#' @param p Lag order.
#' @param kernel Kernel type.
#' @param bandwidth Bandwidth (-1 for automatic).
#'
#' @return List of estimation results.
#' @keywords internal
rbfmvar_estimate <- function(Y, p, kernel, bandwidth) {
  Y <- as.matrix(Y)
  TT <- nrow(Y)
  n <- ncol(Y)

  # Series indexed by original time t = 1..TT (NA where undefined)
  dY <- rbind(NA, diff(Y))
  d2Y <- rbind(NA, diff(dY))
  lagm <- function(M, j) {
    if (j == 0) return(M)
    rbind(matrix(NA_real_, j, ncol(M)), M[seq_len(nrow(M) - j), , drop = FALSE])
  }

  # Estimation sample: the VAR(p) in levels needs y_{t-p}; v_t needs
  # Delta y_{t-2}, hence t >= 4
  t0 <- max(p + 1, 4)
  rows <- t0:TT
  T_eff <- length(rows)
  if (T_eff < 10) {
    stop("Effective sample size too small (T_eff = ", T_eff, ").")
  }

  y_dep <- d2Y[rows, , drop = FALSE]
  Zs <- if (p >= 3) {
    do.call(cbind, lapply(1:(p - 2), function(j) lagm(d2Y, j)[rows, , drop = FALSE]))
  } else matrix(0, T_eff, 0)
  W <- cbind(lagm(dY, 1)[rows, , drop = FALSE], lagm(Y, 1)[rows, , drop = FALSE])
  Z <- cbind(Zs, W)   # full regressor matrix X = (Z, W) of Chang (2000)

  n_gamma <- ncol(Zs)
  n_regs <- ncol(Z)

  # =========================================================================
  # OLS-VAR
  # =========================================================================
  ZtZ <- crossprod(Z)
  ZtZ_inv <- tryCatch(solve(ZtZ), error = function(e) MASS::ginv(ZtZ))
  F_ols <- ZtZ_inv %*% crossprod(Z, y_dep)          # (n_regs x n)
  e_ols <- y_dep - Z %*% F_ols
  Sigma_e <- crossprod(e_ols) / T_eff

  # =========================================================================
  # v_hat (equation 11) and Delta w
  # =========================================================================
  dY1 <- lagm(dY, 1)[rows, , drop = FALSE]
  dY2 <- lagm(dY, 2)[rows, , drop = FALSE]
  N_hat <- t(solve(crossprod(dY2), crossprod(dY2, dY1)))  # dY1 = dY2 N' + u
  v <- cbind(lagm(d2Y, 1)[rows, , drop = FALSE], dY1 - dY2 %*% t(N_hat))
  dW <- cbind(lagm(d2Y, 1)[rows, , drop = FALSE], dY1)

  if (bandwidth < 0) {
    bandwidth <- select_bandwidth_andrews(v, kernel)
  }

  Omega_ev <- estimate_cross_lrv(e_ols, v, kernel, bandwidth)
  Omega_vv <- estimate_cross_lrv(v, v, kernel, bandwidth)
  Omega_vv <- (Omega_vv + t(Omega_vv)) / 2
  Delta_vdw <- estimate_onesided_lrv(v, dW, kernel, bandwidth)

  # The limit of Omega_vv is singular in the stationary direction (Chang,
  # 2000, Remark (e)); a generalised inverse is used if it is singular
  Omega_vv_inv <- tryCatch(solve(Omega_vv), error = function(e) MASS::ginv(Omega_vv))
  B <- Omega_ev %*% Omega_vv_inv                     # (n x 2n)

  # =========================================================================
  # RBFM-VAR correction (equations 12 and 13)
  # =========================================================================
  Y_plus <- y_dep - v %*% t(B)                       # rows: y_t^+'
  Delta_plus <- B %*% Delta_vdw                      # (n x 2n)
  YW <- t(Y_plus) %*% W + T_eff * Delta_plus
  YX <- if (n_gamma > 0) cbind(t(y_dep) %*% Zs, YW) else YW
  F_plus <- t(YX %*% ZtZ_inv)                        # (n_regs x n)

  # =========================================================================
  # Standard errors: Sigma_e (x) (X'X)^{-1} (equation 20)
  # =========================================================================
  SE_mat <- sqrt(outer(diag(Sigma_e), diag(ZtZ_inv)))

  F_ols_t <- t(F_ols)
  F_plus_t <- t(F_plus)

  Gamma_ols <- NULL
  Gamma_plus <- NULL
  if (p >= 3) {
    Gamma_ols <- F_ols_t[, 1:n_gamma, drop = FALSE]
    Gamma_plus <- F_plus_t[, 1:n_gamma, drop = FALSE]
  }

  Pi1_start <- n_gamma + 1
  Pi1_end <- n_gamma + n
  Pi1_ols <- matrix(F_ols_t[, Pi1_start:Pi1_end], n, n)
  Pi1_plus <- matrix(F_plus_t[, Pi1_start:Pi1_end], n, n)

  Pi2_start <- n_gamma + n + 1
  Pi2_end <- n_regs
  Pi2_ols <- matrix(F_ols_t[, Pi2_start:Pi2_end], n, n)
  Pi2_plus <- matrix(F_plus_t[, Pi2_start:Pi2_end], n, n)

  rownames(Pi1_ols) <- colnames(Pi1_ols) <- colnames(Y)
  rownames(Pi1_plus) <- colnames(Pi1_plus) <- colnames(Y)
  rownames(Pi2_ols) <- colnames(Pi2_ols) <- colnames(Y)
  rownames(Pi2_plus) <- colnames(Pi2_plus) <- colnames(Y)
  rownames(Sigma_e) <- colnames(Sigma_e) <- colnames(Y)

  fitted <- Z %*% F_plus
  residuals <- y_dep - fitted
  colnames(residuals) <- colnames(Y)
  colnames(fitted) <- colnames(Y)

  regnames <- character(0)
  if (p >= 3) {
    for (j in 1:(p - 2)) regnames <- c(regnames, paste0("L", j, "D2.", colnames(Y)))
  }
  regnames <- c(regnames, paste0("LD.", colnames(Y)), paste0("L.", colnames(Y)))

  colnames(F_ols) <- colnames(Y)
  rownames(F_ols) <- regnames
  colnames(F_plus) <- colnames(Y)
  rownames(F_plus) <- regnames
  colnames(SE_mat) <- regnames
  rownames(SE_mat) <- colnames(Y)
  dimnames(ZtZ_inv) <- list(regnames, regnames)

  list(
    F_ols = F_ols,
    F_plus = F_plus,
    SE_mat = SE_mat,
    XtX_inv = ZtZ_inv,
    Pi1_ols = Pi1_ols,
    Pi1_plus = Pi1_plus,
    Pi2_ols = Pi2_ols,
    Pi2_plus = Pi2_plus,
    Gamma_ols = Gamma_ols,
    Gamma_plus = Gamma_plus,
    Sigma_e = Sigma_e,
    Omega_ev = Omega_ev,
    Omega_vv = Omega_vv,
    Delta_vdw = Delta_vdw,
    N_hat = N_hat,
    residuals = residuals,
    fitted = fitted,
    T_eff = T_eff,
    bandwidth = bandwidth,
    regnames = regnames
  )
}
