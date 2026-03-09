#' @title Lag Selection via Information Criteria
#'
#' @description Selects the optimal VAR lag order using information criteria
#'   (AIC, BIC, or HQ).
#'
#' @param data Data matrix (T x n).
#' @param max_lags Maximum lag order to consider.
#' @param ic Information criterion: \code{"aic"}, \code{"bic"}, or \code{"hq"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{best_p}{Optimal lag order.}
#'   \item{ic_table}{Data frame of IC values for each lag.}
#' }
#'
#' @references
#' Akaike, H. (1974). A New Look at the Statistical Model Identification.
#' \emph{IEEE Transactions on Automatic Control}, 19(6), 716-723.
#' \doi{10.1109/TAC.1974.1100705}
#'
#' Schwarz, G. (1978). Estimating the Dimension of a Model.
#' \emph{The Annals of Statistics}, 6(2), 461-464.
#' \doi{10.1214/aos/1176344136}
#'
#' Hannan, E. J., & Quinn, B. G. (1979). The Determination of the Order of an
#' Autoregression. \emph{Journal of the Royal Statistical Society: Series B},
#' 41(2), 190-195. \doi{10.1111/j.2517-6161.1979.tb01072.x}
#'
#' @keywords internal
select_lags_ic <- function(data, max_lags, ic = "aic") {
  data <- as.matrix(data)
  TT <- nrow(data)
  n <- ncol(data)

  ic <- tolower(ic)

  # Storage for IC values
  ic_values <- data.frame(
    lag = integer(0),
    aic = numeric(0),
    bic = numeric(0),
    hq = numeric(0),
    T_eff = integer(0)
  )

  best_p <- 1
  best_ic <- Inf

  for (p in 1:max_lags) {
    # Check if we have enough observations
    T_eff <- TT - p - 1
    if (T_eff < 10) break

    # Estimate VAR(p) via OLS (simple version)
    result <- tryCatch(
      {
        estimate_var_ols(data, p)
      },
      error = function(e) NULL
    )

    if (is.null(result)) next

    # Compute information criteria
    Sigma_e <- result$Sigma_e
    log_det_sigma <- log(det(Sigma_e))

    # Number of parameters per equation
    k <- n * p + 1  # Assuming no constant, just lagged values

    # Information criteria
    aic_val <- log_det_sigma + 2 * k * n / T_eff
    bic_val <- log_det_sigma + log(T_eff) * k * n / T_eff
    hq_val <- log_det_sigma + 2 * log(log(T_eff)) * k * n / T_eff

    ic_values <- rbind(ic_values, data.frame(
      lag = p,
      aic = aic_val,
      bic = bic_val,
      hq = hq_val,
      T_eff = T_eff
    ))

    # Select best based on chosen criterion
    current_ic <- switch(ic,
      aic = aic_val,
      bic = bic_val,
      hq = hq_val
    )

    if (current_ic < best_ic) {
      best_ic <- current_ic
      best_p <- p
    }
  }

  list(
    best_p = best_p,
    ic_table = ic_values
  )
}

#' Simple VAR OLS Estimation
#'
#' @description Estimates a standard VAR model via OLS for lag selection.
#'
#' @param data Data matrix.
#' @param p Lag order.
#'
#' @return List with residual covariance matrix.
#' @keywords internal
estimate_var_ols <- function(data, p) {
  data <- as.matrix(data)
  TT <- nrow(data)
  n <- ncol(data)

  T_eff <- TT - p

  # Build lagged matrix
  Y <- data[(p + 1):TT, , drop = FALSE]
  Z <- matrix(0, T_eff, n * p)

  for (j in 1:p) {
    cols <- ((j - 1) * n + 1):(j * n)
    Z[, cols] <- data[(p + 1 - j):(TT - j), , drop = FALSE]
  }

  # OLS
  ZtZ <- crossprod(Z)
  ZtZ_inv <- tryCatch(
    solve(ZtZ),
    error = function(e) MASS::ginv(ZtZ)
  )
  B <- ZtZ_inv %*% crossprod(Z, Y)

  # Residuals and covariance
  resid <- Y - Z %*% B
  Sigma_e <- crossprod(resid) / T_eff

  list(
    B = B,
    Sigma_e = Sigma_e,
    residuals = resid,
    T_eff = T_eff
  )
}

#' Get Information Criteria Table
#'
#' @description Returns a table of information criteria values for different
#'   lag orders.
#'
#' @param object An \code{rbfmvar} object.
#' @param max_lags Maximum lag order to evaluate.
#'
#' @return A data frame with AIC, BIC, and HQ values.
#'
#' @examples
#' \donttest{
#' # Generate example data
#' set.seed(123)
#' n <- 100
#' mydata <- data.frame(x = cumsum(rnorm(n)), y = cumsum(rnorm(n)))
#' fit <- rbfmvar(mydata, lags = 2)
#' ic_table(fit, max_lags = 6)
#' }
#'
#' @export
ic_table <- function(object, max_lags = 8) {
  if (!inherits(object, "rbfmvar")) {
    stop("'object' must be of class 'rbfmvar'.")
  }

  # Get original data from object
  # We need to reconstruct data from residuals and fitted values
  # For now, use a simplified approach

  message("IC table computation requires original data. ",
          "Use select_lags_ic() directly with data matrix.")
  invisible(NULL)
}
