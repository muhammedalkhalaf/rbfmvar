#' @title Out-of-Sample Forecasting
#'
#' @description
#' The generic function \code{forecast} computes forecasts from time series models.
#'
#' @param object A model object.
#' @param ... Additional arguments passed to methods.
#'
#' @return Depends on the method dispatched. See \code{\link{forecast.rbfmvar}}
#'   for the RBFM-VAR method, which returns an object of class
#'   \code{"rbfmvar_forecast"}.
#'
#' @export
forecast <- function(object, ...) {

  UseMethod("forecast")
}

#' @title Out-of-Sample Forecasting for RBFM-VAR
#'
#' @description Generates out-of-sample forecasts from an RBFM-VAR model.
#'
#' @param object An \code{rbfmvar} object from \code{\link{rbfmvar}}.
#' @param h Integer. Forecast horizon (number of periods ahead). Default is 10.
#' @param level Numeric. Confidence level for prediction intervals (0-100).
#'   Default is 95.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' Forecasts of the levels \eqn{y_{T+s}} are generated iteratively from the
#' levels VAR implied by the RBFM-VAR coefficients, starting from the last
#' observations of the data. The forecast error variance at horizon s is
#' \eqn{\sum_{i=0}^{s-1}\Psi_i\Sigma_e\Psi_i'} with \eqn{\Psi_i} the moving
#' average coefficients; prediction intervals assume normal innovations and
#' ignore parameter estimation uncertainty.
#'
#' @return An object of class \code{"rbfmvar_forecast"} containing:
#' \describe{
#'   \item{mean}{Matrix of point forecasts (n x h).}
#'   \item{se}{Matrix of forecast standard errors (n x h).}
#'   \item{lower}{Matrix of lower prediction bounds (n x h).}
#'   \item{upper}{Matrix of upper prediction bounds (n x h).}
#'   \item{horizon}{Forecast horizon.}
#'   \item{level}{Confidence level.}
#'   \item{varnames}{Variable names.}
#' }
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
#' fc <- forecast(fit, h = 10)
#' print(fc)
#' plot(fc)
#'
#' @rdname forecast.rbfmvar
#' @export
forecast.rbfmvar <- function(object, h = 10, level = 95, ...) {
  if (!inherits(object, "rbfmvar")) {
    stop("'object' must be of class 'rbfmvar'.")
  }
  if (h < 1) {
    stop("'h' must be at least 1.")
  }
  if (level <= 0 || level >= 100) {
    stop("'level' must be between 0 and 100.")
  }
  if (is.null(object$data)) {
    stop("'object' does not contain the data; re-estimate with rbfmvar() >= 2.1.0.")
  }

  n <- object$n_vars
  varnames <- object$varnames
  Sigma_e <- object$Sigma_e
  A <- .rbfmvar_levels(object$F_plus, n, object$p_lags)
  P <- length(A)
  Y <- object$data
  TT <- nrow(Y)

  # Point forecasts of the levels from the implied levels VAR
  path <- rbind(Y[(TT - P + 1):TT, , drop = FALSE], matrix(NA_real_, h, n))
  for (s in 1:h) {
    yt <- rep(0, n)
    for (k in 1:P) yt <- yt + A[[k]] %*% path[P + s - k, ]
    path[P + s, ] <- yt
  }
  fc_mean <- t(path[P + 1:h, , drop = FALSE])

  # Forecast error variance: sum_{i < s} Psi_i Sigma_e Psi_i'
  Psi <- .rbfmvar_irf_array(A, Sigma_e, h - 1, ortho = FALSE)
  fc_se <- matrix(0, n, h)
  mse <- matrix(0, n, n)
  for (s in 1:h) {
    Ps <- matrix(Psi[s, , ], n, n)
    mse <- mse + Ps %*% Sigma_e %*% t(Ps)
    fc_se[, s] <- sqrt(diag(mse))
  }

  rownames(fc_mean) <- rownames(fc_se) <- varnames
  colnames(fc_mean) <- colnames(fc_se) <- paste0("h=", 1:h)

  z <- stats::qnorm(1 - (1 - level / 100) / 2)
  result <- list(
    mean = fc_mean,
    se = fc_se,
    lower = fc_mean - z * fc_se,
    upper = fc_mean + z * fc_se,
    horizon = h,
    level = level,
    varnames = varnames
  )

  class(result) <- "rbfmvar_forecast"
  result
}

#' @title Print Method for rbfmvar_forecast Objects
#' @description Prints a summary of an RBFM-VAR forecast.
#' @param x An \code{rbfmvar_forecast} object.
#' @param ... Additional arguments (currently ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.rbfmvar_forecast <- function(x, ...) {
  cat("\nRBFM-VAR Forecast\n")
  cat("=================\n\n")
  cat("Horizon:", x$horizon, "periods\n")
  cat("Confidence level:", x$level, "%\n\n")

  cat("Point Forecasts:\n")
  print(round(x$mean, 4))

  cat("\nStandard Errors:\n")
  print(round(x$se, 4))

  invisible(x)
}

#' @title Plot Method for rbfmvar_forecast Objects
#' @description Plots forecasts from an RBFM-VAR model.
#' @param x An \code{rbfmvar_forecast} object.
#' @param ... Additional arguments passed to \code{plot}.
#' @return No return value, called for side effects (produces a plot).
#' @export
plot.rbfmvar_forecast <- function(x, ...) {
  n <- length(x$varnames)
  h <- x$horizon

  # Set up plotting grid
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  n_row <- ceiling(sqrt(n))
  n_col <- ceiling(n / n_row)
  graphics::par(mfrow = c(n_row, n_col), mar = c(4, 4, 3, 1))

  horizons <- 1:h

  for (i in 1:n) {
    ylim <- range(c(x$lower[i, ], x$upper[i, ]))

    plot(horizons, x$mean[i, ], type = "b", pch = 19,
         xlab = "Horizon", ylab = "Forecast",
         main = x$varnames[i], ylim = ylim, ...)

    graphics::polygon(c(horizons, rev(horizons)),
                      c(x$lower[i, ], rev(x$upper[i, ])),
                      col = grDevices::rgb(0, 0, 1, 0.2), border = NA)

    graphics::abline(h = 0, lty = 2, col = "gray")
  }
}
