#' @title Forecast Error Variance Decomposition
#'
#' @description Computes the forecast error variance decomposition (FEVD) from
#'   an RBFM-VAR model.
#'
#' @param object An \code{rbfmvar} object from \code{\link{rbfmvar}}.
#' @param horizon Integer. Number of periods for the FEVD. Default is 20.
#'
#' @details
#' The FEVD shows the proportion of the forecast error variance of each
#' variable that is attributable to shocks in each of the structural
#' innovations. The decomposition is based on the Cholesky identification
#' scheme, so the ordering of variables matters.
#'
#' At each horizon h, the FEVD sums to 1 (100%) for each variable.
#'
#' @return An object of class \code{"rbfmvar_fevd"} containing:
#' \describe{
#'   \item{fevd}{Array of FEVD values (horizon x n x n). Element [h, i, j] is
#'     the proportion of variable i's forecast error variance at horizon h
#'     explained by shocks in variable j.}
#'   \item{horizon}{FEVD horizon.}
#'   \item{varnames}{Variable names.}
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
#' fv <- fevd(fit, horizon = 20)
#' plot(fv)
#'
#' @export
fevd <- function(object, horizon = 20) {
  if (!inherits(object, "rbfmvar")) {
    stop("'object' must be of class 'rbfmvar'.")
  }

  if (horizon < 1) {
    stop("'horizon' must be at least 1.")
  }

  # First compute IRF
  ir <- irf(object, horizon = horizon, ortho = TRUE, boot = 0)

  n <- length(ir$varnames)
  irf_array <- ir$irf

  # Compute FEVD from orthogonalized IRF
  # FEVD[h, i, j] = sum_{s=0}^h (Psi[s, i, j])^2 / sum_{k=1}^n sum_{s=0}^h (Psi[s, i, k])^2

  fevd_array <- array(0, dim = c(horizon + 1, n, n))
  dimnames(fevd_array) <- dimnames(irf_array)

  for (h in 1:(horizon + 1)) {
    for (i in 1:n) {
      # Total forecast error variance for variable i up to horizon h
      total_var <- 0
      for (k in 1:n) {
        total_var <- total_var + sum(irf_array[1:h, i, k]^2)
      }

      # Contribution from each shock
      if (total_var > 0) {
        for (j in 1:n) {
          fevd_array[h, i, j] <- sum(irf_array[1:h, i, j]^2) / total_var
        }
      } else {
        # Equal contribution if total variance is zero
        fevd_array[h, i, ] <- 1 / n
      }
    }
  }

  result <- list(
    fevd = fevd_array,
    horizon = horizon,
    varnames = ir$varnames
  )

  class(result) <- "rbfmvar_fevd"
  result
}

#' @export
print.rbfmvar_fevd <- function(x, ...) {
  cat("\nForecast Error Variance Decomposition (RBFM-VAR)\n")
  cat("================================================\n\n")
  cat("Horizon:", x$horizon, "periods\n")
  cat("Variables:", paste(x$varnames, collapse = ", "), "\n")
  cat("Cholesky ordering:", paste(x$varnames, collapse = " -> "), "\n")
  cat("\nUse plot() to visualize the FEVD.\n")
  cat("\nFEVD at horizon", x$horizon, ":\n\n")

  fevd_final <- x$fevd[x$horizon + 1, , ]
  rownames(fevd_final) <- x$varnames
  colnames(fevd_final) <- x$varnames
  print(round(fevd_final, 4))

  invisible(x)
}

#' @export
plot.rbfmvar_fevd <- function(x, ...) {
  n <- length(x$varnames)
  horizon <- x$horizon

  # Set up plotting grid
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  n_row <- ceiling(sqrt(n))
  n_col <- ceiling(n / n_row)
  graphics::par(mfrow = c(n_row, n_col), mar = c(4, 4, 3, 1))

  h <- 0:horizon
  colors <- grDevices::rainbow(n, alpha = 0.7)

  for (i in 1:n) {
    fevd_i <- t(x$fevd[, i, ])

    graphics::barplot(fevd_i, col = colors, border = NA,
                      main = paste("FEVD:", x$varnames[i]),
                      xlab = "Horizon", ylab = "Proportion",
                      names.arg = h, las = 1, ...)

    if (i == 1) {
      graphics::legend("topright", legend = x$varnames,
                       fill = colors, bty = "n", cex = 0.8)
    }
  }
}
