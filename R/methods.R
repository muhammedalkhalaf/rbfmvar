#' @title Print Method for rbfmvar Objects
#'
#' @description Prints a summary of an RBFM-VAR estimation.
#'
#' @param x An \code{rbfmvar} object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.rbfmvar <- function(x, ...) {
  cat("\n")
  cat("Residual-Based Fully Modified VAR (RBFM-VAR)\n")
  cat("=============================================\n")
  cat("Chang, Y. (2000). Econometric Theory, 16(6), 905-926.\n\n")

  cat("Variables:     ", paste(x$varnames, collapse = ", "), "\n")
  cat("VAR order (p): ", x$p_lags, "\n")
  cat("Observations:  ", x$nobs, "\n")
  cat("Effective T:   ", x$T_eff, "\n")
  cat("Kernel:        ", format_kernel(x$kernel), "\n")
  cat("Bandwidth:     ", sprintf("%.2f", x$bandwidth), "\n")

  if (x$ic != "none") {
    cat("Lag selection: ", toupper(x$ic), "\n")
  }

  cat("\nUse summary() for detailed coefficient tables.\n")

  invisible(x)
}

#' @title Summary Method for rbfmvar Objects
#'
#' @description Provides detailed summary of RBFM-VAR estimation results.
#'
#' @param object An \code{rbfmvar} object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list of class \code{"summary.rbfmvar"} containing summary information.
#'
#' @export
summary.rbfmvar <- function(object, ...) {
  n <- object$n_vars
  varnames <- object$varnames
  level <- object$level

  F_plus <- object$F_plus
  SE_mat <- object$SE_mat
  regnames <- object$regnames

  # Build coefficient table for each equation
  coef_tables <- list()

  cv <- stats::qnorm(1 - (100 - level) / 200)

  for (i in 1:n) {
    coef_i <- F_plus[, i]
    se_i <- SE_mat[i, ]
    z_i <- coef_i / se_i
    p_i <- 2 * (1 - stats::pnorm(abs(z_i)))
    ci_lo <- coef_i - cv * se_i
    ci_hi <- coef_i + cv * se_i

    tbl <- data.frame(
      Estimate = coef_i,
      Std.Error = se_i,
      z.value = z_i,
      Pr.z = p_i,
      CI.lower = ci_lo,
      CI.upper = ci_hi
    )
    rownames(tbl) <- regnames

    coef_tables[[varnames[i]]] <- tbl
  }

  result <- list(
    call = object$call,
    varnames = varnames,
    n_vars = n,
    p_lags = object$p_lags,
    nobs = object$nobs,
    T_eff = object$T_eff,
    kernel = object$kernel,
    bandwidth = object$bandwidth,
    ic = object$ic,
    level = level,
    coefficients = coef_tables,
    Sigma_e = object$Sigma_e
  )

  class(result) <- "summary.rbfmvar"
  result
}

#' @export
print.summary.rbfmvar <- function(x, digits = 4, ...) {
  cat("\n")
  cat("==============================================================\n")
  cat("   Residual-Based Fully Modified VAR (RBFM-VAR)\n")
  cat("   Chang, Y. (2000). Econometric Theory, 16(6), 905-926.\n")
  cat("==============================================================\n\n")

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Variables:     ", paste(x$varnames, collapse = ", "), "\n")
  cat("VAR order (p): ", x$p_lags, "\n")
  cat("Observations:  ", x$nobs, "\n")
  cat("Effective T:   ", x$T_eff, "\n")
  cat("Kernel:        ", format_kernel(x$kernel), "\n")
  cat("Bandwidth:     ", sprintf("%.2f", x$bandwidth), "\n")

  if (x$ic != "none") {
    cat("Lag selection: ", toupper(x$ic), "\n")
  }

  cat("\n")

  # Print coefficient tables for each equation
  for (eq in x$varnames) {
    cat("--------------------------------------------------------------\n")
    cat("Equation:", eq, "\n")
    cat("--------------------------------------------------------------\n")

    tbl <- x$coefficients[[eq]]

    # Add significance stars
    stars <- character(nrow(tbl))
    stars[tbl$Pr.z < 0.01] <- "***"
    stars[tbl$Pr.z >= 0.01 & tbl$Pr.z < 0.05] <- "**"
    stars[tbl$Pr.z >= 0.05 & tbl$Pr.z < 0.10] <- "*"
    stars[tbl$Pr.z >= 0.10] <- ""

    # Print formatted
    cat(sprintf("%16s %11s %11s %9s %9s    [%d%% Conf. Interval]\n",
                "", "Coef.", "Std.Err.", "z", "P>|z|", x$level))
    cat("----------------+-----------------------------------------------------------\n")

    for (j in 1:nrow(tbl)) {
      cat(sprintf("%16s | %10.6f %10.6f %8.3f %8.4f%s %10.6f  %10.6f\n",
                  substr(rownames(tbl)[j], 1, 16),
                  tbl$Estimate[j],
                  tbl$Std.Error[j],
                  tbl$z.value[j],
                  tbl$Pr.z[j],
                  sprintf("%-3s", stars[j]),
                  tbl$CI.lower[j],
                  tbl$CI.upper[j]))
    }

    cat("\n")
  }

  cat("--------------------------------------------------------------\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05\n\n")

  cat("Error Covariance Matrix (Sigma_e):\n")
  print(round(x$Sigma_e, digits))

  cat("\nNote: Standard errors from Var(vec(F+')) = Sigma_e (x) (Z'Z)^-1.\n")
  cat("      P-values are conservative for nonstationary regressors.\n")

  invisible(x)
}

#' Extract Coefficients from rbfmvar Object
#'
#' @param object An \code{rbfmvar} object.
#' @param type Character. Type of coefficients to extract: \code{"plus"}
#'   (FM+ corrected, default) or \code{"ols"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Coefficient matrix.
#'
#' @export
coef.rbfmvar <- function(object, type = "plus", ...) {
  type <- match.arg(type, c("plus", "ols"))

  if (type == "plus") {
    object$F_plus
  } else {
    object$F_ols
  }
}

#' Extract Residuals from rbfmvar Object
#'
#' @param object An \code{rbfmvar} object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Matrix of residuals.
#'
#' @export
residuals.rbfmvar <- function(object, ...) {
  object$residuals
}

#' Extract Fitted Values from rbfmvar Object
#'
#' @param object An \code{rbfmvar} object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Matrix of fitted values.
#'
#' @export
fitted.rbfmvar <- function(object, ...) {
  object$fitted
}

#' Extract Variance-Covariance Matrix from rbfmvar Object
#'
#' @param object An \code{rbfmvar} object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Error covariance matrix.
#'
#' @export
vcov.rbfmvar <- function(object, ...) {
  object$Sigma_e
}

#' Format Kernel Name for Display
#'
#' @param kernel Kernel name.
#' @return Formatted kernel name.
#' @keywords internal
format_kernel <- function(kernel) {
  switch(tolower(kernel),
         bartlett = "Bartlett",
         parzen = "Parzen",
         qs = "Quadratic Spectral",
         kernel)
}
