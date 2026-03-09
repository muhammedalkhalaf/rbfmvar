#' @title Granger Non-Causality Test
#'
#' @description Tests for Granger non-causality in the RBFM-VAR framework using
#'   a modified Wald statistic. The test is asymptotically chi-squared under
#'   the null hypothesis, regardless of the integration order of the variables.
#'
#' @param object An \code{rbfmvar} object from \code{\link{rbfmvar}}.
#' @param cause Character string. Name of the causing variable.
#' @param effect Character string. Name of the affected variable.
#'
#' @details
#' The Granger non-causality hypothesis is:
#' \deqn{H_0: x \text{ does not Granger-cause } y}
#'
#' This is tested by examining whether the coefficients on lagged values of
#' \code{cause} in the equation for \code{effect} are jointly zero.
#'
#' Under the FM+ framework of Chang (2000), the Wald statistic has an
#' asymptotic chi-squared distribution that provides a conservative (valid)
#' p-value even when variables have unknown integration orders.
#'
#' @return A list of class \code{"rbfmvar_granger"} containing:
#' \describe{
#'   \item{cause}{Name of the causing variable.}
#'   \item{effect}{Name of the affected variable.}
#'   \item{statistic}{Modified Wald statistic.}
#'   \item{df}{Degrees of freedom.}
#'   \item{p.value}{P-value (conservative).}
#'   \item{coefficients}{Restricted coefficients being tested.}
#' }
#'
#' @references
#' Chang, Y. (2000). Vector Autoregressions with Unknown Mixtures of I(0), I(1),
#' and I(2) Components. \emph{Econometric Theory}, 16(6), 905-926.
#' \doi{10.1017/S0266466600166071}
#'
#' Toda, H. Y., & Yamamoto, T. (1995). Statistical Inference in Vector
#' Autoregressions with Possibly Integrated Processes. \emph{Journal of
#' Econometrics}, 66(1-2), 225-250. \doi{10.1016/0304-4076(94)01616-8}
#'
#' @examples
#' # Simulate VAR data
#' set.seed(42)
#' n <- 200
#' e <- matrix(rnorm(n * 3), n, 3)
#' y <- matrix(0, n, 3)
#' colnames(y) <- c("x", "y", "z")
#' for (t in 3:n) {
#'   y[t, "x"] <- 0.5 * y[t-1, "x"] + e[t, 1]
#'   y[t, "y"] <- 0.3 * y[t-1, "y"] + 0.4 * y[t-1, "x"] + e[t, 2]
#'   y[t, "z"] <- 0.2 * y[t-1, "z"] + e[t, 3]
#' }
#'
#' fit <- rbfmvar(y, lags = 2)
#'
#' # Test if x Granger-causes y (should be significant)
#' test1 <- granger_test(fit, cause = "x", effect = "y")
#' print(test1)
#'
#' # Test if z Granger-causes y (should not be significant)
#' test2 <- granger_test(fit, cause = "z", effect = "y")
#' print(test2)
#'
#' @export
granger_test <- function(object, cause, effect) {
  if (!inherits(object, "rbfmvar")) {
    stop("'object' must be of class 'rbfmvar'.")
  }

  varnames <- object$varnames
  p <- object$p_lags

  # Validate variable names
  if (!cause %in% varnames) {
    stop("'cause' variable '", cause, "' not found. ",
         "Available: ", paste(varnames, collapse = ", "))
  }

  if (!effect %in% varnames) {
    stop("'effect' variable '", effect, "' not found. ",
         "Available: ", paste(varnames, collapse = ", "))
  }

  if (cause == effect) {
    stop("'cause' and 'effect' must be different variables.")
  }

  # Get coefficient matrix and standard errors
  F_plus <- object$F_plus
  SE_mat <- object$SE_mat

  # Find equation for 'effect' variable
  effect_idx <- which(varnames == effect)
  cause_idx <- which(varnames == cause)

  # Find regressors involving 'cause' variable
  # These are: L.cause (in Pi2) and LD.cause (in Pi1)
  # And for p >= 3: L1D2.cause, L2D2.cause, etc.

  regnames <- object$regnames
  n <- object$n_vars

  # Indices of 'cause' regressors
  cause_regs <- grep(paste0("\\.", cause, "$"), regnames)

  if (length(cause_regs) == 0) {
    stop("No regressors found for '", cause, "'.")
  }

  # Extract coefficients and SEs for the restriction
  coefs <- F_plus[cause_regs, effect_idx]
  ses <- SE_mat[effect_idx, cause_regs]

  # Build restriction matrix R: we test R * beta = 0
  # For standard Granger test, R is identity selecting cause coefficients
  n_restrictions <- length(coefs)

  # Compute Wald statistic
  # W = beta' * (R * Var(beta) * R')^{-1} * beta

  # Variance of selected coefficients
  # From Var(vec(F+')) = Sigma_e (x) (Z'Z)^{-1}
  # Diagonal elements give variances
  var_coefs <- ses^2

  # Handle potential zero or near-zero variances
  var_coefs <- pmax(var_coefs, .Machine$double.eps)

  # For diagonal covariance, Wald = sum(beta^2 / var(beta))
  wald_stat <- sum(coefs^2 / var_coefs)

  # P-value from chi-squared distribution
  df <- n_restrictions
  p_value <- 1 - stats::pchisq(wald_stat, df = df)

  result <- list(
    cause = cause,
    effect = effect,
    statistic = wald_stat,
    df = df,
    p.value = p_value,
    coefficients = coefs,
    se = ses,
    hypothesis = paste0(cause, " does not Granger-cause ", effect)
  )

  class(result) <- "rbfmvar_granger"
  result
}

#' @export
print.rbfmvar_granger <- function(x, ...) {
  cat("\n")
  cat("Granger Non-Causality Test (Modified Wald)\n")
  cat("==========================================\n\n")
  cat("H0:", x$hypothesis, "\n\n")
  cat("Modified Wald statistic:", sprintf("%.4f", x$statistic), "\n")
  cat("Degrees of freedom:     ", x$df, "\n")
  cat("P-value (conservative): ", sprintf("%.4f", x$p.value), "\n\n")

  # Decision
  if (x$p.value < 0.01) {
    cat("Decision: Reject H0 at 1% (strong evidence of Granger causality)\n")
  } else if (x$p.value < 0.05) {
    cat("Decision: Reject H0 at 5% (evidence of Granger causality)\n")
  } else if (x$p.value < 0.10) {
    cat("Decision: Reject H0 at 10% (weak evidence of Granger causality)\n")
  } else {
    cat("Decision: Cannot reject H0 (no evidence of Granger causality)\n")
  }

  cat("\nNote: P-value is conservative (bounded above by chi-sq)\n")
  invisible(x)
}

#' Granger Causality Matrix
#'
#' @description Computes pairwise Granger causality tests for all variable pairs.
#'
#' @param object An \code{rbfmvar} object.
#'
#' @return A matrix of p-values for all pairwise Granger causality tests.
#'   Row i, column j contains the p-value for "variable j causes variable i".
#'
#' @examples
#' \donttest{
#' # Generate example data
#' set.seed(123)
#' n <- 100
#' mydata <- data.frame(x = cumsum(rnorm(n)), y = cumsum(rnorm(n)))
#' fit <- rbfmvar(mydata, lags = 2)
#' granger_matrix(fit)
#' }
#'
#' @export
granger_matrix <- function(object) {
  if (!inherits(object, "rbfmvar")) {
    stop("'object' must be of class 'rbfmvar'.")
  }

  varnames <- object$varnames
  n <- length(varnames)

  pval_mat <- matrix(NA, n, n)
  rownames(pval_mat) <- colnames(pval_mat) <- varnames

  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        test <- granger_test(object, cause = varnames[j], effect = varnames[i])
        pval_mat[i, j] <- test$p.value
      }
    }
  }

  class(pval_mat) <- c("granger_matrix", "matrix")
  pval_mat
}

#' @export
print.granger_matrix <- function(x, digits = 4, ...) {
  cat("\nGranger Causality P-value Matrix\n")
  cat("================================\n")
  cat("Entry [i,j] = p-value for 'j causes i'\n\n")

  # Format with significance stars
  x_print <- x
  x_print[is.na(x_print)] <- NA
  stars <- matrix("", nrow(x), ncol(x))
  stars[x < 0.01] <- "***"
  stars[x >= 0.01 & x < 0.05] <- "**"
  stars[x >= 0.05 & x < 0.10] <- "*"

  for (i in 1:nrow(x)) {
    cat(sprintf("%-10s", rownames(x)[i]))
    for (j in 1:ncol(x)) {
      if (is.na(x[i, j])) {
        cat(sprintf("%10s", "-"))
      } else {
        cat(sprintf("%7.4f%s", x[i, j], stars[i, j]))
      }
    }
    cat("\n")
  }

  cat("\n*** p<0.01, ** p<0.05, * p<0.10\n")
  invisible(x)
}
