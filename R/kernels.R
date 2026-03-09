#' @title Kernel Functions for Long-Run Variance Estimation
#' @name kernels
#' @description Kernel weight functions for heteroskedasticity and autocorrelation
#'   consistent (HAC) long-run variance estimation.
#'
#' @details
#' These kernel functions are used in the estimation of long-run variance matrices
#' following Andrews (1991) <doi:10.2307/2938229>. The kernels satisfy the conditions
#' for consistent LRV estimation under weak dependence.
#'
#' @references
#' Andrews, D. W. K. (1991). Heteroskedasticity and Autocorrelation Consistent
#' Covariance Matrix Estimation. \emph{Econometrica}, 59(3), 817-858.
#' \doi{10.2307/2938229}
#'
#' Newey, W. K., & West, K. D. (1987). A Simple, Positive Semi-Definite,
#' Heteroskedasticity and Autocorrelation Consistent Covariance Matrix.
#' \emph{Econometrica}, 55(3), 703-708. \doi{10.2307/1913610}
#'
#' @keywords internal
NULL

#' Bartlett (Newey-West) Kernel
#'
#' @param x Numeric scalar or vector of kernel arguments.
#' @return Kernel weights.
#' @keywords internal
kernel_bartlett <- function(x) {
  x <- abs(x)
  ifelse(x < 1, 1 - x, 0)
}

#' Parzen Kernel
#'
#' @param x Numeric scalar or vector of kernel arguments.
#' @return Kernel weights.
#' @keywords internal
kernel_parzen <- function(x) {
  x <- abs(x)
  ifelse(x <= 0.5,
         1 - 6 * x^2 + 6 * x^3,
         ifelse(x <= 1,
                2 * (1 - x)^3,
                0))
}

#' Quadratic Spectral (QS) Kernel
#'
#' @param x Numeric scalar or vector of kernel arguments.
#' @return Kernel weights.
#' @keywords internal
kernel_qs <- function(x) {
  # Avoid division by zero
  ifelse(abs(x) < 1e-10,
         1,
         {
           z <- 6 * pi * x / 5
           25 / (12 * pi^2 * x^2) * (sin(z) / z - cos(z))
         })
}

#' Get Kernel Function by Name
#'
#' @param kernel Character string specifying the kernel type.
#' @return Kernel function.
#' @keywords internal
get_kernel_function <- function(kernel) {
  kernel <- tolower(kernel)
  switch(kernel,
         bartlett = kernel_bartlett,
         parzen = kernel_parzen,
         qs = kernel_qs,
         stop("Unknown kernel: ", kernel, ". Use 'bartlett', 'parzen', or 'qs'."))
}

#' Get Kernel Characteristic Exponent
#'
#' @description Returns the characteristic exponent q used in Andrews (1991)
#'   automatic bandwidth selection.
#'
#' @param kernel Character string specifying the kernel type.
#' @return Characteristic exponent q.
#' @keywords internal
get_kernel_exponent <- function(kernel) {
  kernel <- tolower(kernel)
  switch(kernel,
         bartlett = 1,
         parzen = 2,
         qs = 2,
         stop("Unknown kernel: ", kernel))
}

#' Get Kernel Rate Constant
#'
#' @description Returns the rate constant c_gamma used in Andrews (1991)
#'   automatic bandwidth selection.
#'
#' @param kernel Character string specifying the kernel type.
#' @return Rate constant.
#' @keywords internal
get_kernel_constant <- function(kernel) {
  kernel <- tolower(kernel)
  switch(kernel,
         bartlett = 1.1447,
         parzen = 2.6614,
         qs = 1.3221,
         stop("Unknown kernel: ", kernel))
}
