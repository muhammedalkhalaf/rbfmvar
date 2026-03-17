\name{rbfmvar}
\alias{rbfmvar}
\title{Residual-Based Fully Modified VAR Estimation}
\usage{
rbfmvar(Y, lags = 2L, max_lags = 8L,
  ic = c("none", "aic", "bic", "hq"),
  kernel = c("bartlett", "parzen", "qs"),
  bandwidth = -1, granger = NULL)
}
\arguments{
\item{Y}{Numeric matrix. Each column is a time-series variable
  (T rows x n columns).}

\item{lags}{Integer. VAR lag order \eqn{p}. Default is \code{2}.}

\item{max_lags}{Integer. Maximum lag searched when \code{ic} is not
  \code{"none"}. Default is \code{8}.}

\item{ic}{Character. Information criterion for lag selection:
  \code{"aic"}, \code{"bic"}, \code{"hq"}, or \code{"none"} (default).}

\item{kernel}{Character. Kernel for long-run variance estimation:
  \code{"bartlett"} (default), \code{"parzen"}, or \code{"qs"}.}

\item{bandwidth}{Numeric. Bandwidth parameter. Use \code{-1} (default)
  for automatic Andrews (1991) selection.}

\item{granger}{Character or \code{NULL}. Granger causality specification
  in the form \code{"y1 -> y2"}. Default is \code{NULL}.}
}
\value{
A list of class \code{"rbfmvar"} containing: \code{F_ols},
\code{F_plus}, \code{SE_mat}, \code{Sigma_e}, \code{Pi1_ols},
\code{Pi2_ols}, \code{Pi1_plus}, \code{Pi2_plus}, \code{Gamma_ols},
\code{Omega_ev}, \code{Omega_vv}, \code{Delta_vdw}, \code{bandwidth},
\code{p_lags}, \code{nobs}, \code{n_vars}, \code{T_eff}, \code{kernel},
\code{granger}, \code{irf}, and \code{residuals}.
}
\description{
Estimates a Residual-Based Fully Modified VAR (RBFM-VAR) that handles an
unknown mixture of I(0), I(1), and I(2) variables without prior knowledge
of integration orders (Chang, 2000).
}
\details{
The RBFM-VAR rewrites the VAR(\eqn{p}) in ECM-like form with second
differences on the LHS. OLS estimates are corrected using a nonparametric
long-run variance estimate (Bartlett, Parzen, or quadratic spectral kernel
with automatic Andrews 1991 bandwidth selection), yielding mixed-normal
limiting distributions for the corrected coefficient matrix \eqn{F^+}.
}
\examples{
set.seed(42)
n  <- 80
y1 <- cumsum(rnorm(n))
y2 <- 0.5 * y1 + rnorm(n, sd = 0.5)
Y  <- cbind(y1, y2)
res <- rbfmvar(Y, lags = 2)
print(res)
}
\references{
Chang, Y. (2000). Vector autoregressions with unknown mixtures of I(0),
I(1), and I(2) components. \emph{Econometric Theory}, 16(6), 905--926.
\doi{10.1017/S0266466600166017}

Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation
consistent covariance matrix estimation. \emph{Econometrica}, 59(3),
817--858.
}
