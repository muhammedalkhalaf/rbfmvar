# rbfmvar: Residual-Based Fully Modified Vector Autoregression

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/rbfmvar)](https://CRAN.R-project.org/package=rbfmvar)
<!-- badges: end -->

## Overview

**rbfmvar** implements the Residual-Based Fully Modified Vector Autoregression
(RBFM-VAR) estimator following [Chang (2000)](https://doi.org/10.1017/S0266466600166071).
The RBFM-VAR procedure extends Phillips (1995) FM-VAR to handle any unknown
mixture of I(0), I(1), and I(2) components without prior knowledge of the
number or location of unit roots.

## Key Features

- **Robust to unknown integration orders**: Handles I(0), I(1), and I(2) variables automatically
- **Asymptotically valid inference**: Chi-squared Wald statistics for hypothesis testing
- **Automatic lag selection**: Via AIC, BIC, or HQ information criteria
- **Multiple kernels for LRV estimation**: Bartlett, Parzen, Quadratic Spectral
- **Andrews (1991) automatic bandwidth selection**
- **Granger non-causality testing**
- **Impulse response functions** with bootstrap confidence intervals
- **Forecast error variance decomposition**
- **Out-of-sample forecasting**

## Installation

Install from CRAN:

```r
install.packages("rbfmvar")
```

Or install the development version from GitHub:

```r
# install.packages("devtools")
```

## Usage

### Basic Estimation

```r
library(rbfmvar)

# Simulate some VAR data
set.seed(123)
n <- 200
y <- matrix(0, n, 3)
colnames(y) <- c("gdp", "inflation", "interest")
e <- matrix(rnorm(n * 3), n, 3)

for (t in 3:n) {
  y[t, ] <- 0.5 * y[t-1, ] + 0.2 * y[t-2, ] + e[t, ]
}

# Estimate RBFM-VAR with automatic lag selection
fit <- rbfmvar(y, max_lags = 6, ic = "bic", kernel = "bartlett")
summary(fit)
```

### Granger Causality Testing

```r
# Test if gdp Granger-causes inflation
granger_test(fit, cause = "gdp", effect = "inflation")

# Pairwise Granger causality matrix
granger_matrix(fit)
```

### Impulse Response Functions

```r
# Compute IRF with bootstrap confidence intervals
ir <- irf(fit, horizon = 20, boot = 500, ci = 95)
plot(ir)
```

### Forecast Error Variance Decomposition

```r
fv <- fevd(fit, horizon = 20)
plot(fv)
```

### Forecasting

```r
fc <- forecast(fit, h = 12)
plot(fc)
```

## Methodology

The RBFM-VAR model is based on Chang (2000), which develops a fully modified
VAR estimation procedure that is robust to unknown integration orders. The key
innovation is using second differences to eliminate I(2) trends while applying
FM corrections to handle endogeneity from I(1) regressors.

The model is specified as:

$$\Delta^2 y_t = \sum_{j=1}^{p-2} \Gamma_j \Delta^2 y_{t-j} + \Pi_1 \Delta y_{t-1} + \Pi_2 y_{t-1} + e_t$$

The FM+ estimator corrects for the asymptotic bias from the correlation between
regression errors and innovations in integrated regressors, achieving:

- Zero mean mixed normal limiting distribution
- Chi-square Wald statistics for hypothesis testing
- Consistent estimation regardless of integration orders

## References

- Chang, Y. (2000). Vector Autoregressions with Unknown Mixtures of I(0), I(1),
  and I(2) Components. *Econometric Theory*, 16(6), 905-926.
  [doi:10.1017/S0266466600166071](https://doi.org/10.1017/S0266466600166071)

- Phillips, P. C. B. (1995). Fully Modified Least Squares and Vector
  Autoregression. *Econometrica*, 63(5), 1023-1078.
  [doi:10.2307/2171721](https://doi.org/10.2307/2171721)

- Andrews, D. W. K. (1991). Heteroskedasticity and Autocorrelation Consistent
  Covariance Matrix Estimation. *Econometrica*, 59(3), 817-858.
  [doi:10.2307/2938229](https://doi.org/10.2307/2938229)

## License

GPL (>= 3)

## Author

