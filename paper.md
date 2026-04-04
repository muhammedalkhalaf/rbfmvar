---
title: 'rbfmvar: Residual-Based Fully Modified Vector Autoregression'
tags:
  - R
  - econometrics
  - VAR
  - cointegration
  - Granger causality
  - impulse response
  - time series
authors:
  - name: Muhammad Abdullah Alkhalaf
    orcid: 0009-0002-2677-9246
    corresponding: true
    email: muhammedalkhalaf@gmail.com
    affiliation: 1
affiliations:
  - name: Rufyq Elngeh for Academic and Business Services, Riyadh, Saudi Arabia
    index: 1
date: 25 March 2026
bibliography: paper.bib
---

# Summary

`rbfmvar` implements the Residual-Based Fully Modified Vector Autoregression (RBFM-VAR) estimator of @Chang2000. The RBFM-VAR procedure extends the Phillips [-@Phillips1995] FM-VAR approach to handle any unknown mixture of I(0), I(1), and I(2) components without requiring prior knowledge of the number or location of unit roots. This robustness is achieved by applying a fully modified correction to VAR residuals rather than to levels, yielding asymptotically chi-squared Wald statistics for Granger non-causality tests regardless of integration order. The package provides automatic lag selection via information criteria (AIC, BIC, HQ), long-run variance estimation using Bartlett, Parzen, or Quadratic Spectral kernels with @Andrews1991 automatic bandwidth selection, impulse response functions (IRF) with bootstrap confidence intervals, forecast error variance decomposition (FEVD), and out-of-sample forecasting.

# Statement of Need

Standard VAR inference breaks down when the system contains integrated or cointegrated variables: OLS Wald statistics for Granger causality tests have non-standard distributions, and pre-tests for integration order are required before inference. The Toda-Yamamoto [-@TodaYamamoto1995] approach augments the VAR lag length but requires knowledge of the maximum integration order. The RBFM-VAR of @Chang2000 sidesteps both issues by correcting for serial correlation and endogeneity in the residuals, producing standard chi-squared inference without pre-testing. While RBFM-VAR implementations exist in GAUSS and EViews, no R package has provided this estimator. `rbfmvar` fills this gap, enabling robust causal inference in multivariate time series with mixed integration orders entirely within R.

# Usage

## Fitting the RBFM-VAR Model

```r
library(rbfmvar)

# Fit RBFM-VAR with automatic lag selection
fit <- rbfmvar(data = macro_ts, lag_max = 4, ic = "BIC",
               kernel = "Bartlett")
summary(fit)
print(fit)
```

## Granger Non-Causality Tests

```r
# Wald test: does x Granger-cause y? (asymptotic chi-squared)
gc_result <- granger(fit, cause = "x", effect = "y")
print(gc_result)
```

## Impulse Response Functions

```r
# IRF with bootstrap confidence intervals
irf_result <- irf(fit, impulse = "x", response = "y",
                  h = 20, nboot = 1000, ci = 0.95)
plot(irf_result)
```

## Forecast Error Variance Decomposition

```r
fevd_result <- fevd(fit, h = 20)
print(fevd_result)
plot(fevd_result)
```

## Out-of-Sample Forecasting

```r
fc <- forecast(fit, h = 8)
print(fc)
```

# Implementation

`rbfmvar` is implemented in R with a dependency on `MASS` for matrix operations. Lag selection in `lag_selection()` evaluates AIC, BIC, and Hannan-Quinn criteria over a grid of candidate orders. Long-run variance estimation in `lrv()` uses the `kernels()` module which implements Bartlett, Parzen, and Quadratic Spectral kernels with @Andrews1991 data-driven bandwidth. The core `rbfmvar()` function estimates the VAR by OLS, extracts residuals, applies the FM correction following @Chang2000, and constructs the RBFM coefficient matrix. Granger causality Wald statistics in `granger()` are asymptotically chi-squared under the null regardless of integration order. IRFs are computed via the moving-average representation and bootstrapped using a residual bootstrap with 1000 replications by default. FEVD decomposition in `fevd()` uses the Cholesky factorisation of the residual covariance matrix.

# References
