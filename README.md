# rbfmvar

**Residual-Based Fully Modified VAR Estimation**

[![CRAN status](https://www.r-pkg.org/badges/version/rbfmvar)](https://CRAN.R-project.org/package=rbfmvar)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

`rbfmvar` implements the **Residual-Based Fully Modified VAR (RBFM-VAR)**
estimator of Chang (2000) for vector autoregressive models containing an
**unknown mixture of I(0), I(1), and I(2) variables**.

The procedure:
1. Rewrites the VAR(p) in ECM form with second differences on the LHS.
2. Estimates OLS coefficients and residuals.
3. Applies a nonparametric **long-run variance correction** (Bartlett, Parzen, or QS kernel with automatic Andrews 1991 bandwidth selection).
4. Produces FM-corrected estimates with mixed-normal limiting distributions under any integration mix.

## Installation

```r
install.packages("rbfmvar")
```

## Quick Start

```r
library(rbfmvar)

set.seed(42)
n  <- 100
y1 <- cumsum(rnorm(n))
y2 <- 0.5 * y1 + rnorm(n, sd = 0.5)
Y  <- cbind(y1, y2)

res <- rbfmvar(Y, lags = 2)
print(res)

# With Granger causality test
res2 <- rbfmvar(cbind(y1 = y1, y2 = y2), lags = 2, granger = "y1 -> y2")
print(res2$granger)
```

## Reference

Chang, Y. (2000). Vector autoregressions with unknown mixtures of I(0), I(1), and I(2) components. *Econometric Theory*, 16(6), 905–926. <https://doi.org/10.1017/S0266466600166017>

## Author

Muhammad Alkhalaf <muhammedalkhalaf@gmail.com>
