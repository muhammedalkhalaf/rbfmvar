# rbfmvar 2.0.0

## Initial Release

* Implements Residual-Based Fully Modified VAR (RBFM-VAR) estimator following
  Chang (2000).

* Core estimation:
  - `rbfmvar()`: Main estimation function for RBFM-VAR models.
  - Handles unknown mixtures of I(0), I(1), and I(2) variables.
  - FM+ bias correction for asymptotically valid inference.

* Lag selection:
  - Automatic lag selection via AIC, BIC, or HQ information criteria.
  - `ic_table()`: Display information criteria comparison.

* Long-run variance estimation:
  - Bartlett (Newey-West), Parzen, and Quadratic Spectral kernels.
  - Andrews (1991) automatic bandwidth selection.

* Inference:
  - `granger_test()`: Granger non-causality testing with modified Wald statistics.
  - `granger_matrix()`: Pairwise Granger causality tests.
  - Asymptotically chi-squared inference regardless of integration orders.

* Impulse response analysis:
  - `irf()`: Orthogonalized impulse response functions.
  - Bootstrap confidence intervals via Kilian (1998) method.

* Forecast error variance decomposition:
  - `fevd()`: Cholesky-identified variance decomposition.

* Forecasting:
  - `forecast()`: Out-of-sample forecasting with prediction intervals.

* Methods:
  - `print()`, `summary()`, `plot()` methods for all major objects.
  - `coef()`, `residuals()`, `fitted()`, `vcov()` extractors.

## References

Chang, Y. (2000). Vector Autoregressions with Unknown Mixtures of I(0), I(1),
and I(2) Components. *Econometric Theory*, 16(6), 905-926.
doi:10.1017/S0266466600166071
