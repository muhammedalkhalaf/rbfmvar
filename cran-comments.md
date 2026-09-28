## rbfmvar 2.1.0

This release (the current CRAN version is 2.0.2) corrects the computations below.

* The estimator now implements the RBFM-VAR correction of Chang (2000, equations 10 to 13). Previous versions used the OLS residuals in place of the process v_t of equation (11), computed the long-run covariance of the residuals instead of the cross covariance with v_t, and spread the row means of a correction matrix over all nonstationary regressors, so the reported "FM+" coefficients were not the RBFM-VAR estimates. The correction is now (Y'Z, Y+'W + T Delta+)(X'X)^{-1} with Y+' = Y' - Omega_ev Omega_vv^{-1} V', Delta+ = Omega_ev Omega_vv^{-1} Delta_vDw, v_t = (D2 y_{t-1}, D y_{t-1} - N D y_{t-2}) and Dw_t = (D2 y_{t-1}, D y_{t-1}). A test reproduces the formula directly. With the data generating process of Chang (2000, Section 5), the OLS biases match Table 1 of the paper and the RBFM-VAR biases have the same signs and orders of magnitude; the paper does not state the kernel or bandwidth it used, so exact replication is not possible.
* The estimation sample now starts at t = max(p + 1, 4) (the levels VAR(p) needs y_{t-p} and v_t needs D y_{t-2}).
* New internal `estimate_cross_lrv()`; `estimate_onesided_lrv()` sums k(j/K) T^{-1} sum_t a_{t+j} b_t' as in the convention of Chang (2000, p. 909); the Quadratic Spectral kernel now uses all lags (its weights do not vanish beyond the bandwidth). The automatic bandwidth is Andrews (1991) applied to v_t.
* `granger_test()`: the modified Wald statistic ignored the covariances between the tested coefficients (it summed squared t ratios). It now uses Sigma_e[i, i] (X'X)^{-1} for the tested block (Chang, 2000, eq. 20); the chi-square p-value is conservative (Theorem 2).
* `irf()`: the responses used only Pi_1 as a VAR(1) matrix, and the "bootstrap" intervals added random noise to the point estimates. The responses are now those of the levels VAR implied by the estimates, and the intervals come from a recursive-design residual bootstrap that re-estimates the model. `fevd()` inherits the correction.
* `forecast()`: forecasts were Pi_1 times the last fitted second difference. They are now level forecasts from the implied levels VAR, with the forecast error variance sum Psi_i Sigma_e Psi_i'.
* `ic_table()` returned nothing; it now reports the criteria from the stored data. Lag selection compares all lag orders on a common sample and counts n p parameters per equation (there is no constant).
* The fitted object now also contains `data`, `XtX_inv` and `N_hat`.

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
