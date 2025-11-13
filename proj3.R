# Practical 3 — Smooth deconvolution
#
# Yunhan Zhang: s2176155. Xiyu Wu: s2799746. Tianyu Wang: s2794991
# Yunhan did question 2&3. Xiyu did question 4&5. Tianyu did question 1&6.
# Each member did roughly the same amount of work.
# Github repo: https://github.com/2413605283/ESP_37.git (branch: practice-3)

# Overview
#   This R file performs smooth deconvolution of England COVID-19 daily
#   death data to recover the underlying infection curve over time.
#
#   The model represents the latent infection curve f(t) with a smooth
#   B-spline basis and links it to daily deaths through a realistic
#   infection-to-death delay distribution. An extended time grid and a
#   convolution structure allow early infections to influence later deaths,
#   and a penalized Poisson likelihood is used to fit the model smoothly.
#
#   The model is fitted by BFGS using an analytically derived gradient,
#   the smoothing parameter is chosen by BIC, and uncertainty in f(t)
#   is assessed using a non-parametric bootstrap. Final plots show the
#   observed and fitted deaths on the real time axis, and the estimated
#   infection curve together with its 95% confidence band, aligned so
#   that the relative timing between infection and death is clearly visible.


library(splines)
engcov <- read.table("engcov.txt", header = TRUE)
#engcov <- read.table("/Users/koo/Desktop/engcov.txt", header = TRUE)


# Task 1

# Function: 
#    evaluate_matrices
#
# Purpose:
#    Build three matrices, tilde(X), X and S, needed for deconvolution model.
#
# Input:
#    K: (integer) number of B-spline basis functions which are used to represent f(t).
#    n: (integer) number of days.
#
# Output:
#    X_tilde: (n + 30)* K dimension matrix, B-spline design matrix for f(t).
#    X: n * K dimension matrix, design matrix that maps the B-spline coefficients to the expected daily deaths.
#    S: K * K dimension matrix, second-difference penalty matrix that is used for smoothing penalty.
#    pi: length-L delay pmf pi(1),...,pi(L).
#
# How it works:
#    1. Define a time grid of length m = n + 30 for the infection curve f(t).
#    2. Place K + 4 spaced knots and build a cubic B-spline basis.
#    3. Construct the infection-to-death delay pmf pi(j), then normalize pd.
#    4. Build a sparse lower-triangular convolution matrix W.
#    4. Form the death-side design matrix X.


evaluate_matrices <- function(K, n) {
  
  # Extended grid length: n observed days and 30 days before first death
  m <- n + 30
  t_extented <- 1:m
  
  # Evaluate evenly-spaced knots and B-spline design matrix
  knots <- seq(1, m, length.out = K + 4)
  X_tilde <- splineDesign(knots, t_extented, outer.ok = TRUE)
  
  # Infection to death delay pmf pi(j)
  d <- 1:80
  edur <- 3.151;
  sdur <- 0.469
  # Evaluate log-normal density valuesand renormalize
  pd <- dlnorm(d, edur, sdur);
  pd <- pd / sum(pd)
  
  # Effective maximum lag on this grid
  nlag <- min(length(pd), m - 1)
  
  # Sparse lower-triangular convolution matrix W
  diags <- lapply(1:nlag, function(ell) rep(pd[ell], m - ell))
  W <- Matrix::bandSparse(m, m, k = -(1:nlag), diagonals = diags)
  
  # Apply convolution to all spline basis columns
  Y_all <- W %*% X_tilde 
  
  # Death-side design matrix
  X <- as.matrix(Y_all[31:(30 + n), , drop = FALSE])
  
  # Second difference penalty matrix
  S <- crossprod(diff(diag(K), differences = 2))
  
  list(X_tilde = X_tilde, X = X, S = S, pi = pd)
}

# Apply Task 1 to get the data
y <- as.numeric(engcov$nhs)
K <- 30
mats <- evaluate_matrices(K = K, n = length(y))
X       <- mats$X
S       <- mats$S
X_tilde <- mats$X_tilde
lambda  <- 5e-5


#Task 2
# Function: nll_gamma()
#
# Purpose:
#   Compute the penalized negative log-likelihood (NLL) for the Poisson deconvolution model.
#   The NLL measures how well the model's predicted deaths (μ) match the observed deaths (y),
#   and plus a smoothness penalty term to avoid overfit.
#
#
# Inputs:
#   gamma:  a numeric vector of length K, log-parameters (γ) used to ensure β = exp(γ) > 0.
#   y:      a numeric vector of observed daily deaths (response variable).
#   X:      n × K model matrix linking infection curve f(t) to daily deaths μ.
#   S:      K × K smoothing penalty matrix (second-difference operator).
#   lambda: scalar smoothing parameter λ controlling the smoothness of f(t).
#
# Output:
#   A single numeric value — the penalized negative log-likelihood to be minimized.
#
# How it works:
#   (1) Compute β = exp(γ) to ensure β > 0.
#   (2) Compute expected deaths μ = Xβ.
#   (3) Evaluate the Poisson NLL (Σ_i [ μ_i − y_i log(μ_i) ]) ignoring constants (log(y_i!)).
#   (4) Compute the smoothness penalty (λ / 2) β^T %*% S %*% β.
#   (5) Return the total NLL = likelihood term + penalty term.
#
nll_gamma <- function(gamma, y, X, S, lambda) {
  
  beta <- exp(gamma)                          # transform parameters to keep β > 0
  mu <- drop(X %*% beta)                      # expected deaths from convolution model
  mu  <- pmax(mu, 1e-12)                      # avoid log(0) or division by zero
  pois <- sum(mu - y * log(mu))               # Poisson negative log-likelihood part
  pen  <- 0.5 * lambda * drop(t(beta) %*% S %*% beta)  # smoothness penalty
  pois + pen                                  # return total penalized NLL
}


# Function: grad_gamma()
#
# Purpose:
#   Compute the exact gradient (first derivative vector) of the penalized
#   negative log-likelihood with respect to γ.
#   This gradient is used by the BFGS optimizer for efficient model fitting.
#
#
# Inputs:
#   gamma:  a numeric vector of log-coefficients γ.
#   y:      numeric vector of observed daily deaths.
#   X:      n × K model matrix linking infections to deaths.
#   S:      K × K penalty matrix (second-difference operator).
#   lambda: scalar smoothing parameter λ controlling penalty strength.
#
# Output:
#   A numeric vector (length K): gradient of the penalized NLL with respect to γ.
#
# How it works:
#   (1) Compute β = exp(γ) and μ = Xβ.
#   (2) Compute gradient w.r.t β:
#       ∇_β L = X^T (1 − y / μ) + λ S β
#   (3) Apply chain rule: ∇_γ L = diag(β) ∇_β L, since dβ/dγ = β.
#   (4) Return the resulting vector.
#
grad_gamma <- function(gamma, y, X, S, lambda) {
  
  beta <- exp(gamma)                          # parameter transformation
  mu  <- drop(X %*% beta)                     # predicted deaths
  mu  <- pmax(mu, 1e-12)
  
  g_beta_poiss <- crossprod(X, 1 - (y / mu))  # Poisson part: X^T (1 − y/μ)
  g_beta_pen   <- lambda * (S %*% beta)       # penalty part: λ S β
  g_beta <- drop(g_beta_poiss) + g_beta_pen   # combine both parts
  
  g_gamma <- beta * g_beta                    # apply chain rule: diag(β) × g_beta
  g_gamma                             
}

# Test gradient using finite differences
gamma0 <- rep(log(1e-3), ncol(S))      # test point for γ
nll0  <- nll_gamma(gamma0, y, X, S, lambda)   # NLL at γ0
g_exact <- grad_gamma(gamma0, y, X, S, lambda) # analytic gradient at γ0

eps <- 1e-5                           # finite-difference step size
fd <- numeric(ncol(S))                # store FD approximations

for (i in 1:ncol(S)) {
  g1 <- gamma0; g1[i] <- g1[i] + eps  # γ0 with +ε at position i
  g2 <- gamma0; g2[i] <- g2[i] - eps  # γ0 with −ε at position i
  # central-difference approximation to ∂L/∂γi
  fd[i] <- (nll_gamma(g1, y, X, S, lambda) -
              nll_gamma(g2, y, X, S, lambda)) / (2 * eps)
}

# compare finite-difference and analytic gradients (first 10 entries)
abs_err <- abs(fd - g_exact)
# simple pass/fail check
if (max(abs_err) < 1) {
  cat("Gradient check: PASS\n")
} else {
  cat("Gradient check: FAIL (max abs error =", max(abs_err), ")\n")
}

## Task 3 (Sanity-check fit)
## ---------------------------------------

day <- engcov$julian  # 时间轴（day of year 2020）

# Objective & gradient (from Task 2)
obj <- function(g) nll_gamma(g, y, X, S, lambda)
gr  <- function(g) grad_gamma(g, y, X, S, lambda)

# ---- Fit with BFGS ----
fit <- optim(gamma0, fn = obj, gr = gr, method = "BFGS")

# Extract estimates
beta_hat <- exp(fit$par)
mu_hat   <- drop(X %*% beta_hat)           # fitted deaths
f_hat    <- drop(X_tilde %*% beta_hat)     # estimated infection curve

# Diagnostics
cat("Convergence:", fit$convergence, "(0 is good)\n")
cat("Final nll:", fit$value, "\n")
cat("||grad||_inf:", max(abs(gr(fit$par))), "\n")

# Time axis for infection curve: from day[1] - 30 to day[n]
day_inf <- day[1] - 31 + seq_along(f_hat)

# ---- Plots ----
op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# (a) Observed vs Fitted deaths (vs real time)
plot(day, y, pch = 16, cex = 0.4, col = "red",
     xlab = "Day of year 2020", ylab = "Deaths",
     main = sprintf("Observed vs Fitted Deaths (λ = %.1e)", lambda))
lines(day, mu_hat, lwd = 2, col = "black")
legend("topright", legend = c("Observed", "Fitted"), 
       pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2),
       col = c("red", "black"), bty = "n")

# (b) Infection curve with correct timing
plot(day_inf, f_hat, type = "l", lwd = 2, col = "blue",
     xlab = "Day of year 2020", ylab = "Infections (relative)",
     main = "Estimated infection curve f(t)")
abline(v = day[1], lty = 2)  # optional reference line

par(op)
## ----------------------------
## Task 4 — Select smoothing parameter λ by BIC
##
## Overview:
## Fitting a Poisson deconvolution model with a second-order difference 
## penalty across a set of candidate λ values.  
## For every fit compute:
##     BIC(λ) = -2·logLik(β̂_λ) + log(n)·EDF(λ)
## where EDF(λ) = trace(H_λ^{-1} H_0),
##       H_λ = X'WX + λS,  and  W = diag(y / μ̂_λ^2).
##
## Goal:
## Find the λ that minimizes the BIC criterion and keep track of
## all intermediate results for diagnostics and plotting.
##
## Design ideas:
## 1) “Fit–then–score”:  fit the model for each λ using the penalized
##    negative log-likelihood, then compute BIC from the unpenalized
##    likelihood and the EDF.
## 2) Hessian is computed in β-space (not γ-space)
## 3) Numerical stability and efficiency:By employing a warm-start strategy (utilizing the solution from the previous λ
## as the initial value for the next λ), the matrix W is not explicitly constructed, 
## with row scaling implemented to achieve X'WX.
## 4) Store both the best fit and the full BIC path for later use.
## ============================================================

## ---- Unpenalized Poisson log-likelihood ----
## Inputs:  mu (expected deaths), y (observed deaths)
## Output:  log-likelihood  Σ[y·logμ − μ]  ignoring constant Σ log(y!) because BIC is not require

loglik_poisson <- function(mu, y) sum(y * log(pmax(mu, 1e-12)) - mu) # pmax can prevent log(0)

## ---- Function factory binding λ ----
## Purpose: Create objective and gradient functions for a given λ.
## Inputs:  y, X, S,(same in task 2) lambda(different λ)
## Output:  a list with {obj, gr} functions usable by optim()
## Why:     Avoid repeatedly defining closures inside the main loop.
make_obj_gr <- function(y, X, S, lambda) {
  obj <- function(g) nll_gamma(g, y, X, S, lambda)
  gr  <- function(g) grad_gamma(g, y, X, S, lambda)
  list(obj = obj, gr = gr)
}

## grid for log-lambda
## Follow the rule: log λ ∈ [−13, −7], 50 points.
loglam_grid <- seq(-13, -7, length.out = 50)

bic_path <- data.frame(
  log_lambda = loglam_grid,
  lambda     = exp(loglam_grid),
  nll        = NA_real_,                  # The Penalized Negative Log-Likelihood under Optimal λ
  loglik     = NA_real_,                  # Unpenalized Log-Likelihood
  edf        = NA_real_,                  # The effective degrees of freedom of the model
  bic        = NA_real_,                  # BIC(λ) = -2·logLik(β̂_λ) + log(n)·EDF(λ
  conv       = NA_integer_                # Convergence flag of optimizer (0=converged)
)

## Warm start(using previous λ solution as initial value)
gamma_start <- gamma0
best <- list(bic = Inf)

## Main loop:fit, compute EDF AND BIC for each λ
for (k in seq_along(loglam_grid)) {
  lam <- bic_path$lambda[k]
  og  <- make_obj_gr(y, X, S, lam)        # creat obj/gr for current λ
  
  # The model is fitted through the minimization of the penalized negative log-likelihood function.
  # Obtain the value of γ̂_λ、β̂_λ、μ̂_λ
  fit_k <- optim(gamma_start, fn = og$obj, gr = og$gr, method = "BFGS",
                 control = list(maxit = 1000, reltol = 1e-8))
  
  gamma_hat <- fit_k$par
  beta_hat  <- exp(gamma_hat)
  mu_hat    <- drop(X %*% beta_hat)
  
  # compute EDF, H0 = X' W X, Hλ = H0 + λS, with W = diag(y / mu^2)
  w_vec <- as.numeric(y / pmax(mu_hat, 1e-12)^2)
  # compute crossprod(X, W X) without forming big diag: WX = X * w
  WX  <- X * w_vec
  H0  <- crossprod(X, WX)
  Hlam <- H0 + lam * S
  
  # EDF = trace(Hλ^{-1} H0)
  # Use solve(Hlam, H0) stably(than solve(Hlam) %*% H0)
  EDF <- sum(diag(solve(Hlam, H0)))
  
  # compute BIC
  ll  <- loglik_poisson(mu_hat, y)
  BIC <- -2 * ll + log(length(y)) * EDF
  
  # record metrics for diagnostics and plotting
  bic_path$nll[k]    <- fit_k$value
  bic_path$loglik[k] <- ll
  bic_path$edf[k]    <- EDF
  bic_path$bic[k]    <- BIC
  bic_path$conv[k]   <- fit_k$convergence
  
  # warm start next λ with current optimum
  gamma_start <- gamma_hat
  
  # keep best BIC solution
  if (BIC < best$bic) {
    best <- list(bic = BIC, lambda = lam, log_lambda = log(lam),
                 fit = fit_k, gamma = gamma_hat, beta = beta_hat,
                 mu = mu_hat, edf = EDF)
  }
}

# output best λ and full search path
best_lambda   <- best$lambda
best_results4 <- list(best = best, path = bic_path)

cat(sprintf("Task 4 — Best lambda = %.3e (logλ = %.3f), BIC = %.3f, EDF = %.2f\n",
            best_lambda, log(best_lambda), best$bic, best$edf))

## ============================================================
## Task 5 — Assess uncertainty of f(t) using nonparametric bootstrap
##
## Overview:
## Utilizing non-parametric Bootstrap (resampling) methodology to 
## evaluate the uncertainty of the infection curve f̂(t)
##
## Methodology:
## 1. The data consist of n daily death observations (y₁, …, yₙ).
##    In each iteration, Bootstrap employs a "with replacement" sampling method to draw n samples from these days.
## 2. This resampling is equivalent to reweighting the Poisson log-likelihood by
##    integer weights wᵢ = 0, 1, 2, … counting how many times each day is resampled.
## 3. For each replicate, re-fit the penalized Poisson model using the same λ,
##    obtain β̂*, compute f̂*(t) = X̃ β̂*, and store the result.
## 4. After B (200) replicates, summarize the bootstrap distribution of f̂*(t)
##    to form 95% confidence intervals for each time point.
##
## Design notes:
## - This implementation uses *reweighting*, so X, S and π(j) do not need
##   to be recomputed for each bootstrap sample.
## - Warm-starts (initializing γ with previous fit) speed up convergence.
## - The result is a list containing the point estimate f̂, the bootstrap replicates,
##   and the lower/upper 95% confidence limits.
## ============================================================

## ---- Weighted negative log-likelihood ----
## Inputs:
##   gamma  : log-coefficients γ ensuring β = exp(γ) > 0
##   y      : numeric vector of observed daily deaths
##   X      : n × K model matrix linking infections to deaths
##   S      : K × K penalty matrix enforcing smoothness
##   lambda : fixed smoothing parameter λ (chosen from Task 4)
##   w      : bootstrap weights (integer counts of how often each day is resampled)
##
## Output:
##   Penalized negative log-likelihood with reweighting
##
## Explanation:
##   Reweighting means we multiply each observation's contribution by wᵢ.
##   This avoids reconstructing resampled datasets and reusing X, S directly.
nll_gamma_w <- function(gamma, y, X, S, lambda, w) {
  beta <- exp(gamma)
  mu   <- drop(X %*% beta)
  mu   <- pmax(mu, 1e-12)
  pois <- sum(w * (mu - y * log(mu)))
  pen  <- 0.5 * lambda * drop(t(beta) %*% S %*% beta)
  pois + pen
}

## ---- Weighted gradient of penalized NLL ----
## Purpose:
##   Compute exact gradient (∇γ L) for the weighted case.
##
## Explanation:
##   The only difference from grad_gamma() is that each observation’s
##   Poisson term is multiplied by its bootstrap weight wᵢ.
grad_gamma_w <- function(gamma, y, X, S, lambda, w) {
  beta <- exp(gamma)
  mu   <- drop(X %*% beta)
  mu   <- pmax(mu, 1e-12)                   # prevent log(0)
  # X' (w * (1 - y/mu))  —— elementwise weight w on each data term
  g_beta_poiss <- crossprod(X, w * (1 - (y / mu)))
  g_beta_pen   <- lambda * (S %*% beta)
  g_beta <- drop(g_beta_poiss) + g_beta_pen
  beta * g_beta
}

# use best lambda from Task 4; start from its optimizer result if available
lambda_boot <- if (exists("best_lambda")) best_lambda else lambda
gamma_init  <- if (exists("best") && !is.null(best$gamma)) best$gamma else gamma0

B <- 200L                            # bootstrap replication count
n <- length(y)
nf <- nrow(X_tilde)
f_boot <- matrix(NA_real_, nrow = nf, ncol = B)       # store

# main loop: bootstrap fitting
for (b in 1:B) {
  # step 1: Generate non-parametric Bootstrap weights.
  # w_b[i]: The number of times day i got picked
  w_b <- tabulate(sample.int(n, size = n, replace = TRUE), nbins = n)
  
  # step 2: Definition of the Weighted Objective Function and Gradient Function
  obj_w <- function(g) nll_gamma_w(g, y, X, S, lambda_boot, w_b)
  gr_w  <- function(g) grad_gamma_w(g, y, X, S, lambda_boot, w_b)
  
  # step 3: The weighted fitting process is implemented utilizing the BFGS algorithm.
  fit_b <- optim(gamma_init, fn = obj_w, gr = gr_w, method = "BFGS",
                 control = list(maxit = 1000, reltol = 1e-8))
  
  # step 4: Estimated infection curve computation f̂*(t)
  beta_b <- exp(fit_b$par)
  f_b    <- drop(X_tilde %*% beta_b)
  f_boot[, b] <- f_b
  
  # step 5: warm start, The current results shall be utilized as the initial values for the subsequent iteration.
  gamma_init <- fit_b$par
  if (b %% 20 == 0) cat("Bootstrap replicate:", b, " / ", B, "\n")
}

# summarize: point estimate (from best) + 95% CI
beta_star <- if (exists("best") && !is.null(best$beta)) best$beta else exp(task3_result$fit$par)
f_hat_star <- drop(X_tilde %*% beta_star)

f_ci_lo <- apply(f_boot, 1, quantile, probs = 0.025, na.rm = TRUE)
f_ci_hi <- apply(f_boot, 1, quantile, probs = 0.975, na.rm = TRUE)

# store the result for task 6
task5_results <- list(
  lambda = lambda_boot,
  f_hat  = f_hat_star,
  f_boot = f_boot,
  f_ci   = cbind(lo = f_ci_lo, hi = f_ci_hi)
)

cat("Task 5 — Bootstrap completed: B =", B, "\n")


# Task 6

# Purpose:
#    Plot the daily deaths and the fitted death curve (under the selected λ), and show the estimated infection
#    curve f(t) together with its 95% confidence band.
#    
# Input:
#   Uses objects that are created in previous tasks:
#   - engcov$julian: day-of-year index for 2020
#   - y: observed daily deaths
#   - best: list containing the best-λ fit (including best$mu, best$lambda)
#   - task5_results: list with point estimate f_hat and bootstrap CI for f(t)
#
# Output:
#   A two-panel plot - Observed vs fitted daily deaths over day-of-year and 
#   estimated infection curve f(t) with 95% bootstrap confidence interval.

day <- engcov$julian
n   <- length(y)

# Fitted deaths under the λ from Task 4
mu_best <- best$mu

# Point estimate and 95% CI (from Task 5)
f_estimated <- task5_results$f_hat
f_lower    <- task5_results$f_ci[, "lo"]
f_upper    <- task5_results$f_ci[, "hi"]

# Time axis for f(t)
nf      <- length(f_estimated)
# Over from day[1] - 30 to day[n]
day_inf <- day[1] - 31 + seq_len(nf)

# Use a common x-axis range
x_lim <- range(c(day_inf, day))

# Plot in 2-row layout
op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# 1. Daily deaths and fitted deaths over day-of-year
plot(day, y, pch = 16, cex = 0.4, col = "red",
     xlab = "Day of year 2020",
     ylab = "Daily deaths",
     main = sprintf("Daily COVID-19 deaths in English hospitals (λ = %.1e)", best$lambda),
     xlim = x_lim)
lines(day, mu_best, lwd = 2)
legend("topright", legend = c("Observed", "Fitted"),
       pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2),
       col = c("red", "black"), bty = "n")

# 2. Infection curve f(t) with 95% confidence band
plot(day_inf, f_estimated, type = "n",
     xlab = "Day of year 2020",
     ylab = "Daily new infections (relative units)",
     main = "Estimated infection curve f(t) with 95% bootstrap CI",
     xlim = x_lim)

# 95% confidence band as a shaded polygon
polygon(
  x = c(day_inf, rev(day_inf)),
  y = c(f_lower, rev(f_upper)),
  border = NA,
  col = rgb(0.7, 0.7, 0.7, 0.5)
)

# Add the point estimate of f(t)
lines(day_inf, f_estimated, lwd = 2)

# Reference line that marks the first day with observed deaths
abline(v = day[1], lty = 2)

legend("topright", legend = c("f(t) estimate", "95% CI"),
       lty = c(1, NA), lwd = c(2, NA),
       pch = c(NA, 15),
       col = c("black", rgb(0.5, 0.5, 0.5, 0.7)),
       pt.cex = 1.5,
       bty = "n")

par(op)