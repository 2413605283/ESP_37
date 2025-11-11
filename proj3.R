library(splines)
evaluate_matrices <- function(K, n = 80) {
  # time grid
  t <- 1:n
  
  # 1. Build B-spline basis
  knots <- seq(1, n, length.out = K + 4)
  X_tilde <- splineDesign(knots, t, outer.ok = TRUE)
  
  # 2. Build infection-to-death pmf π(j)
  d <- 1:80
  edur <- 3.151; sdur <- 0.469
  pd <- dlnorm(d, edur, sdur)
  pd <- pd / sum(pd)
  
  # 3. Build convolution matrix X
  L <- length(pd)
  X <- matrix(0, n, K)
  for (j in seq_len(L)) {
    if (j < n)
      X[(j+1):n, ] <- X[(j+1):n, ] + pd[j] * X_tilde[1:(n-j), ]
  }
  
  # 4. Build penalty matrix S
  S <- crossprod(diff(diag(K), diff = 2))
  
  # 5. Return results
  list(X_tilde = X_tilde, X = X, S = S, pi = pd)
}

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
#   plus a smoothness penalty term that discourages overly wiggly infection curves f(t).
#
#   L(γ) = Σ_i [ μ_i − y_i log(μ_i) ] + (λ / 2) βᵀ S β
#   where  β = exp(γ), μ = X β
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
#   (4) Compute the smoothness penalty (λ / 2) βᵀSβ.
#   (5) Return the total NLL = likelihood term + penalty term.
#
nll_gamma <- function(gamma, y, X, S, lambda) {
  
  beta <- exp(gamma)                          # (1) transform parameters to keep β > 0
  mu   <- drop(X %*% beta)                    # (2) expected deaths from convolution model
  mu   <- pmax(mu, 1e-12)                     # avoid log(0) or division by zero
  pois <- sum(mu - y * log(mu))               # (3) Poisson negative log-likelihood part
  pen  <- 0.5 * lambda * drop(t(beta) %*% S %*% beta)  # (4) smoothness penalty
  pois + pen                                  # (5) return total penalized NLL
}


# Function: grad_gamma()
#
# Purpose:
#   Compute the exact gradient (first derivative vector) of the penalized
#   negative log-likelihood with respect to γ (the log-parameters).
#   This gradient is used by the BFGS optimizer for efficient model fitting.
#
#   ∇_γ L = diag(β) [ Xᵀ(1 − y / μ) + λ S β ]
#   where β = exp(γ), μ = X β
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
#       ∇_β L = Xᵀ(1 − y / μ) + λ S β
#       – the first term measures model–data mismatch (Poisson part),
#         and the second term is the smoothness penalty derivative.
#   (3) Apply chain rule: ∇_γ L = diag(β) ∇_β L, since dβ/dγ = β.
#   (4) Return the resulting vector.
#
grad_gamma <- function(gamma, y, X, S, lambda) {
  
  beta <- exp(gamma)                          # (1) parameter transformation
  mu   <- drop(X %*% beta)                    # predicted deaths
  mu   <- pmax(mu, 1e-12)                     # avoid division by zero
  
  g_beta_poiss <- crossprod(X, 1 - (y / mu))  # (2a) Poisson part: Xᵀ(1 − y/μ)
  g_beta_pen   <- lambda * (S %*% beta)       # (2b) penalty part: λ S β
  g_beta <- drop(g_beta_poiss) + g_beta_pen   # combine both parts
  
  g_gamma <- beta * g_beta                    # (3) apply chain rule: diag(β) × g_beta
  g_gamma                                    # (4) return gradient vector
}


# Test gradient (your concise version)
gamma0 <- rep(log(1e-3), ncol(S))
nll0 <- nll_gamma(gamma0, y, X, S, lambda)
g_exact <- grad_gamma(gamma0, y, X, S, lambda)
eps <- 1e-5
fd <- numeric(ncol(S))
for (i in 1:ncol(S)) {
  g1 <- gamma0; g1[i] <- g1[i] + eps
  g2 <- gamma0; g2[i] <- g2[i] - eps
  fd[i] <- (nll_gamma(g1, y, X, S, lambda) - nll_gamma(g2, y, X, S, lambda)) / (2 * eps)
}
# compute absolute error between finite-diff and analytic gradients
abs_err <- abs(fd[1:10] - g_exact[1:10])
table <- cbind("finite diff" = fd[1:10],
               "grad_exact" = g_exact[1:10],
               "abs_err" = abs_err)
print(table)


# simple diagnostic
if (max(abs_err) < 1) {
  cat("Gradient check: PASS \n")
} else {
  cat("Gradient check: FAIL (max abs error =", max(abs_err), ")\n")
}

## Task 3 (Sanity-check fit)
## ---------------------------------------

# Objective & gradient (from Task 2)
obj <- function(g) nll_gamma(g, y, X, S, lambda)
gr  <- function(g) grad_gamma(g, y, X, S, lambda)

# ---- Fit with BFGS ----
fit <- optim(gamma0, fn = obj, gr = gr, method = "BFGS")

# Extract estimates
beta_hat <- exp(fit$par)
mu_hat   <- drop(X %*% beta_hat)           # fitted deaths
f_hat    <- drop(X_tilde %*% beta_hat)      # estimated infection curve (relative units)

# Basic diagnostics
cat("Convergence:", fit$convergence, "(0 is good)\n")
cat("Final nll:", fit$value, "\n")
cat("||grad||_inf:", max(abs(gr(fit$par))), "\n")

# ---- Plots ----
op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# (a) Observed vs Fitted deaths
plot(y, type = "h", lwd = 2, xlab = "Day", ylab = "Deaths",
     main = sprintf("Observed vs Fitted Deaths (λ = %.1e)", lambda))
lines(mu_hat, lwd = 2)
legend("topleft", legend = c("Observed", "Fitted"), lty = 1, lwd = 2, bty = "n")

# (b) Estimated infection curve f(t)
plot(f_hat, type = "l", lwd = 2, xlab = "Day", ylab = "Infections (relative)",
     main = "Estimated f(t)")

# ---- Return a convenient list (optional) ----
task3_result <- list(
  fit = fit,
  beta = beta_hat,
  mu = mu_hat,
  f = f_hat,
  X = X, X_tilde = X_tilde, S = S,
  lambda = lambda
)
