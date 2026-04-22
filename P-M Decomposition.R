set.seed(1)

softmax <- function(eta) {
  # eta: M x N
  m <- apply(eta, 1, max)
  ex <- exp(eta - m)
  ex / rowSums(ex)
}

joint_loglik <- function(X, gamma, Z, B, Y) {
  # X: MxK (volume covs), gamma: K
  # Z: MxK2 (composition covs), B: K2xN
  # Y: MxN counts (y_{ij})
  n <- rowSums(Y)          # n_i = sum_{j=1}^N y_{ij}
  eta <- drop(X %*% gamma)
  lam <- exp(eta)
  
  P <- softmax(Z %*% B)    # MxN, p_{ij}
  
  ll <- sum(n * eta - lam) +
    sum(Y * log(P + 1e-300)) -
    sum(lgamma(Y + 1))      # const wrt params but included for completeness
  ll
}

fit_softmax_mle <- function(Z, Y, maxit = 200) {
  # maximize sum_{i=1}^M sum_{j=1}^N y_{ij} log softmax(ZB)_{ij}
  M <- nrow(Y); N <- ncol(Y); K2 <- ncol(Z)
  
  nll <- function(bvec) {
    B <- matrix(bvec, K2, N)
    eta <- Z %*% B
    logden <- apply(eta, 1, function(r) log(sum(exp(r - max(r)))) + max(r))
    logP <- eta - logden
    -sum(Y * logP)
  }
  
  b0 <- rep(0, K2 * N)
  opt <- optim(b0, nll, method = "BFGS", control = list(maxit = maxit))
  list(Bhat = matrix(opt$par, K2, N), opt = opt)
}

# ---------- simulate ----------
M <- 5000   # samples i=1,...,M
K <- 3
N <- 4      # categories j=1,...,N

X <- cbind(1, matrix(rnorm(M * (K - 1)), M, K - 1))  # intercept + 2 covs
gamma_true <- c(0.2, -0.4, 0.3)

eta <- drop(X %*% gamma_true)
lam <- exp(eta)
n <- rpois(M, lam)

K2 <- 2
Z <- X[, 1:K2, drop = FALSE]  # intercept + 1 cov (for composition)
B_true <- matrix(rnorm(K2 * N, sd = 0.5), K2, N)

P <- softmax(Z %*% B_true)
Y <- t(sapply(1:M, function(i) rmultinom(1, size = n[i], prob = P[i, ])))
stopifnot(all(rowSums(Y) == n))

# ---------- separability check ----------
B0 <- B_true
B1 <- B0 + matrix(rnorm(length(B0), sd = 0.1), nrow(B0), ncol(B0))

gamma0 <- gamma_true
gamma2 <- gamma_true + c(0.3, -0.2, 0.1)

d1 <- joint_loglik(X, gamma0, Z, B1, Y) - joint_loglik(X, gamma0, Z, B0, Y)
d2 <- joint_loglik(X, gamma2, Z, B1, Y) - joint_loglik(X, gamma2, Z, B0, Y)

cat("Delta ll from changing B at gamma0:", d1, "\n")
cat("Delta ll from changing B at gamma2:", d2, "\n")
cat("Difference (should be ~0):", d1 - d2, "\n\n")

# ---------- fit gamma from n only ----------
fit_pois <- glm(n ~ X[,2] + X[,3], family = poisson())
cat("gamma_true:", gamma_true, "\n")
cat("gamma_hat :", coef(fit_pois), "\n\n")

# ---------- fit B from Y only ----------
fit_soft <- fit_softmax_mle(Z, Y, maxit = 300)
cat("softmax optim converged:", fit_soft$opt$convergence == 0, "\n")
cat("B_true first row:", B_true[1, ], "\n")
cat("B_hat  first row:", fit_soft$Bhat[1, ], "\n")