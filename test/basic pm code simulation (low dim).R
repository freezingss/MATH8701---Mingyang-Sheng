# ==============================================================================
# Corrected Factor Analyzer: PME-EM Architecture (Matrix Dimension Fix Applied)
# ==============================================================================
set.seed(2026)

# ------------------------------------------------------------------------------
# 1. Simulation Dimensions and True Factor Structure
# ------------------------------------------------------------------------------
J <- 25         # Number of groups
N_j <- 50       # Number of individuals per group
Q <- 50         # Number of categories
K <- 2          # Latent factor dimension
Y_total <- 250  # Total multinomial trials per individual

# True idiosyncratic error variance (D = sigma2 * I)
true_sigma2 <- 0.25

cat("--- Step 1: Constructing true parameters satisfying upper-triangular constraint and D ---\n")

true_mu <- rnorm(Q, mean = 0, sd = 0.5)

# Factor loading matrix B (Q x K) satisfying upper-triangular constraint
# i.e., B[q, k] = 0 for k < q; diagonal entries strictly positive
true_B <- matrix(rnorm(Q * K, mean = 0, sd = 0.6), nrow = Q, ncol = K)
for (q in 1:K) {
  if (q > 1) {
    true_B[q, 1:(q-1)] <- 0   # Zero out lower-triangular entries for identifiability
  }
  true_B[q, q] <- abs(true_B[q, q]) + 0.4  # Enforce strictly positive diagonal
}

# Group-level latent factors F (J x K), f_j ~ N(0, I_K)
true_F <- matrix(rnorm(J * K, mean = 0, sd = 1), nrow = J, ncol = K)

# Idiosyncratic residuals E — NOTE: obs-level (J x N_j x Q) array
# This does not match the group-level epsilon_j in the draft; kept as-is for this version
true_E <- array(rnorm(J * N_j * Q, mean = 0, sd = sqrt(true_sigma2)), dim = c(J, N_j, Q))

# Generate observed count tensor Y
Y <- array(0, dim = c(J, N_j, Q))
for (j in 1:J) {
  for (i in 1:N_j) {
    eta  <- true_mu + true_B %*% true_F[j, ] + true_E[j, i, ]
    prob <- exp(eta) / sum(exp(eta))
    Y[j, i, ] <- rmultinom(1, size = Y_total, prob = prob)
  }
}
Y_plus <- apply(Y, c(1, 2), sum)  # Individual-level multinomial totals y_{ij+}

# ------------------------------------------------------------------------------
# 2. EM Algorithm Initialization
# ------------------------------------------------------------------------------
cat("--- Step 2: Initializing parameters ---\n")

# ==============================================================================
# Initialize mu via marginal MLE (Chiquet 2019):
# Compute per-category marginal log-probabilities under the marginal Poisson model
# ==============================================================================
total_margin_counts <- apply(Y, 3, sum)  # Total count per category q across all (j, i)
marginal_probs      <- total_margin_counts / sum(total_margin_counts)  # Marginal probability p_q

# Analytic MLE under log-linear model with offset
est_mu <- log(marginal_probs) + log(Q)
est_mu <- est_mu - mean(est_mu)  # Center to remove translation indeterminacy

est_F      <- matrix(0, nrow = J, ncol = K)       # Initialize latent factors at zero
est_E      <- array(0, dim = c(J, N_j, Q))        # Initialize idiosyncratic residuals at zero
est_sigma2 <- 0.5                                  # Initial guess for residual variance D
est_Vinv   <- array(0, dim = c(K, K, J))          # Posterior precision matrices (per group)

# Initialize B with upper-triangular zero constraint and positive diagonal
est_B <- matrix(0.1, nrow = Q, ncol = K)
for (q in 1:K) {
  if (q > 1) est_B[q, 1:(q-1)] <- 0
  est_B[q, q] <- 0.5
}

# Initialize PME profile parameters xi_{ij} = log(y_{ij+}) - log(Q)
est_xi <- log(Y_plus) - log(Q)

max_iter <- 15
cat("Starting EM iterations with upper-triangular B constraint and D estimation...\n\n")

# ------------------------------------------------------------------------------
# 3. PME-EM Main Loop
# ------------------------------------------------------------------------------
for (iter in 1:max_iter) {
  
  # ============================================================
  # E-STEP: Alternating posterior mode estimation for (f_j, e_{ij})
  # via Block Coordinate Descent
  # ============================================================
  for (j in 1:J) {
    fj <- est_F[j, ]
    
    for (alt in 1:3) {
      
      # --- Block 1: Fix f_j, update idiosyncratic residuals e_{ijq} (scalar Newton) ---
      for (i in 1:N_j) {
        for (q in 1:Q) {
          for (inner_e in 1:2) {
            lambda_ijq     <- exp(est_xi[j, i] + est_mu[q] + sum(est_B[q, ] * fj) + est_E[j, i, q])
            grad_e         <- Y[j, i, q] - lambda_ijq - est_E[j, i, q] / est_sigma2
            hess_e         <- -lambda_ijq - 1 / est_sigma2
            est_E[j, i, q] <- est_E[j, i, q] - grad_e / hess_e
          }
        }
      }
      
      # --- Block 2: Fix e_{ij}, update group latent factor f_j (Newton in R^K) ---
      for (inner_f in 1:2) {
        grad_f <- rep(0, K)
        Hess_f <- matrix(0, K, K)
        for (i in 1:N_j) {
          lambda_j <- as.vector(exp(est_xi[j, i] + est_mu + est_B %*% fj + est_E[j, i, ]))
          grad_f   <- grad_f + as.vector(t(est_B) %*% (Y[j, i, ] - lambda_j))
          Hess_f   <- Hess_f + t(est_B) %*% (lambda_j * est_B)
        }
        grad_f <- grad_f - fj          # Add prior gradient: -f_j (standard normal prior)
        Hess_f <- Hess_f + diag(K)     # Add prior Hessian: I_K
        fj     <- fj + solve(Hess_f, grad_f)
      }
    }
    est_F[j, ] <- fj
    
    # Compute group-level posterior precision matrix (Laplace approximation)
    # H_j = I_K + sum_i B' diag(w_{ij}) B, where w_{ijq} = lambda_{ijq} / (1 + sigma2 * lambda_{ijq})
    H_j <- diag(K)
    for (i in 1:N_j) {
      lambda_j <- as.vector(exp(est_xi[j, i] + est_mu + est_B %*% fj + est_E[j, i, ]))
      weight   <- lambda_j / (1 + est_sigma2 * lambda_j)
      H_j      <- H_j + t(est_B) %*% (weight * est_B)
    }
    est_Vinv[, , j] <- solve(H_j)  # Store posterior covariance V_j = H_j^{-1}
  }
  
  # Compute variance correction tensors S and var_eps for the Q-function
  # S[j,i,q]:       total conditional variance of eta_{ijq} under posterior
  # var_eps[j,i,q]: conditional variance of epsilon_{jq} under posterior
  var_eps <- array(0, dim = c(J, N_j, Q))
  S       <- array(0, dim = c(J, N_j, Q))
  for (j in 1:J) {
    fj     <- est_F[j, ]
    Vinv_j <- est_Vinv[, , j]
    for (i in 1:N_j) {
      lambda_ji <- as.vector(exp(est_xi[j, i] + est_mu + est_B %*% fj + est_E[j, i, ]))
      for (q in 1:Q) {
        l_val <- lambda_ji[q]
        denom <- 1 + est_sigma2 * l_val
        B_V_B <- as.numeric(t(est_B[q, ]) %*% Vinv_j %*% est_B[q, ])
        
        S[j, i, q]       <- est_sigma2 / denom + (1 / denom)^2 * B_V_B
        var_eps[j, i, q] <- est_sigma2 / denom + (est_sigma2 * l_val / denom)^2 * B_V_B
      }
    }
  }
  
  # ============================================================
  # M-STEP: Category-by-category Newton update for (mu_q, B_q)
  # Upper-triangular constraint enforced: B[q, k] = 0 for k < q
  # ============================================================
  for (q in 1:Q) {
    # Free (estimable) columns of B for row q under upper-triangular constraint
    free_B_idx <- if (q <= K) q:K else 1:K
    n_free_B   <- length(free_B_idx)
    theta_q    <- c(est_mu[q], est_B[q, free_B_idx])
    dim_q      <- 1 + n_free_B
    
    for (inner_m in 1:3) {
      mu_c   <- theta_q[1]
      B_full <- rep(0, K)
      if (n_free_B > 0) B_full[free_B_idx] <- theta_q[2:dim_q]
      
      g_q <- rep(0, dim_q)
      H_q <- matrix(0, dim_q, dim_q)
      
      for (j in 1:J) {
        fj <- est_F[j, ]
        for (i in 1:N_j) {
          # Jensen-corrected intensity: tilde_lambda = exp(xi + mu_q + B_q'f_j + e_{ijq} + 0.5*S)
          tilde_lambda <- exp(est_xi[j, i] + mu_c + sum(B_full * fj) + est_E[j, i, q] + 0.5 * S[j, i, q])
          
          g_mu <- Y[j, i, q] - tilde_lambda
          g_B  <- (Y[j, i, q] - tilde_lambda) * fj[free_B_idx]
          g_q  <- g_q + c(g_mu, g_B)
          
          H_q[1, 1]                         <- H_q[1, 1] - tilde_lambda
          H_q[1, 2:dim_q]                   <- H_q[1, 2:dim_q] - tilde_lambda * fj[free_B_idx]
          H_q[2:dim_q, 1]                   <- H_q[2:dim_q, 1] - tilde_lambda * fj[free_B_idx]
          H_q[2:dim_q, 2:dim_q]             <- H_q[2:dim_q, 2:dim_q] - tilde_lambda * outer(fj[free_B_idx], fj[free_B_idx])
        }
      }
      theta_q <- theta_q - solve(H_q, g_q)
    }
    est_mu[q] <- theta_q[1]
    if (n_free_B > 0) est_B[q, free_B_idx] <- theta_q[2:dim_q]
  }
  
  # ============================================================
  # D UPDATE: Closed-form MLE for residual variance sigma2
  # sigma2 = mean_q mean_j mean_i [ e_{ijq}^2 + Var(e_{ijq} | y) ]
  # ============================================================
  est_sigma2 <- mean(est_E^2 + var_eps)
  
  # ============================================================
  # SIGN ALIGNMENT: Anchor diagonal of B to be positive
  # Flip column k of B and F simultaneously if B[k,k] < 0
  # ============================================================
  for (k in 1:K) {
    if (est_B[k, k] < 0) {
      est_B[, k] <- -est_B[, k]
      est_F[, k] <- -est_F[, k]
    }
  }
  
  # ============================================================
  # PME PROFILING UPDATE: Refresh delta_{ij} (absorbs log-sum-exp coupling)
  # delta_{ij} = log(y_{ij+}) - log(sum_q exp(eta_{ijq}))
  # Jensen correction 0.5*S included in sum_exp to match Q-function expectation
  # ============================================================
  for (j in 1:J) {
    fj <- est_F[j, ]
    for (i in 1:N_j) {
      sum_exp <- 0
      for (q in 1:Q) {
        sum_exp <- sum_exp + exp(est_mu[q] + sum(est_B[q, ] * fj) + est_E[j, i, q] + 0.5 * S[j, i, q])
      }
      est_xi[j, i] <- log(Y_plus[j, i]) - log(sum_exp)
    }
  }
  
  # ============================================================
  # Monitoring: Pearson correlation and Frobenius norm distance
  # ============================================================
  corr_mu  <- cor(true_mu, est_mu)
  corr_B   <- cor(as.vector(true_B), as.vector(est_B))
  corr_F   <- cor(as.vector(true_F), as.vector(est_F))
  
  fnorm_mu <- sqrt(sum((true_mu - est_mu)^2))
  fnorm_B  <- sqrt(sum((true_B  - est_B)^2))
  fnorm_F  <- sqrt(sum((true_F  - est_F)^2))
  
  cat(sprintf("[Iter %02d]\n", iter))
  cat(sprintf("  Pearson corr  -> mu: %.4f | B: %.4f | F: %.4f\n", corr_mu, corr_B, corr_F))
  cat(sprintf("  Frobenius     -> mu: %.4f | B: %.4f | F: %.4f\n", fnorm_mu, fnorm_B, fnorm_F))
  cat(sprintf("  Residual var D -> true: %.2f | estimated: %.4f\n", true_sigma2, est_sigma2))
  cat("--------------------------------------------------------------------------\n")
}

cat("\n--- Done: PME-EM factor model simulation complete ---\n")