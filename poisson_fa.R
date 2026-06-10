# =============================================================================
# Poisson FA
# =============================================================================

set.seed(8008308)

# =============================================================================
# 1. DATA GENERATION
# =============================================================================
generate_data <- function(N=800, Q=10, P=2, J=30, K=2,
                          sigma2_true=0.3, M_mean=150, seed=2026) { 
  # M_mean :The Mean of counts in the Poisson distribution
  set.seed(seed)
  mu_true <- rnorm(Q, 0, 0.5)
  mu_true <- mu_true - mean(mu_true)
  B_true <- matrix(0, Q, K)
  for (k in 1:K) {
    B_true[k, k] <- runif(1, 0.5, 1.2)
    if (k < Q) B_true[(k+1):Q, k] <- rnorm(Q - k, 0, 0.6)
  }
  phi_true <- matrix(0, P, Q)
  if (P > 1) phi_true[2:P, ] <- matrix(rnorm((P-1)*Q, 0, 0.3), P-1, Q)
  F_true <- matrix(rnorm(K * J), K, J)
  alpha_true <- matrix(0, Q, J)
  for (j in 1:J)
    alpha_true[, j] <- B_true %*% F_true[, j] + rnorm(Q, 0, sqrt(sigma2_true))
  group <- sample(1:J, N, replace = TRUE)
  X <- cbind(1, matrix(rnorm(N * (P-1)), N, P-1))
  eta <- matrix(0, N, Q)
  for (i in 1:N)
    eta[i, ] <- mu_true + X[i, ] %*% phi_true + alpha_true[, group[i]]
  M <- rpois(N, M_mean)
  xi <- log(M) - log(rowSums(exp(eta)))
  Y <- matrix(0, N, Q)
  for (i in 1:N) { p <- exp(eta[i,]+xi[i]); Y[i,] <- rmultinom(1,M[i],p/sum(p)) }
  list(Y=Y, X=X, group=group, M=M,
       mu_true=mu_true, B_true=B_true, sigma2_true=sigma2_true,
       phi_true=phi_true, F_true=F_true, alpha_true=alpha_true,
       N=N, Q=Q, P=P, J=J, K=K)
}

# =============================================================================
# 2. PME OFFSET
# =============================================================================
compute_xi <- function(M, eta_full) log(M) - log(rowSums(exp(eta_full)))

# =============================================================================
# 3. LAPLACE ON alpha_j (Q-dimensional)
# =============================================================================
laplace_alpha_j <- function(Y_j_mat, X_j, M_j, mu, phi, Sigma_inv,
                             max_iter=50, tol=1e-9) {
  N_j <- nrow(Y_j_mat); Q <- ncol(Y_j_mat)
  fixed <- matrix(rep(mu, each=N_j), N_j, Q) + X_j %*% phi

  log_post <- function(alpha) {
    eta <- sweep(fixed, 2, alpha, "+")
    lp_m <- apply(eta, 1, max)
    lse <- lp_m + log(rowSums(exp(eta - lp_m)))
    sum(Y_j_mat * eta) - sum(M_j * lse) -
      0.5 * as.numeric(t(alpha) %*% Sigma_inv %*% alpha)
  }

  grad_alpha <- function(alpha) {
    eta <- sweep(fixed, 2, alpha, "+")
    lp_m <- apply(eta, 1, max)
    exp_e <- exp(eta - lp_m)
    pi_mat <- exp_e / rowSums(exp_e)
    as.vector(colSums(Y_j_mat - sweep(pi_mat, 1, M_j, "*")) -
                Sigma_inv %*% alpha)
  }

  hess_alpha <- function(alpha) {
    eta <- sweep(fixed, 2, alpha, "+")
    lp_m <- apply(eta, 1, max)
    exp_e <- exp(eta - lp_m)
    pi_mat <- exp_e / rowSums(exp_e)
    M_pi <- sweep(pi_mat, 1, M_j, "*")
    wpi <- sweep(pi_mat, 1, sqrt(M_j), "*")
    -diag(colSums(M_pi), Q) + t(wpi) %*% wpi - Sigma_inv
  }

  alpha <- rep(0, Q); lp_cur <- log_post(alpha)
  for (iter in 1:max_iter) {
    g <- grad_alpha(alpha)
    H <- hess_alpha(alpha)
    # Fix 1: fallback uses -g (not +g) so alpha - step*(-g) = alpha + step*g (uphill)
    delta <- tryCatch(solve(H, g), error = function(e) -g * 0.01)
    step <- 1.0
    for (bt in 1:15) {
      fn <- alpha - step * delta
      if (log_post(fn) >= lp_cur - 1e-12) break
      step <- step * 0.5
    }
    alpha <- alpha - step * delta
    lp_cur <- log_post(alpha)
    if (max(abs(step * delta)) < tol) break
  }
  H_mode <- hess_alpha(alpha)
  S <- tryCatch(solve(-H_mode), error = function(e) diag(Q) * 0.1)
  S <- (S + t(S)) / 2
  list(alpha_hat = alpha, S_hat = S, lp_mode = lp_cur)
}

# =============================================================================
# 4. E-STEP: Laplace for all J groups
# =============================================================================
run_estep <- function(J, group, Y, X, M, mu, phi, Sigma_inv, Sigma) {
  Q <- length(mu)
  alpha_hat <- matrix(0, Q, J)
  S_hat <- vector("list", J)
  lp_total  <- 0
  for (j in 1:J) {
    idx <- which(group == j)
    if (!length(idx)) {
      S_hat[[j]] <- tryCatch(solve(Sigma_inv), error = function(e) diag(Q) * 0.5)
      next
    }
    lap <- laplace_alpha_j(Y[idx,,drop=FALSE], X[idx,,drop=FALSE],
                                     M[idx], mu, phi, Sigma_inv)
    alpha_hat[,j] <- lap$alpha_hat
    S_hat[[j]] <- lap$S_hat
    lp_total <- lp_total + lap$lp_mode
  }
  list(alpha_hat = alpha_hat, S_hat = S_hat, lp_total = lp_total)
}

# =============================================================================
# 5. RUBIN-THAYER (spherical D = sigma2 * I)
# =============================================================================
rubin_thayer_spherical <- function(Sigma_obs, K, B_init=NULL, sigma2_init=0.3,
                                    max_iter=500, tol=1e-10) {
  Q <- nrow(Sigma_obs)
  if (is.null(B_init)) {
    sv  <- svd(Sigma_obs, nu=K, nv=K)
    lam <- pmax(sv$d[1:K] - sigma2_init, 0.05)
    B <- sv$u %*% diag(sqrt(lam), K)
  } else { B <- B_init }
  sigma2 <- max(sigma2_init, 1e-6)
  for (iter in 1:max_iter) {
    B_old <- B; s_old <- sigma2
    Sigma <- B %*% t(B) + sigma2 * diag(Q)
    Si <- tryCatch(solve(Sigma), error = function(e) diag(Q) / sigma2)
    beta <- t(B) %*% Si
    Theta <- diag(K) - beta %*% B + beta %*% Sigma_obs %*% t(beta)
    B_new <- Sigma_obs %*% t(beta) %*% solve(Theta)
    sigma2_new <- max(sum(diag(Sigma_obs - B_new %*% beta %*% Sigma_obs)) / Q, 1e-6)
    B <- B_new; sigma2 <- sigma2_new
    if (max(abs(B - B_old)) < tol && abs(sigma2 - s_old) < tol) break
  }
  list(B = B, sigma2 = sigma2)
}

# =============================================================================
# 6. PLT CONSTRAINT
# =============================================================================
apply_PLT <- function(B) {
  K <- ncol(B)
  qr_obj <- qr(t(B[1:K, , drop=FALSE]))
  B_rot  <- B %*% qr.Q(qr_obj)
  for (k in 1:K) if (B_rot[k,k] < 0) B_rot[,k] <- -B_rot[,k]
  for (k in 1:K) if (k > 1) B_rot[1:(k-1), k] <- 0
  B_rot
}

# =============================================================================
# 7. MAIN EM
# =============================================================================
poisson_fa_em <- function(dat, K, max_iter=80, tol=1e-4, verbose=TRUE) {

  Y <- dat$Y 
  X <- dat$X 
  group <- dat$group 
  M <- dat$M
  N <- dat$N
  Q <- dat$Q
  P <- dat$P
  J <- dat$J

  avg_prop <- colMeans(Y / pmax(rowSums(Y), 1))
  mu       <- log(avg_prop + 1e-8); mu <- mu - mean(mu)
  phi      <- matrix(0, P, Q)

  # Group-level initialization using totals + pseudo-counts (avoids NaN)
  gm <- matrix(0, J, Q)
  for (j in 1:J) {
    idx <- which(group == j)
    if (length(idx) > 0) {
      gs     <- colSums(Y[idx, , drop=FALSE])
      gm[j,] <- log((gs + 1e-5) / (sum(gs) + Q * 1e-5))
    }
  }
  gm_c <- sweep(gm, 2, colMeans(gm))
  sv0 <- svd(t(gm_c), nu=K, nv=K)
  B <- sv0$u %*% diag(pmax(sv0$d[1:K] * 0.5, 0.1), K)
  B <- apply_PLT(B)
  sigma2 <- 0.5   # intentionally wrong init to test sigma2 update

  alpha_hat <- matrix(0, Q, J)
  logev_vec <- numeric(max_iter)

  for (em_iter in 1:max_iter) {

    # Current prior
    Sigma <- B %*% t(B) + sigma2 * diag(Q)
    Sigma_inv <- tryCatch(solve(Sigma), error = function(e) diag(Q) / sigma2)

    # ===== E-STEP =====
    es <- run_estep(J, group, Y, X, M, mu, phi, Sigma_inv, Sigma)
    alpha_hat <- es$alpha_hat
    S_hat <- es$S_hat

    # ===== M-STEP (a): mu and phi =====
    eta_full <- matrix(0, N, Q)
    for (i in 1:N)
      eta_full[i,] <- mu + X[i,] %*% phi + alpha_hat[, group[i]]
    xi_all <- compute_xi(M, eta_full)
    for (q in 1:Q) {
      fit <- tryCatch(
        glm.fit(x=X, y=Y[,q],
                offset = xi_all + alpha_hat[q, group],
                family = poisson(link="log"),
                control = glm.control(maxit=30, epsilon=1e-8)),
        error = function(e) NULL)
      if (!is.null(fit) && fit$converged) {
        phi[,q]  <- fit$coefficients
        mu[q] <- fit$coefficients[1]
        phi[1,q] <- 0
      }
    }

    # ===== M-STEP (b): B and sigma2 =====
    Sigma_obs <- matrix(0, Q, Q)
    for (j in 1:J)
      Sigma_obs <- Sigma_obs + outer(alpha_hat[,j], alpha_hat[,j]) + S_hat[[j]]
    Sigma_obs <- (Sigma_obs / J + t(Sigma_obs / J)) / 2

    rt <- rubin_thayer_spherical(Sigma_obs, K, B_init=B, sigma2_init=sigma2)
    B <- apply_PLT(rt$B)
    sigma2 <- rt$sigma2

    # ===== Laplace log-evidence convergence criterion =====
    ln_det_Sigma <- as.numeric(determinant(Sigma, logarithm=TRUE)$modulus)
    ln_det_S_total <- sum(sapply(S_hat, function(S)
      as.numeric(determinant(S, logarithm=TRUE)$modulus)))
    logev <- es$lp_total - 0.5 * J * ln_det_Sigma + 0.5 * ln_det_S_total
    logev_vec[em_iter] <- logev

    if (verbose)
      cat(sprintf("Iter %3d | LogEvidence=%12.2f | delta=%9.5f | sigma2=%.5f\n",
                  em_iter, logev,
                  if (em_iter > 1) logev - logev_vec[em_iter-1] else NA_real_,
                  sigma2))

    if (em_iter > 3 &&
        abs(logev_vec[em_iter] - logev_vec[em_iter-1]) < tol) {
      cat(sprintf("Converged at iter %d\n", em_iter)); break
    }
  }

  # Post-convergence factor scores
  V_f <- tryCatch(solve(diag(K) + t(B) %*% B / sigma2), error = function(e) diag(K))
  f_hat <- V_f %*% t(B) %*% alpha_hat / sigma2

  list(mu=mu, phi=phi, B=B, sigma2=sigma2,
       alpha_hat=alpha_hat, f_hat=f_hat, S_hat=S_hat,
       logev_vec=logev_vec[1:em_iter])
}

# =============================================================================
# 8. EVALUATION
# =============================================================================
subspace_dist <- function(A, B) {
  Qa <- qr.Q(qr(A)); Qb <- qr.Q(qr(B))
  sv <- svd(t(Qa) %*% Qb)$d
  sqrt(sum(sin(acos(pmin(sv, 1)))^2))
}

# =============================================================================
# 9. RUN
# =============================================================================
cat("=== Data ===\n")
dat <- generate_data(N=800, Q=50, P=2, J=30, K=2, sigma2_true=0.3)
cat(sprintf("N=%d Q=%d P=%d J=%d K=%d sigma2_true=%.2f\n",
            dat$N, dat$Q, dat$P, dat$J, dat$K, dat$sigma2_true))

cat("\n=== Final Laplace-EM (all 4 fixes applied) ===\n")
fit <- poisson_fa_em(dat, K=2, max_iter=60, tol=1e-4, verbose=TRUE)

cat("\n=== Results ===\n")
pg_cor <- sapply(1:dat$J, function(j)
  cor(fit$alpha_hat[,j], dat$alpha_true[,j]))
cat(sprintf("alpha | mean cor=%.4f | min=%.4f | n_neg=%d\n",
            mean(pg_cor), min(pg_cor), sum(pg_cor < 0)))
cat(sprintf("B dist| subspace=%.4f (worst=%.3f)\n",
            subspace_dist(fit$B, dat$B_true), sqrt(dat$K)))
cat(sprintf("sigma2| est=%.4f  true=%.4f\n", fit$sigma2, dat$sigma2_true))
if (dat$P > 1)
  cat(sprintf("phi | cor=%.4f\n",
              cor(as.vector(fit$phi[2:dat$P,]),
                  as.vector(dat$phi_true[2:dat$P,]))))

cat("\nLogEvidence trace (all values, should be monotone):\n")
print(round(fit$logev_vec, 2))

cat("\n=== Monotonicity check ===\n")
deltas <- diff(fit$logev_vec)
cat(sprintf("n_negative_deltas = %d / %d\n", sum(deltas < 0), length(deltas)))
cat(sprintf("max negative delta = %.6f\n", min(deltas)))

# =============================================================================
# F Recovery Check Uses Procrustes
# =============================================================================
# Find optimal orthogonal R minimizing ||f_true - R %*% f_hat||_F
sv <- svd(dat$F_true %*% t(fit$f_hat)) # SVD of f_true %*% f_hat'
R <- sv$u %*% t(sv$v) # optimal rotation (K x K, orthogonal)
f_aligned <- R %*% fit$f_hat # aligned estimate (K x J)

# Per-factor Pearson correlation after alignment
f_cors <- sapply(1:dat$K, function(k)
  cor(f_aligned[k, ], dat$F_true[k, ]))

cat("Factor score recovery (after Procrustes):\n")
for (k in 1:dat$K)
  cat(sprintf("  Factor %d: r = %.4f\n", k, f_cors[k]))

# Overall relative Frobenius error after alignment
f_rel_err <- norm(f_aligned - dat$F_true, "F") / norm(dat$F_true, "F")
cat(sprintf("Relative Frobenius error: %.4f\n", f_rel_err))

# Back to aligned values
rbind(
  f_true   = dat$F_true[, 1:5],
  f_aligned = f_aligned[, 1:5]
)

# =============================================================================
# Diagnostic Plots for Poisson Factor Analysis EM
# =============================================================================

library(ggplot2)
library(patchwork)

# =============================================================================
# PLOT 1: Laplace log-evidence convergence
# =============================================================================

df1 <- data.frame(
  iter  = seq_along(fit$logev_vec),
  logev = fit$logev_vec,
  delta = c(NA_real_, diff(fit$logev_vec))
)
# Flag any step where the evidence DECREASES by more than numerical noise
df1$violation <- !is.na(df1$delta) & df1$delta < -1e-6

p1 <- ggplot(df1, aes(x = iter, y = logev)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(aes(color = violation), size = 2.4, show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
  labs(
    title = "Laplace log-evidence over EM iterations",
    subtitle = sprintf(
      "%d monotonicity violation(s) \u2014 red dots indicate steps with decrease > 1e-6",
      sum(df1$violation, na.rm = TRUE)),
    x = "EM iteration",
    y = "Log-evidence"
  ) +
  theme_bw(base_size = 11)

# =============================================================================
# PLOT 3: Random effects recovery
# =============================================================================

# Per-group Pearson correlation between estimated and true alpha
pg_cor <- sapply(1:dat$J, function(j)
  cor(fit$alpha_hat[, j], dat$alpha_true[, j]))

# Build long-format dataframe: one row per (q, j) entry
# as.vector() on a Q x J matrix goes column-major, so group j occupies
# rows ((j-1)*Q + 1) : (j*Q), matching rep(pg_cor, each = dat$Q)
df3 <- data.frame(
  true  = as.vector(dat$alpha_true),
  est = as.vector(fit$alpha_hat),
  cor_j = rep(pg_cor, each = dat$Q)
)

p3 <- ggplot(df3, aes(x = true, y = est, color = cor_j)) +
  geom_point(alpha = 0.35, size = 0.7) +
  geom_abline(intercept = 0, slope = 1,
              color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_color_viridis_c(
    name = "Group\ncorr.",
    limits = c(0, 1),
    option = "plasma"
  ) +
  # Annotation: summary statistics in top-left corner
  annotate("text",
           x = -Inf, y = Inf, hjust = -0.1, vjust = 1.6,
           label = sprintf("mean r = %.3f  |  min r = %.3f",
                           mean(pg_cor), min(pg_cor)),
           size = 3.2) +
  labs(
    title = "Random effects recovery",
    x = expression(paste("True  ", alpha[jq]^"*")),
    y = expression(paste("Estimated  ", hat(alpha)[jq]))
  ) +
  theme_bw(base_size = 11)

# =============================================================================
# PLOT 4: Multi-seed robustness (20 independent replications)
# =============================================================================

cat("Running 20 independent replications for robustness check.\n")

df4 <- do.call(rbind, lapply(1:20, function(s) {
  
  # Generate a fresh dataset under the same DGP with a different seed
  dat_s <- generate_data(N = 800, Q = 50, P = 2, J = 30,
                         K = 2, sigma2_true = 0.3, seed = s)
  
  # Run EM; suppress per-iteration console output
  invisible(capture.output(
    fit_s <- poisson_fa_em(dat_s, K = 2, max_iter = 60,
                           tol = 1e-4, verbose = FALSE)
  ))
  
  data.frame(
    seed = s,
    # Subspace distance between estimated and true factor loading column spaces
    sub_dist = subspace_dist(fit_s$B, dat_s$B_true),
    # Estimated idiosyncratic variance
    sigma2 = fit_s$sigma2,
    # Mean Pearson correlation between alpha_hat_j and alpha_true_j across groups
    alpha_cor = mean(sapply(1:dat_s$J, function(j)
      cor(fit_s$alpha_hat[, j], dat_s$alpha_true[, j])))
  )
}))

# --- Panel A: B subspace distance ---
# sqrt(K) is the theoretical worst-case value (orthogonal subspaces)
p4a <- ggplot(df4, aes(x = "", y = sub_dist)) +
  geom_boxplot(width = 0.45, fill = "#AEC6CF", outlier.shape = NA) +
  geom_jitter(width = 0.07, alpha = 0.7, size = 1.8) +
  geom_hline(yintercept = sqrt(dat$K), color = "red",
             linetype = "dashed", linewidth = 0.7) +
  labs(
    title = "B subspace distance",
    x = NULL,
    y = expression(d[sub](hat(B), B^"*")),
    caption = sprintf("Dashed: worst-case bound sqrt(K) = %.3f", sqrt(dat$K))
  ) +
  theme_bw(base_size = 11) +
  theme(axis.ticks.x = element_blank())

# --- Panel B: sigma2 estimate ---
p4b <- ggplot(df4, aes(x = "", y = sigma2)) +
  geom_boxplot(width = 0.45, fill = "#FADADD", outlier.shape = NA) +
  geom_jitter(width = 0.07, alpha = 0.7, size = 1.8) +
  geom_hline(yintercept = dat$sigma2_true, color = "red",
             linetype = "dashed", linewidth = 0.7) +
  labs(
    title = expression(hat(sigma)^2 ~ "estimate"),
    x = NULL,
    y = expression(hat(sigma)^2),
    caption = sprintf("Dashed: true value %.1f", dat$sigma2_true)
  ) +
  theme_bw(base_size = 11) +
  theme(axis.ticks.x = element_blank())

# --- Panel C: mean alpha correlation ---
p4c <- ggplot(df4, aes(x = "", y = alpha_cor)) +
  geom_boxplot(width = 0.45, fill = "#B5EAD7", outlier.shape = NA) +
  geom_jitter(width = 0.07, alpha = 0.7, size = 1.8) +
  geom_hline(yintercept = 1.0, color = "red",
             linetype = "dashed", linewidth = 0.7) +
  ylim(NA, 1.05) +
  labs(
    title = expression(bar(r)(hat(alpha), alpha^"*")),
    x = NULL,
    y = "Mean group correlation",
    caption = "Dashed: perfect recovery"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.ticks.x = element_blank())

p4 <- (p4a | p4b | p4c) +
  plot_annotation(
    title = "Algorithm robustness across 20 independent replications",
    subtitle = sprintf("DGP fixed: N=%d, Q=%d, J=%d, K=%d, sigma2_true=%.1f",
                       dat$N, dat$Q, dat$J, dat$K, dat$sigma2_true)
  )

print(p1)
print(p3)
print(p4)

# =============================================================================
# PLOT 5: Factor score recovery after Procrustes alignment
# =============================================================================
df_f <- data.frame(
  true = c(dat$F_true[1, ],  dat$F_true[2, ]),
  aligned = c(f_aligned[1, ],   f_aligned[2, ]),
  factor = rep(paste0("Factor ", 1:dat$K), each = dat$J)
)

# Correlation labels positioned in top-left of each facet
cor_labels <- data.frame(
  factor = paste0("Factor ", 1:dat$K),
  label = sprintf("r = %.4f", f_cors)
)

p5 <- ggplot(df_f, aes(x = true, y = aligned)) +
  geom_point(alpha = 0.65, size = 2.2, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "black", linewidth = 0.6) +
  # Per-facet correlation annotation
  geom_text(data = cor_labels,
            aes(label = label),
            x = -Inf, y = Inf,
            hjust = -0.15, vjust = 1.8,
            size = 3.5, inherit.aes = FALSE) +
  facet_wrap(~factor, scales = "free") +
  labs(
    title = "Factor score recovery after Procrustes alignment",
    # subtitle = sprintf("Relative Frobenius error = %.4f  (reflects posterior shrinkage, not directional error)",
    #                    f_rel_err),
    x = expression(paste("True  ", f[jk]^"*")),
    y = expression(paste("Estimated  ", hat(f)[jk]))
  ) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

print(p5)
