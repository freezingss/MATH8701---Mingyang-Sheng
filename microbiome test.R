# =============================================================================
# REAL DATASET: dietswap (microbiome package)
# =============================================================================

# --- Package setup ---
library(microbiome)
library(ggplot2)
library(patchwork)

# =============================================================================
# 1. DATA PREPARATION
# =============================================================================

data(dietswap)

# OTU count matrix
otu <- otu_table(dietswap)
Y_raw <- if (taxa_are_rows(otu)) t(as.matrix(otu)) else as.matrix(otu)

# Filter: keep taxa present in >= 5 samples to removes very rare taxa
keep_taxa <- colSums(Y_raw > 0) >= 5
Y_mat <- Y_raw[, keep_taxa]

# Sample metadata
meta <- data.frame(sample_data(dietswap))

# Covariate: nationality (African = 1, Finnish = 0)
# Binary indicator; no scaling needed for 0/1
X_mat <- cbind(1, as.numeric(meta$nationality == "AFR"))

# Group = subject ID — natural repeated-measures grouping
group_vec <- as.integer(as.factor(as.character(meta$subject)))

M_vec <- rowSums(Y_mat)

real_data <- list(
  Y = Y_mat,
  X = X_mat,
  group = group_vec,
  M = M_vec,
  N = nrow(Y_mat),
  Q = ncol(Y_mat),
  P = ncol(X_mat),
  J = length(unique(group_vec))
)

cat("=== Real Dataset Summary (dietswap) ===\n")
cat(sprintf("N = %d | Q = %d | P = %d | J = %d | mean N_j = %.1f\n",
            real_data$N, real_data$Q, real_data$P, real_data$J,
            real_data$N / real_data$J))
cat(sprintf("Total read count range: [%d, %d]\n", min(M_vec), max(M_vec)))

# =============================================================================
# 2. FIT MODEL
# =============================================================================

cat("\n=== Fitting Poisson Factor Analysis EM ===\n")
# fit_real <- poisson_fa_em(real_data, K = 2, max_iter = 120, tol = 1e-4, verbose = TRUE)
fit_real <- poisson_fa_em(real_data, K=3, max_iter=120, tol=1e-4, verbose=TRUE)
# fit_real <- poisson_fa_em(real_data, K=4, max_iter=120, tol=1e-4, verbose=TRUE)

# =============================================================================
# 3. CONVERGENCE DIAGNOSTICS
# =============================================================================

cat("\n=== Convergence Diagnostics ===\n")
deltas <- diff(fit_real$logev_vec)
cat(sprintf("Total EM iterations: %d\n",  length(fit_real$logev_vec)))
cat(sprintf("Monotonicity violations (>1e-6):  %d\n",  sum(deltas < -1e-6)))
cat(sprintf("Final log-evidence: %.4f\n", tail(fit_real$logev_vec, 1)))

df_ev <- data.frame(
  Iteration   = seq_along(fit_real$logev_vec),
  LogEvidence = fit_real$logev_vec
)
p_conv <- ggplot(df_ev, aes(x = Iteration, y = LogEvidence)) +
  geom_line(color = "darkblue", linewidth = 0.8) +
  geom_point(color = "darkblue", size = 1.8) +
  labs(title = "Laplace log-evidence over EM iterations (dietswap)",
       x = "EM iteration", y = "Log-evidence") +
  theme_bw(base_size = 11)
print(p_conv)

# =============================================================================
# 4. GOODNESS-OF-FIT: RECONSTRUCT PREDICTED COUNTS
# =============================================================================

X <- real_data$X
Y_obs <- real_data$Y
group <- real_data$group
M <- real_data$M

# (b) mu broadcast to N x Q
mu_mat <- matrix(fit_real$mu, nrow = nrow(X), ncol = length(fit_real$mu), byrow = TRUE)

# (a) correct linear predictor: mu + X*phi + alpha_hat[group]
# phi[1, ] == 0 so X*phi contributes only the non-intercept covariates
Eta <- mu_mat + X %*% fit_real$phi + t(fit_real$alpha_hat)[group, ]

# (c) stable softmax
lp_m <- apply(Eta, 1, max)
Pi_pred <- exp(Eta - lp_m) / rowSums(exp(Eta - lp_m))
Y_pred<- M * Pi_pred

# =============================================================================
# 5. GOODNESS-OF-FIT METRICS
# =============================================================================

cat("\n=== Goodness-of-Fit Metrics ===\n")

overall_r2 <- cor(as.vector(Y_obs), as.vector(Y_pred))^2
cat(sprintf("Overall matrix pseudo-R2: %.4f\n", overall_r2))

taxa_names <- colnames(Y_obs)
taxa_r2 <- sapply(seq_len(ncol(Y_obs)), function(q)
  cor(Y_obs[, q], Y_pred[, q])^2)
names(taxa_r2) <- taxa_names

cat("\nTop 5 best-fitted taxa (R2):\n")
print(round(sort(taxa_r2, decreasing = TRUE)[1:5], 4))
cat("\nTop 5 worst-fitted taxa (R2):\n")
print(round(sort(taxa_r2, decreasing = FALSE)[1:5], 4))

# =============================================================================
# 6. GOODNESS-OF-FIT PLOTS
# =============================================================================

df_gof <- data.frame(
  Observed  = as.vector(Y_obs),
  Predicted = as.vector(Y_pred),
  Taxa = rep(taxa_names, each = nrow(Y_obs))
)

# --- Plot A: overall scatter ---
p_gof_all <- ggplot(df_gof, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.3, color = "steelblue", size = 1.2) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(title = "Overall posterior predictive check (dietswap)",
       subtitle = sprintf("Overall pseudo-R\u00b2 = %.4f", overall_r2),
       x = "Observed counts", y = "Predicted counts") +
  theme_bw(base_size = 11)
print(p_gof_all)

# --- Plot B: faceted by top 6 most abundant taxa ---
top6 <- names(sort(colSums(Y_obs), decreasing = TRUE)[1:6])
p_gof_taxa <- ggplot(subset(df_gof, Taxa %in% top6),
                     aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.7, color = "chocolate", size = 2) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "red") +
  facet_wrap(~Taxa, scales = "free") +
  labs(title    = "Predictive check: top 6 most abundant taxa",
       subtitle = "Free scales; dashed line = perfect prediction",
       x = "Observed counts", y = "Predicted counts") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 7))
print(p_gof_taxa)

# --- Plot C: per-taxa R2 ranked bar chart ---
df_r2 <- data.frame(
  Taxa = names(taxa_r2),
  R2 = as.numeric(taxa_r2)
)
df_r2 <- df_r2[order(df_r2$R2, decreasing = TRUE), ]
df_r2$Taxa <- factor(df_r2$Taxa, levels = df_r2$Taxa)

p_r2_bar <- ggplot(df_r2, aes(x = Taxa, y = R2)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_hline(yintercept = mean(df_r2$R2), color = "red",
             linetype = "dashed", linewidth = 0.7) +
  labs(title = "Per-taxa pseudo-R\u00b2 (ranked)",
       subtitle = sprintf("Mean R\u00b2 = %.4f  | dashed = mean", mean(df_r2$R2), na.rm = TRUE),
       x = NULL, y = expression(R^2)) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
print(p_r2_bar)

# --- Plot D: Score Scatter Plot between K ---
library(cowplot)

nat_by_subject <- tapply(as.character(meta$nationality),
                         group_vec, function(x) x[1])

df_f3 <- data.frame(
  f1  = fit_real$f_hat[1, ],
  f2  = fit_real$f_hat[2, ],
  f3  = fit_real$f_hat[3, ],
  nat = as.factor(nat_by_subject)
)

cols <- c("AFR" = "#E07B54", "FIN" = "#4A90D9")

make_panel <- function(x, y, xlab, ylab) {
  ggplot(df_f3, aes(x = .data[[x]], y = .data[[y]], color = nat)) +
    geom_point(size = 2.5) +
    scale_color_manual(values = cols,
                       labels = c("AFR" = "African", "FIN" = "Finnish")) +
    labs(x = xlab, y = ylab, color = "Nationality") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none")
}

p12 <- make_panel("f1", "f2", "Factor 1", "Factor 2")
p13 <- make_panel("f1", "f3", "Factor 1", "Factor 3")
p23 <- make_panel("f2", "f3", "Factor 2", "Factor 3")

legend <- cowplot::get_legend(
  p12 + theme(legend.position = "right")
)

p_all <- (p12 | p13 | p23 | legend) +
  plot_annotation(
    title = "Subject-level factor scores (K = 3)",
    subtitle = "Coloured by nationality"
  )

print(p_all)
ggsave("Rplot_fscores_K3.pdf", p_all, width = 9, height = 3)