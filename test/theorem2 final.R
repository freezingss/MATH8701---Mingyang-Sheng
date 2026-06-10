# ====================================================================>
# Feasibility Test for Theorem 2: Calibrated Latent Structure Recovery
# ====================================================================>
set.seed(2026)

# 1. 基础参数规格
J   <- 50       # 组别数量 (Groups)
n_j <- 30       # 每组样本量 (Cluster size)
N   <- J * n_j  # 总观测单位数
Q   <- 100      # 类别维度 (调到100以加快仿真运行速度)
K   <- 4        # 潜在因子维度
p   <- 5        # 协变量维度

# 2. 真实参数生成 (校准标准差，消除0计数导致的数值分离)
phi_true      <- matrix(rnorm(p * Q, mean = 0, sd = 0.2), nrow = p, ncol = Q)
phi_true[, Q] <- 0  

mu_true      <- rnorm(Q, mean = 0, sd = 0.2)
mu_true[Q]   <- 0   

# 生成满足严格下三角约束的 B_true
B_true <- matrix(rnorm(Q * K, mean = 0, sd = 0.4), nrow = Q, ncol = K)
for (q in 1:K) {
  if (q < K) {
    B_true[q, (q+1):K] <- 0  
  }
}

f_true   <- matrix(rnorm(K * J, 0, 1), nrow = K, ncol = J)
D_sd     <- runif(Q, min = 0.2, max = 0.5) 
D_sd[Q]  <- 0 

# 3. 合成群组效应
U_group  <- matrix(0, nrow = Q, ncol = J)
for (j in 1:J) {
  e_j          <- rnorm(Q, mean = 0, sd = D_sd)
  U_group[, j] <- B_true %*% f_true[, j] + e_j
}
U_group[Q, ] <- 0

# 4. 生成高质量非稀疏计数数据
Y_matrix <- matrix(0, nrow = N, ncol = Q)
V_array  <- matrix(0, nrow = N, ncol = p)
group_id <- rep(1:J, each = n_j)
Xi_true  <- numeric(N) 
M_vector <- numeric(N) 

for (i in 1:N) {
  j            <- group_id[i]
  v_i          <- rnorm(p, 0, 1)
  V_array[i, ] <- v_i
  M_vector[i]  <- sample(500:800, 1) # 略微提升总试验次数，强化统计表现
  
  eta_i        <- mu_true + U_group[, j] + as.vector(crossprod(phi_true, v_i))
  Xi_true[i]   <- log(M_vector[i]) - log(sum(exp(eta_i)))
  prob_i       <- exp(eta_i) / sum(exp(eta_i))
  Y_matrix[i, ]<- rmultinom(1, size = M_vector[i], prob = prob_i)
}

# 5. PME Stream 提取复合群组效应
phi_est <- matrix(0, nrow = p, ncol = Q)
U_combined_est <- matrix(0, nrow = Q, ncol = J) 
mu_est  <- numeric(Q)

cat("Step 1: Running PME Stream to extract structures...\n")
for (q in 1:(Q-1)) {
  y_q <- Y_matrix[, q]
  fit <- glm(y_q ~ V_array + factor(group_id) - 1, family = poisson(link = "log"), offset = Xi_true)
  phi_est[, q] <- coef(fit)[1:p]
  U_combined_est[q, ] <- coef(fit)[(p+1):(p+J)]
}

# 6. 带截距项的结构解耦
cat("Step 2: Decoupling and recovering B and D...\n")
B_est    <- matrix(0, nrow = Q, ncol = K)
D_sd_est <- numeric(Q)
F_matrix <- t(f_true) 

for (q in 1:(Q-1)) {
  active_factors <- min(q, K) 
  X_f <- F_matrix[, 1:active_factors, drop = FALSE]
  y_u <- U_combined_est[q, ]
  
  fit_structure <- lm(y_u ~ X_f)
  mu_est[q] <- coef(fit_structure)[1]
  B_est[q, 1:active_factors] <- coef(fit_structure)[-1]
  D_sd_est[q] <- sd(residuals(fit_structure))
}

# 7. 精度评估
mae_phi <- mean(abs(phi_est - phi_true))
mae_B   <- mean(abs(B_est - B_true))
mae_D   <- mean(abs(D_sd_est - D_sd))

cat("\n======================= EVALUATION =======================\n")
cat(sprintf("MAE for Fixed Effects (phi):   %.4f\n", mae_phi))
cat(sprintf("MAE for Factor Loadings (B):   %.4f\n", mae_B))
cat(sprintf("MAE for Unique Variances (D):  %.4f\n", mae_D))
cat("==========================================================\n")

cat("\nTop 5 rows of B Matrix Comparison (True vs Estimated):\n")
cat("--- True B ---\n"); print(round(B_true[1:5, ], 3))
cat("--- Estimated B ---\n"); print(round(B_est[1:5, ], 3))

plot(D_sd, D_sd_est, main = "D Matrix Comparison: True vs Estimated",
     xlab = "True D", ylab = "Estimated D", pch = 19, col = "blue")
abline(0, 1, col = "red", latch = 2)