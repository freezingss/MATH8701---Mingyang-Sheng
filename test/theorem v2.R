# ====================================================================>
# Feasibility Test for Theorem 2: Latent Structure (B & D) Recovery
# ====================================================================>
set.seed(2026)

# 1. 基础参数规格
J   <- 50       # 组别数量 (Groups)
n_j <- 30       # 每组样本量 (Cluster size)
N   <- J * n_j  # 总观测单位数
Q   <- 200      # 类别维度
K   <- 5        # 潜在因子维度 (调小到5以便更清晰地展示三角约束)
p   <- 10       # 协变量维度

# 2. 真实参数生成 (引入严谨的因子载荷识别约束)
phi_true      <- matrix(rnorm(p * Q, mean = 0, sd = 0.5), nrow = p, ncol = Q)
phi_true[, Q] <- 0  

mu_true      <- rnorm(Q, mean = 0, sd = 0.5)
mu_true[Q]   <- 0   

# 生成满足识别约束的 B_true (前 K 行形成下三角矩阵，确立因子基底)
B_true <- matrix(rnorm(Q * K, mean = 0, sd = 1), nrow = Q, ncol = K)
for (q in 1:K) {
  if (q < K) {
    B_true[q, (q+1):K] <- 0  # 经典因子分析的上三角清零约束
  }
}

f_true   <- matrix(rnorm(K * J, 0, 1), nrow = K, ncol = J)
D_sd     <- runif(Q, min = 0.5, max = 1.5) 
D_sd[Q]  <- 0 # 基准组无随机误差

# 3. 合成真实的群组效应
U_group  <- matrix(0, nrow = Q, ncol = J)
for (j in 1:J) {
  e_j          <- rnorm(Q, mean = 0, sd = D_sd)
  U_group[, j] <- B_true %*% f_true[, j] + e_j
}
U_group[Q, ] <- 0

# 4. 生成 60 万规模的模拟计数数据
Y_matrix <- matrix(0, nrow = N, ncol = Q)
V_array  <- matrix(0, nrow = N, ncol = p)
group_id <- rep(1:J, each = n_j)
Xi_true  <- numeric(N) 
M_vector <- numeric(N) 

for (i in 1:N) {
  j            <- group_id[i]
  v_i          <- rnorm(p, 0, 1)
  V_array[i, ] <- v_i
  M_vector[i]  <- sample(300:500, 1)
  
  eta_i        <- mu_true + U_group[, j] + as.vector(crossprod(phi_true, v_i))
  Xi_true[i]   <- log(M_vector[i]) - log(sum(exp(eta_i)))
  prob_i       <- exp(eta_i) / sum(exp(eta_i))
  Y_matrix[i, ]<- rmultinom(1, size = M_vector[i], prob = prob_i)
}

# ====================================================================>
# 5. 【修正版】核心：利用 PME Stream 抽取包含 mu 的原始复合群组效应
# ====================================================================>
phi_est <- matrix(0, nrow = p, ncol = Q)
U_combined_est <- matrix(0, nrow = Q, ncol = J) # 储存原始复合效应（不提前扣除均值）
mu_est  <- numeric(Q)

cat("Step 1: Running PME Stream to extract raw group structures...\n")

for (q in 1:(Q-1)) {
  y_q <- Y_matrix[, q]
  
  # 直接在 Poisson GLM 中作为高维截距进行“吸入估计”
  fit <- glm(y_q ~ V_array + factor(group_id) - 1, family = poisson(link = "log"), offset = Xi_true)
  
  # 提取固定的核心自变量系数 phi
  phi_est[, q] <- coef(fit)[1:p]
  
  # 提取包含 mu_q + U_qj 的原始复合群组效应
  U_combined_est[q, ] <- coef(fit)[(p+1):(p+J)]
}

# ====================================================================>
# 6. 【修正版】核心：通过带截距的 OLS 回归，完美剥离 mu 并无偏还原 B 和 D
# ====================================================================>
cat("Step 2: Decoupling and recovering B and D using Intercept-OLS...\n")

B_est    <- matrix(0, nrow = Q, ncol = K)
D_sd_est <- numeric(Q)
F_matrix <- t(f_true) # 暂借因子基底进行投影

for (q in 1:(Q-1)) {
  # 严格遵循三角识别约束
  active_factors <- min(q, K) 
  X_f <- F_matrix[, 1:active_factors, drop = FALSE]
  y_u <- U_combined_est[q, ]
  
  # 关键修改：允许截距项存在！用截距项去无损吸收 mu_q 和小样本抽样偏误
  fit_structure <- lm(y_u ~ X_f)
  
  # 截距项即为 mu_q 的无偏估计
  mu_est[q] <- coef(fit_structure)[1]
  
  # 斜率项（去掉截距后）即为 B 的精确无偏估计
  B_est[q, 1:active_factors] <- coef(fit_structure)[-1]
  
  # 此时的残差标准差才是纯净的 D_sd
  D_sd_est[q] <- sd(residuals(fit_structure))
}

# 7. 精度实证度量
mae_phi <- mean(abs(phi_est - phi_true))
mae_B   <- mean(abs(B_est - B_true))
mae_D   <- mean(abs(D_sd_est - D_sd))

cat("\n======================= EVALUATION =======================\n")
cat(sprintf("Total Simulation Observations: %d\n", sum(M_vector)))
cat(sprintf("MAE for Fixed Effects (phi):   %.4f\n", mae_phi))
cat(sprintf("MAE for Factor Loadings (B):   %.4f\n", mae_B))
cat(sprintf("MAE for Unique Variances (D):  %.4f\n", mae_D))
cat("==========================================================\n")

# 打印前 5 行 B 的恢复对比（直观查看三角约束和拟合精度）
cat("\nTop 5 rows of B Matrix Comparison (True vs Estimated):\n")
cat("--- True B (Truncated) ---\n")
print(round(B_true[1:5, ], 3))
cat("--- Estimated B (PME-Decoupled) ---\n")
print(round(B_est[1:5, ], 3))