# ====================================================================
# Monte Carlo Simulation: Multinomial Count Data Framework (Size > 1)
# ====================================================================
set.seed(2026)

# 1. Parameter Specifications
J   <- 50       # Number of groups
n_j <- 30       # Cluster size
N   <- J * n_j  # Number of observation units
Q   <- 200       # Category dimension
K   <- 8        # Latent factor dimension
p   <- 20        # Covariate dimension

# 2. True Parameter Generation
phi_true      <- matrix(rnorm(p * Q, mean = 0, sd = 1), nrow = p, ncol = Q)
phi_true[, Q] <- 0  # Identification constraint for baseline category

mu_true      <- rnorm(Q, mean = 0, sd = 0.5)
mu_true[Q]   <- 0   # Identification constraint for baseline intercept

B_true   <- matrix(rnorm(Q * K, 0, 1), nrow = Q, ncol = K)
f_true   <- matrix(rnorm(K * J, 0, 1), nrow = K, ncol = J)
D_sd     <- runif(Q, min = 0.5, max = 1.5) 

# 3. Generating Group-Specific Total Random Effects
U_group  <- matrix(0, nrow = Q, ncol = J)
for (j in 1:J) {
  e_j          <- rnorm(Q, mean = 0, sd = D_sd)
  U_group[, j] <- B_true %*% f_true[, j] + e_j
}

# 4. Synthesizing Observable Count Matrix (M_ij > 1)
Y_matrix <- matrix(0, nrow = N, ncol = Q)
V_array  <- matrix(0, nrow = N, ncol = p)
group_id <- rep(1:J, each = n_j)
Xi_true  <- numeric(N) 
M_vector <- numeric(N) # Store individual total trial counts

for (i in 1:N) {
  j            <- group_id[i]
  v_i          <- rnorm(p, 0, 1)
  V_array[i, ] <- v_i
  
  # Define individual-specific large choice capacities (e.g., between 15 and 40 trials)
  M_vector[i]  <- sample(300:500, 1)
  
  # Compute raw utility field
  eta_i        <- mu_true + U_group[, j] + as.vector(crossprod(phi_true, v_i))
  
  # Extract PME profile parameter incorporating log(M_i)
  Xi_true[i]   <- log(M_vector[i]) - log(sum(exp(eta_i)))
  
  prob_i       <- exp(eta_i) / sum(exp(eta_i))
  
  # Generate non-negative integer counts instead of binary vectors
  Y_matrix[i, ]<- rmultinom(1, size = M_vector[i], prob = prob_i)
}

# 5. Estimation via Profile Poisson GLM Stream
phi_est <- matrix(0, nrow = p, ncol = Q)

for (q in 1:(Q-1)) { 
  y_q        <- Y_matrix[, q]
  offset_val <- numeric(N)
  
  for (i in 1:N) {
    j             <- group_id[i]
    offset_val[i] <- U_group[q, j] + Xi_true[i]
  }
  
  # Run the decoupled conditional Poisson GLM on integer counts
  fit          <- glm(y_q ~ V_array, family = poisson(link = "log"), offset = offset_val)
  phi_est[, q] <- coef(fit)[-1] 
}
phi_est[, Q] <- 0 

# 6. Empirical Evaluation Metrics
mae_error <- mean(abs(phi_est - phi_true), na.rm = TRUE)
cat("====================================================\n")
cat(sprintf("Total Choice Trials (M_total): %d\n", sum(M_vector)))
cat(sprintf("Empirical Mean Absolute Error (MAE) for phi: %.4f\n", mae_error))
cat("====================================================\n")

cat("\nSample Comparison (True vs Estimated for Category 1):\n")
print(rbind(True = phi_true[, 1], Estimated = phi_est[, 1]))

# ====================================================================>
# 7. 期刊级参数恢复可视化 (Theorem 1 Validation)
# ====================================================================>
# 确保安装并加载了 ggplot2 和 patchwork（用于合并图）
# install.packages(c("ggplot2", "patchwork"))
library(ggplot2)
library(patchwork)

# 1. 将 4000 个高维参数拉平并转化为清洁的数据框
plot_df <- data.frame(
  True_Value      = as.vector(phi_true),
  Estimated_Value = as.vector(phi_est)
)
plot_df$Error <- plot_df$Estimated_Value - plot_df$True_Value

# 计算拟合优度 R2 和 RMSE 作为图上标注
r2_val   <- cor(plot_df$True_Value, plot_df$Estimated_Value)^2
rmse_val <- sqrt(mean(plot_df$Error^2))

# 图 A：45度对角线回归契合图
p1 <- ggplot(plot_df, aes(x = True_Value, y = Estimated_Value)) +
  # 绘制六边形热力密度（4000个点很密，用hex或透明度点避免重叠）
  geom_point(alpha = 0.4, color = "#1f77b4", size = 1.5) +
  # 完美的 45 度参考线
  geom_abline(intercept = 0, slope = 1, color = "#d62728", linetype = "dashed", size = 1) +
  theme_bw(base_size = 14) +
  labs(
    x = "True Parameter (Baseline)",
    y = "Estimated Parameter (PME-Stream)",
    title = "(a) Parameter Recovery Alignment"
  ) +
  # 在图上优雅地打上 R² 和 MAE 标签
  annotate("text", x = min(plot_df$True_Value)*0.7, y = max(plot_df$Estimated_Value)*0.8, 
           label = sprintf("R² = %.4f\nMAE = %.4f\nRMSE = %.4f", r2_val, mae_error, rmse_val),
           size = 4.5, fontface = "bold", hjust = 0, bg.r = 0.15, bg.color = "white") +
  theme(plot.title = element_text(face = "bold", size = 14))

# 图 B：估计残差的分布图（证明无偏性）
p2 <- ggplot(plot_df, aes(x = Error)) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "#2ca02c", alpha = 0.6, color = "white") +
  geom_density(color = "#1f77b4", size = 1) +
  geom_vline(xintercept = 0, color = "#d62728", linetype = "solid", size = 0.8) +
  theme_bw(base_size = 14) +
  labs(
    x = "Estimation Error (Est - True)",
    y = "Density",
    title = "(b) Distribution of Errors"
  ) +
  theme(plot.title = element_text(face = "bold", size = 14))

# 终极合并展示
combined_plot <- p1 + p2 + plot_layout(ncol = 2)
print(combined_plot)

# ===================== TRY 2: Do not use true value =====================
# ====================================================================
# Monte Carlo Simulation: Multinomial Count Data Framework (Size > 1)
# Strictly Non-Cheating Framework for Theorem 1 Validation
# ====================================================================
set.seed(2026)

# 1. Parameter Specifications
J   <- 50       # Number of groups
n_j <- 30       # Cluster size
N   <- J * n_j  # Number of observation units
Q   <- 200      # Category dimension
K   <- 8        # Latent factor dimension
p   <- 20       # Covariate dimension

# 2. True Parameter Generation
phi_true      <- matrix(rnorm(p * Q, mean = 0, sd = 1), nrow = p, ncol = Q)
phi_true[, Q] <- 0  # Identification constraint for baseline category

mu_true      <- rnorm(Q, mean = 0, sd = 0.5)
mu_true[Q]   <- 0   # Identification constraint for baseline intercept

B_true   <- matrix(rnorm(Q * K, 0, 1), nrow = Q, ncol = K)
f_true   <- matrix(rnorm(K * J, 0, 1), nrow = K, ncol = J)
D_sd     <- runif(Q, min = 0.5, max = 1.5) 

# 3. Generating Group-Specific Total Random Effects
U_group  <- matrix(0, nrow = Q, ncol = J)
for (j in 1:J) {
  e_j          <- rnorm(Q, mean = 0, sd = D_sd)
  U_group[, j] <- B_true %*% f_true[, j] + e_j
}

# 4. Synthesizing Observable Count Matrix
Y_matrix <- matrix(0, nrow = N, ncol = Q)
V_array  <- matrix(0, nrow = N, ncol = p)
group_id <- rep(1:J, each = n_j)
M_vector <- numeric(N) 

for (i in 1:N) {
  j            <- group_id[i]
  v_i          <- rnorm(p, 0, 1)
  V_array[i, ] <- v_i
  
  # 遵循 Claude 建议：降低 Trial 计数至合理区间，制造更真实的统计噪声
  M_vector[i]  <- sample(50:150, 1)
  
  eta_i        <- mu_true + U_group[, j] + as.vector(crossprod(phi_true, v_i))
  prob_i       <- exp(eta_i) / sum(exp(eta_i))
  Y_matrix[i, ]<- rmultinom(1, size = M_vector[i], prob = prob_i)
}

# ====================================================================
# 5. Estimation via Profile Poisson GLM Stream (STRICTLY NON-CHEATING)
# ====================================================================
# 【彻底修正】严禁使用任何 U_group 或 Xi_true，全部从 0 开始盲跑
phi_est <- matrix(0, nrow = p, ncol = Q)
U_est   <- matrix(0, nrow = Q, ncol = J) # 自由更新的分组效应矩阵
xi_est  <- rep(0, N)                     # 动态剖面拦截项

# 预先构建设计矩阵（包含协变量与分组 Dummy），大幅压榨 glm.fit 的性能
X_design <- cbind(V_array, model.matrix(~ factor(group_id) - 1))

max_iter <- 15
tol      <- 1e-3

cat("开始执行真实的交替剖面似然估计循环（无先知信息借调）...\n")
for (iter in 1:max_iter) {
  phi_old <- phi_est
  
  # --- 步骤 A: PME Profiling Step (动态集中消去拦截项) ---
  # 仅基于当前轮次的 phi_est 和 U_est 计算线性预测子
  Eta <- V_array %*% phi_est + t(U_est[, group_id])
  # 动态投影出这一轮的 xi 估计值
  xi_est <- log(M_vector) - log(rowSums(exp(Eta)))
  
  # --- 步骤 B: 类别解耦的 M-step (通过底层 glm.fit 提速) ---
  for (q in 1:(Q - 1)) {
    # 将当前轮次的 xi_est 作为唯一的已知 offset 喂给模型
    fit <- glm.fit(X_design, Y_matrix[, q], family = poisson(link = "log"), offset = xi_est)
    
    # 完美拆分：前 p 个是结构参数 phi，后 J 个是该类别的分组截距
    phi_est[, q] <- fit$coefficients[1:p]
    U_est[q, ]   <- fit$coefficients[(p + 1):(p + J)]
  }
  
  # 强行施加基准类别（Baseline Category Q）的识别约束
  phi_est[, Q] <- 0
  U_est[Q, ]   <- 0
  
  # 计算参数估计的平均绝对变化量，判断收敛
  param_diff <- mean(abs(phi_est - phi_old), na.rm = TRUE)
  cat(sprintf("Iteration %d/%d - Parameter Shift: %.6f\n", iter, max_iter, param_diff))
  
  if (param_diff < tol) {
    cat("算法顺利收敛！\n\n")
    break
  }
}

# ====================================================================
# 6. Empirical Evaluation Metrics
# ====================================================================
# 此时计算出来的 MAE 才是真正具备统计学说服力的结果
mae_error <- mean(abs(phi_est - phi_true), na.rm = TRUE)
cat("====================================================\n")
cat(sprintf("Total Choice Trials (M_total): %d\n", sum(M_vector)))
cat(sprintf("REAL Empirical Mean Absolute Error (MAE) for phi: %.4f\n", mae_error))
cat("====================================================\n")

# ====================================================================
# 7. 期刊级参数恢复可视化 (Theorem 1 Validation)
# ====================================================================
library(ggplot2)
library(patchwork)

plot_df <- data.frame(
  True_Value      = as.vector(phi_true),
  Estimated_Value = as.vector(phi_est)
)
plot_df$Error <- plot_df$Estimated_Value - plot_df$True_Value

r2_val   <- cor(plot_df$True_Value, plot_df$Estimated_Value)^2
rmse_val <- sqrt(mean(plot_df$Error^2))

p1 <- ggplot(plot_df, aes(x = True_Value, y = Estimated_Value)) +
  geom_point(alpha = 0.3, color = "#1f77b4", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, color = "#d62728", linetype = "dashed", size = 1) +
  theme_bw(base_size = 14) +
  labs(
    x = "True Parameter (Phi)",
    y = "Estimated Parameter (Non-Cheating PME)",
    title = "(a) Parameter Recovery Alignment"
  ) +
  annotate("text", x = min(plot_df$True_Value)*0.7, y = max(plot_df$Estimated_Value)*0.8, 
           label = sprintf("R² = %.4f\nMAE = %.4f\nRMSE = %.4f", r2_val, mae_error, rmse_val),
           size = 4.5, fontface = "bold", hjust = 0) +
  theme(plot.title = element_text(face = "bold", size = 12))

p2 <- ggplot(plot_df, aes(x = Error)) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "#2ca02c", alpha = 0.6, color = "white") +
  geom_density(color = "#1f77b4", size = 1) +
  geom_vline(xintercept = 0, color = "#d62728", linetype = "solid", size = 0.8) +
  theme_bw(base_size = 14) +
  labs(
    x = "Estimation Error (Est - True)",
    y = "Density",
    title = "(b) Distribution of Errors"
  ) +
  theme(plot.title = element_text(face = "bold", size = 12))

combined_plot <- p1 + p2 + plot_layout(ncol = 2)
print(combined_plot)
