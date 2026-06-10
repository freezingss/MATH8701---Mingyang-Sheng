# ==============================================================================
# 正统因子分析终极闭环版：基于 PME-EM 架构（含文献双重优化与双指标监控）
# 识别规范：Effects Coding (Little 2006) 扩展版 —— 平均绝对值全局锚定
# ==============================================================================
set.seed(2026)

# ------------------------------------------------------------------------------
# 1. 模拟数据生成维度与正统因子结构设定
# ------------------------------------------------------------------------------
J <- 25         # 群组数量 (Groups)
N_j <- 12       # 每个群组内部的个体数 (Individuals per group)
Q <- 10         # 类别数量 (Categories)
K <- 2          # 潜在因子维度 (Latent Factor Dimension)
Y_total <- 250  # 每个个体的总观测计数

# 定义正统因子分析中的 Error Variance 标量 (D = sigma2 * I)
true_sigma2 <- 0.25 

cat("--- 步骤 1: 正在构建满足上三角约束与残差方差 D 的真实参数 ---\n")

true_mu <- rnorm(Q, mean = 0, sd = 0.5)

# 构建满足【上三角约束】的因子载荷矩阵 B (维度: Q x K)
true_B <- matrix(rnorm(Q * K, mean = 0, sd = 0.6), nrow = Q, ncol = K)
for (q in 1:K) {
  if (q > 1) {
    true_B[q, 1:(q-1)] <- 0 # 下三角区域强行清零以确立几何可识别性
  }
  true_B[q, q] <- abs(true_B[q, q]) + 0.4 # 确保对角线严格为正
}

# 生成隐因子 F (维度: J x K)
true_F <- matrix(rnorm(J * K, mean = 0, sd = 1), nrow = J, ncol = K)

# ==============================================================================
# 数理同步：对真实参数施行 Effects Coding 尺度变换，使后续 F-norm 监控具备绝对可比性
# ==============================================================================
for (k in 1:K) {
  c_true <- mean(abs(true_B[, k]))
  if (c_true < 1e-6) c_true <- 1e-6
  true_B[, k] <- true_B[, k] / c_true
  true_F[, k] <- true_F[, k] * c_true
}

# 生成个体-类别层面的特异性残差 E (维度: J x N_j x Q)
true_E <- array(rnorm(J * N_j * Q, mean = 0, sd = sqrt(true_sigma2)), dim = c(J, N_j, Q))

# 初始化观测张量 Y
Y <- array(0, dim = c(J, N_j, Q))
for (j in 1:J) {
  for (i in 1:N_j) {
    eta <- true_mu + true_B %*% true_F[j, ] + true_E[j, i, ]
    prob <- exp(eta) / sum(exp(eta))
    Y[j, i, ] <- rmultinom(1, size = Y_total, prob = prob)
  }
}
Y_plus <- apply(Y, c(1, 2), sum)

# ------------------------------------------------------------------------------
# 2. EM 算法参数初始化
# ------------------------------------------------------------------------------
cat("--- 步骤 2: 基于文献 (Chiquet 2019) 边际独立 MLE 初始化参数 ---\n")

# 基于边际计数的对数线性解析解进行初始化
total_margin_counts <- apply(Y, 3, sum) 
marginal_probs <- total_margin_counts / sum(total_margin_counts) 

est_mu <- log(marginal_probs) + log(Q) 
est_mu <- est_mu - mean(est_mu) # 中心化消除平移不确定性

est_F      <- matrix(0, nrow = J, ncol = K)   
est_E      <- array(0, dim = c(J, N_j, Q))    # 个体残差模式估计
est_sigma2 <- 0.5                             # 初始残差方差 D 的猜测值
est_Vinv   <- array(0, dim = c(K, K, J))     

# 初始化载荷矩阵 B 并严格绑定上三角零约束
est_B <- matrix(0.1, nrow = Q, ncol = K)
for (q in 1:K) {
  if (q > 1) est_B[q, 1:(q-1)] <- 0
  est_B[q, q] <- 0.5
}

est_xi <- log(Y_plus) - log(Q) 

max_iter <- 15
cat("开始执行引入全局效应编码约束且锁定几何识别性的 EM 迭代...\n\n")

# ------------------------------------------------------------------------------
# 3. EM 算法主循环
# ------------------------------------------------------------------------------
for (iter in 1:max_iter) {
  
  # ==================== E-STEP: 隐因子与解耦残差的交替后验模式估计 ====================
  for (j in 1:J) {
    fj <- est_F[j, ]
    
    # 块坐标下降法 (Block Coordinate Descent) 寻找联合后验 Mode
    for (alt in 1:3) {
      # (步骤 1) 锁定因子 f，极速独立刷新各类别的特异性残差 e
      for (i in 1:N_j) {
        for (q in 1:Q) {
          for (inner_e in 1:2) {
            lambda_ijq <- exp(est_xi[j, i] + est_mu[q] + sum(est_B[q, ] * fj) + est_E[j, i, q])
            grad_e <- Y[j, i, q] - lambda_ijq - est_E[j, i, q] / est_sigma2
            hess_e <- -lambda_ijq - 1 / est_sigma2
            est_E[j, i, q] <- est_E[j, i, q] - grad_e / hess_e
          }
        }
      }
      # (步骤 2) 锁定残差 e，更新低维空间下的群组隐因子 f
      for (inner_f in 1:2) {
        grad_f <- rep(0, K)
        Hess_f <- matrix(0, K, K)
        for (i in 1:N_j) {
          lambda_j <- as.vector(exp(est_xi[j, i] + est_mu + est_B %*% fj + est_E[j, i, ]))
          grad_f <- grad_f + as.vector(t(est_B) %*% (Y[j, i, ] - lambda_j))
          Hess_f <- Hess_f + t(est_B) %*% (lambda_j * est_B)
        }
        grad_f <- grad_f - fj
        Hess_f |> (`+`)(diag(K)) -> Hess_f
        fj <- fj + solve(Hess_f, grad_f)
      }
    }
    est_F[j, ] <- fj
    
    # 基于完全对角化积分，解析计算融合残差后的精确群组后验精度矩阵
    H_j <- diag(K)
    for (i in 1:N_j) {
      lambda_j <- as.vector(exp(est_xi[j, i] + est_mu + est_B %*% fj + est_E[j, i, ]))
      weight   <- lambda_j / (1 + est_sigma2 * lambda_j)
      H_j      <- H_j + t(est_B) %*% (weight * est_B)
    }
    est_Vinv[, , j] <- solve(H_j)
  }
  
  # 计算构建 Q-函数所需的校正方差张量
  var_eps <- array(0, dim = c(J, N_j, Q))
  S       <- array(0, dim = c(J, N_j, Q))
  for (j in 1:J) {
    fj <- est_F[j, ]
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
  
  # ==================== M-STEP: 包含上三角物理约束的解耦更新 ====================
  for (q in 1:Q) {
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
          tilde_lambda <- exp(est_xi[j, i] + mu_c + sum(B_full * fj) + est_E[j, i, q] + 0.5 * S[j, i, q])
          
          g_mu <- Y[j, i, q] - tilde_lambda
          g_B  <- (Y[j, i, q] - tilde_lambda) * fj[free_B_idx]
          g_q  <- g_q + c(g_mu, g_B)
          
          H_q[1, 1] <- H_q[1, 1] - tilde_lambda
          H_q[1, 2:dim_q] <- H_q[1, 2:dim_q] - tilde_lambda * fj[free_B_idx]
          H_q[2:dim_q, 1] <- H_q[2:dim_q, 1] - tilde_lambda * fj[free_B_idx]
          H_q[2:dim_q, 2:dim_q] <- H_q[2:dim_q, 2:dim_q] - tilde_lambda * outer(fj[free_B_idx], fj[free_B_idx])
        }
      }
      theta_q <- theta_q - solve(H_q, g_q)
    }
    est_mu[q] <- theta_q[1]
    if (n_free_B > 0) est_B[q, free_B_idx] <- theta_q[2:dim_q]
  }
  
  # ==================== D MATRIX UPDATE: 残差方差闭式更新 ====================
  est_sigma2 <- mean(est_E^2 + var_eps)
  
  # ==================== SCALE ANCHORING: 效应编码全局绝对值归一化 ====================
  for (k in 1:K) {
    c_k <- mean(abs(est_B[, k]))
    if (c_k < 1e-6) c_k <- 1e-6 
    
    est_B[, k] <- est_B[, k] / c_k
    est_F[, k] <- est_F[, k] * c_k
  }
  
  # ==================== SIGN ALIGNMENT: 符号轴向锚定 ====================
  for (k in 1:K) {
    if (est_B[k, k] < 0) {
      est_B[, k] <- -est_B[, k]
      est_F[, k] <- -est_F[, k]
    }
  }
  
  # ==================== XI UPDATE: 记账流动控制更新 ====================
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
  
  # ==================== 4. 真正全指标精度全景监控 (含 F 范数距离) ====================
  corr_mu <- cor(true_mu, est_mu)
  corr_B  <- cor(as.vector(true_B), as.vector(est_B))
  corr_F  <- cor(as.vector(true_F), as.vector(est_F))
  
  fnorm_mu <- sqrt(sum((true_mu - est_mu)^2))
  fnorm_B  <- sqrt(sum((true_B - est_B)^2))
  fnorm_F  <- sqrt(sum((true_F - est_F)^2))
  
  cat(sprintf("[迭代 %02d]\n", iter))
  cat(sprintf("  ▶ 皮尔逊相关度  -> 截距 mu: %.4f | 载荷 B: %.4f | 隐因子 f: %.4f\n", corr_mu, corr_B, corr_F))
  cat(sprintf("  ▶ Frobenius范数 -> 截距 mu: %.4f | 载荷 B: %.4f | 隐因子 f: %.4f\n", fnorm_mu, fnorm_B, fnorm_F))
  cat(sprintf("  ▶ 残差方差项 D  -> 真实值: %.2f   | 估计值: %.4f\n", true_sigma2, est_sigma2))
  cat("--------------------------------------------------------------------------\n")
}

cat("\n--- 运行结束: 正统因子模型确认闭环 ---\n")