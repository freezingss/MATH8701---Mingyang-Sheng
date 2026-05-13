lambda_schedule <- function(t, lambda0, kappa, Delta) {
  kappa^t * lambda0 + (1 - kappa^t) / (1 - kappa) * Delta
}

free_B_mask <- function(Q, K) {
  outer(seq_len(Q), seq_len(K), function(q, k) k <= q)
}

logpost_f <- function(f, yj, wj, mu, B, d2) {
  eta      <- mu + as.vector(B %*% f)
  log_lam  <- eta + log(pmax(wj, 1e-12)) + 0.5 * d2
  lam_mean <- exp(pmin(log_lam, 30))          
  val <- sum(yj * eta - lam_mean) - 0.5 * sum(f * f)
  if (!is.finite(val)) return(-1e15)
  val
}

grad_logpost_f <- function(f, yj, wj, mu, B, d2) {
  eta      <- mu + as.vector(B %*% f)
  log_lam  <- eta + log(pmax(wj, 1e-12)) + 0.5 * d2
  lam_mean <- exp(pmin(log_lam, 30))
  as.vector(crossprod(B, yj - lam_mean)) - f
}

hess_logpost_f <- function(f, yj, wj, mu, B, d2) {
  eta      <- mu + as.vector(B %*% f)
  log_lam  <- eta + log(pmax(wj, 1e-12)) + 0.5 * d2
  lam_mean <- exp(pmin(log_lam, 30))
  -crossprod(B, B * lam_mean) - diag(length(f))
}

nr_mode_f <- function(f0, yj, wj, mu, B, d2, maxit = 50, tol = 1e-8) {
  f <- f0
  for (iter in seq_len(maxit)) {
    g <- grad_logpost_f(f, yj, wj, mu, B, d2)
    if (sqrt(sum(g * g)) < tol) break
    H    <- hess_logpost_f(f, yj, wj, mu, B, d2)
    step <- tryCatch(solve(H, g), error = function(e) as.vector(solve(as.matrix(H), g)))
    if (!all(is.finite(step))) break          
    m0   <- logpost_f(f, yj, wj, mu, B, d2)
    alph <- 1
    for (ls in seq_len(20)) {
      f_new <- f - alph * step
      if (logpost_f(f_new, yj, wj, mu, B, d2) > m0 - 1e-8 * alph * sum(g * step)) {
        f <- f_new
        break
      }
      alph <- alph * 0.5
    }
  }
  H_final <- hess_logpost_f(f, yj, wj, mu, B, d2)
  H_final <- H_final - 1e-6 * diag(length(f))
  list(f_hat = f, H = H_final)
}

m1_row_q <- function(q, y_mat, w_mat, f_tilde, d2, lam, mu0, b0) {
  K    <- ncol(f_tilde)
  kmax <- min(q, K)
  yc   <- y_mat[, q]
  wc   <- w_mat[, q]
  if (kmax < 1) {
    mu_q <- log((sum(yc) + 1e-6) / (sum(wc * exp(0.5 * d2[q])) + 1e-6))
    return(list(mu_q = mu_q, b = numeric(0)))
  }
  A   <- f_tilde[, seq_len(kmax), drop = FALSE]
  par <- c(0, rep(0, kmax))

  fn_obj <- function(p) {
    mq     <- p[1]
    b      <- p[-1]
    Ab     <- as.vector(A %*% b)
    lin    <- mq + Ab
    log_lq <- mq + 0.5 * d2[q] + Ab
    lam_q  <- exp(pmin(log_lq, 30)) * wc
    sm     <- sum(yc * lin - lam_q)
    pen    <- lam * sqrt(sum(b * b) + 1e-12)
    if (!is.finite(sm)) return(1e15)
    -(sm - pen)
  }

  gr_obj <- function(p) {
    eps  <- 1e-6
    v    <- numeric(length(p))
    f0   <- fn_obj(p)
    for (i in seq_along(p)) {
      pp     <- p; pp[i] <- pp[i] + eps
      v[i]   <- (fn_obj(pp) - f0) / eps
    }
    v
  }

  o <- tryCatch(
    optim(par, fn_obj, gr = gr_obj, method = "L-BFGS-B",
          lower   = c(-20, rep(-5, kmax)),
          upper   = c( 20, rep( 5, kmax)),
          control = list(maxit = 500, factr = 1e9)),
    error = function(e) list(par = par) 
  )
  list(mu_q = o$par[1], b = o$par[-1])
}

Sq <- function(q, mu, B, f_tilde, Hinv_list, s_mat) {
  J   <- nrow(f_tilde)
  Bq  <- B[q, , drop = FALSE]
  s   <- 0
  for (j in seq_len(J)) {
    t1  <- as.numeric(Bq %*% Hinv_list[[j]] %*% t(Bq))
    lin <- mu[q] + sum(Bq * f_tilde[j, ])
    s   <- s + t1 + (lin - s_mat[j, q])^2
  }
  s
}

ell_smooth_j <- function(yj, wj, mu, B, d2, fj) {
  eta <- mu + log(pmax(wj, 1e-12)) + 0.5 * d2 + as.vector(B %*% fj)
  sum(yj * eta - exp(pmin(eta, 30)))
}

ell_smooth_all <- function(y_mat, w_mat, mu, B, d2, Fm) {
  J <- nrow(y_mat); s <- 0
  for (j in seq_len(J))
    s <- s + ell_smooth_j(y_mat[j,], w_mat[j,], mu, B, d2, Fm[j,])
  s
}

row_l2_sum_B <- function(B, K) {
  Q <- nrow(B); p <- 0
  for (q in seq_len(Q)) {
    km <- min(q, K)
    if (km >= 1) p <- p + sqrt(sum(B[q, seq_len(km)]^2))
  }
  p
}

algorithm2_fit <- function(y_mat, w_mat, K = 2,
                           lambda0 = 0.5, kappa = 0.9, Delta = 0.1,
                           max_iter = 40, tol = 1e-5,
                           verbose = TRUE, trace = FALSE) {
  y_mat <- as.matrix(y_mat)
  w_mat <- as.matrix(w_mat)
  J <- nrow(y_mat); Q <- ncol(y_mat)

  mu <- rep(0, Q)
  B  <- matrix(0, Q, K)
  m  <- free_B_mask(Q, K)
  B[m] <- rnorm(sum(m), sd = 0.05)
  d2 <- rep(0.25, Q)
  s_mat <- sweep(matrix(0, J, Q), 2, mu + 0.5 * d2, "+")

  if (trace) {
    tr  <- data.frame(t = integer(max_iter), lam = numeric(max_iter),
                      lam_geom = numeric(max_iter), lam_floor = numeric(max_iter),
                      dmu = numeric(max_iter), dB = numeric(max_iter),
                      dd2 = numeric(max_iter), ell_smooth = numeric(max_iter),
                      pen = numeric(max_iter), ell_penalized = numeric(max_iter))
    ntr <- 0L
  }

  for (t in seq_len(max_iter)) {
    tau       <- t - 1L
    lam       <- lambda_schedule(tau, lambda0, kappa, Delta)
    lam_geom  <- lambda0 * kappa^tau
    lam_floor <- Delta * (1 - kappa^tau) / (1 - kappa)

    # E-step
    Fm   <- matrix(0, J, K)
    Hinv <- vector("list", J)
    for (j in seq_len(J)) {
      o       <- nr_mode_f(rep(0, K), y_mat[j,], w_mat[j,], mu, B, d2)
      Fm[j,]  <- o$f_hat
      Hinv[[j]] <- tryCatch(solve(-o$H),           
                            error = function(e) diag(K))
    }

    mu_old <- mu; Bo <- B; d2_old <- d2

    # M1
    for (q in seq_len(Q)) {
      km <- min(q, K)
      b0 <- if (km >= 1) B[q, seq_len(km)] else numeric(0)
      r  <- m1_row_q(q, y_mat, w_mat, Fm, d2, lam, mu[q], b0)
      mu[q] <- r$mu_q
      if (km >= 1) {
        B[q, seq_len(km)] <- r$b
        if (K > km) B[q, (km + 1):K] <- 0
      }
    }

    # M2
    for (q in seq_len(Q))
      d2[q] <- max(Sq(q, mu, B, Fm, Hinv, s_mat) / J, 1e-8)

    s_mat <- sweep(Fm %*% t(B), 2, mu + 0.5 * d2, "+")

    dm  <- max(abs(mu - mu_old))
    dB  <- max(abs(B  - Bo))
    dd2 <- max(abs(d2 - d2_old))

    if (trace) {
      es  <- ell_smooth_all(y_mat, w_mat, mu, B, d2, Fm)
      pen <- row_l2_sum_B(B, K)
      ntr <- ntr + 1L
      tr[ntr, ] <- list(t, lam, lam_geom, lam_floor, dm, dB, dd2,
                        es, pen, es - lam * pen)
    }

    if (verbose) message(sprintf("iter %d  lam=%.4f  |dmu|=%.2e  |dB|=%.2e", t, lam, dm, dB))
    if (max(dm, dB) < tol) break
  }

  out <- list(mu = mu, B = B, d2 = d2, f = Fm,
              y_mat = y_mat, w_mat = w_mat, K = K)
  if (trace) out$trace <- tr[seq_len(ntr), , drop = FALSE]
  out
}

plot_algorithm2_convergence <- function(fit, main = "Algorithm 2 diagnostics") {
  if (is.null(fit$trace) || !nrow(fit$trace))
    stop("fit$trace missing; re-run with trace = TRUE")
  tr  <- fit$trace
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))

  plot(tr$t, tr$ell_smooth, type = "l", col = "darkblue",
       xlab = "Iteration", ylab = "log-lik (smooth)",
       main = "Group-level Poisson log-likelihood")
  lines(tr$t, tr$ell_penalized, col = "coral", lty = 2)
  legend("bottomright", c("ell_smooth", "ell_penalized"),
         col = c("darkblue","coral"), lty = c(1,2), bty = "n", cex = 0.8)

  plot(tr$t, tr$dmu, type = "l", col = "black",
       xlab = "Iteration", ylab = "max abs change",
       main = "Parameter increments",
       ylim = range(c(tr$dmu, tr$dB, tr$dd2), na.rm = TRUE))
  lines(tr$t, tr$dB,  col = "darkgreen")
  lines(tr$t, tr$dd2, col = "purple")
  legend("topright", c("dmu","dB","dd2"),
         col = c("black","darkgreen","purple"), lty = 1, bty = "n", cex = 0.8)

  plot(tr$t, tr$lam, type = "l", col = "black",
       xlab = "Iteration", ylab = expression(lambda^(t)),
       main = "Penalty schedule")
  lines(tr$t, tr$lam_geom,  col = "red")
  lines(tr$t, tr$lam_floor, col = "blue")
  legend("bottomright",
         c("lambda","kappa^t*lambda0","floor"),
         col = c("black","red","blue"), lty = 1, bty = "n", cex = 0.75)

  plot(tr$t, tr$pen, type = "l", col = "brown",
       xlab = "Iteration", ylab = "sum row L2",
       main = "Sum of row L2 norms of B")
  mtext(main, outer = TRUE, cex = 1.1, font = 2)
  invisible(tr)
}

simulate_algo2_dataset <- function(J = 60, Q = 12, K = 3, w_const = 5,
                                   mu = NULL, B = NULL, d2 = NULL, seed = 1) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(mu)) mu <- rnorm(Q, 0, 0.2)
  if (is.null(d2)) d2 <- rep(0.08, Q)
  if (is.null(B)) {
    B <- matrix(0, Q, K)
    for (q in seq_len(Q))
      for (k in seq_len(min(q, K))) B[q, k] <- runif(1, 0.15, 0.45)
  }
  f     <- matrix(rnorm(J * K), J, K)
  w_mat <- matrix(w_const, J, Q)
  y_mat <- matrix(0L, J, Q)
  for (j in seq_len(J)) {
    eta <- mu + log(w_mat[j,]) + 0.5 * d2 + as.vector(B %*% f[j,])
    for (q in seq_len(Q)) y_mat[j, q] <- rpois(1, exp(eta[q]))
  }
  list(y_mat = y_mat, w_mat = w_mat,
       truth = list(mu = mu, B = B, d2 = d2, f = f, J = J, Q = Q, K = K))
}

run_algo2_demo <- function() {
  dat <- simulate_algo2_dataset(J = 60, Q = 12, K = 3, seed = 7)
  fit <- algorithm2_fit(
    dat$y_mat, dat$w_mat,
    K = 3, lambda0 = 0.35, kappa = 0.92, Delta = 0.06,
    max_iter = 80, trace = TRUE, verbose = TRUE
  )
  plot_algorithm2_convergence(fit)
  invisible(list(dat = dat, fit = fit))
}

run_algo2_demo_no_penalty <- function() {
  dat <- simulate_algo2_dataset(J = 60, Q = 12, K = 3, seed = 7)
  fit <- algorithm2_fit(
    dat$y_mat, dat$w_mat,
    K = 3, lambda0 = 0, kappa = 0.92, Delta = 0,
    max_iter = 80, trace = TRUE, verbose = TRUE
  )
  plot_algorithm2_convergence(fit)
  invisible(list(dat = dat, fit = fit))
}
