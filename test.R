x <- c(2,4,5,9,6,3,7,5,4)
group <- rep(c("A","B"), times = c(5,4))
t.test(x ~ group)
# ANOVA
model <- lm(x ~group)
anova(model)

# ----------------------------
# Monte Carlo integration helper
# ----------------------------
mc_int_01 <- function(f, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- runif(n, 0, 1)
  fx <- f(x)
  est <- mean(fx)              # since interval length is 1
  se  <- sd(fx) / sqrt(n)      # Monte Carlo standard error
  ci  <- est + c(-1, 1) * 1.96 * se
  list(n = n, estimate = est, se = se, ci_low = ci[1], ci_high = ci[2])
}

# ----------------------------
# Q2: integral_0^1 cos(2*pi*x) dx  (exact = 0)
# ----------------------------
f2 <- function(x) cos(2*pi*x)

res2_100  <- mc_int_01(f2, n = 100,  seed = 1)
res2_1000 <- mc_int_01(f2, n = 1000, seed = 1)

exact2 <- 0

cat("Q2: Integral_0^1 cos(2*pi*x) dx\n")
cat("Exact value =", exact2, "\n\n")

print(as.data.frame(rbind(
  res2_100,
  res2_1000
)))

cat("\nAbsolute errors:\n")
cat("n=100  :", abs(res2_100$estimate  - exact2), "\n")
cat("n=1000 :", abs(res2_1000$estimate - exact2), "\n\n")


# ----------------------------
# Q3: integral_0^1 cos(2*pi*x^2) dx  (no closed form)
# ----------------------------
f3 <- function(x) cos(2*pi*x^2)

res3_100  <- mc_int_01(f3, n = 100,  seed = 1)
res3_1000 <- mc_int_01(f3, n = 1000, seed = 1)

cat("Q3: Integral_0^1 cos(2*pi*x^2) dx (no closed-form exact answer)\n\n")

print(as.data.frame(rbind(
  res3_100,
  res3_1000
)))