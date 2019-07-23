test_that("linear regression works", {
  set.seed(1)
  n <- 1000; p <- 10; groupsize <- 5
  index <- ceiling(1:p / groupsize)
  beta <- c(-1:1)
  x <- matrix(rnorm(n * p), ncol = p, nrow = n)
  y <- as.vector(x[,1:3] %*% beta + 0.1 * rnorm(n))
  fit0 <- lm.fit(x, y)
  fit1 <- asgl(x, y, index, family = "gaussian", lambda_min = 1e-9)
  coef_lm <- as.vector(coef(fit0)); coef_asgl <- fit1$beta[,20]
  expect_equal(coef_lm, coef_asgl, tolerance = 1e-8)
})

test_that("logistic regression works", {
  set.seed(1)
  n <- 1000; p <- 10; groupsize <- 5
  index <- ceiling(1:p / groupsize)
  beta <- c(-1:1)
  x <- matrix(rnorm(n * p), ncol = p, nrow = n)
  eta <- -1 + x[, 1:3] %*% beta
  prob <- exp(eta) / (1 + exp(eta))
  y <- rbinom(n, 1, prob)
  fit0 <- glm.fit(cbind(1, x), y, family = binomial())
  fit1 <- asgl(x, y, index, family = "binomial", lambda_min = 1e-9)
  coef_glm <- as.vector(coef(fit0))
  coef_asgl <- c(fit1$intercept[20], fit1$beta[,20])
  expect_equal(coef_glm, coef_asgl, tolerance = 1e-2)
})

