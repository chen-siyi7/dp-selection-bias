
##  Reproduces Tables 1, 2, 3, 4 of the manuscript and verifies
##  Proposition 4 numerically.

## Computes bar_a_n (the closed-form posterior mean of a(G) under
## DP conjugacy applied to the latent rates of the observed
## species) and the induced plug-in for f_0.
##
naive_dp <- function(obs, alpha = 5, kappa_0 = 1, beta_0 = 1) {
  K_n   <- length(obs)
  L_G0  <- (beta_0 / (beta_0 + 1))^kappa_0
  w     <- ((beta_0 + 1) / (beta_0 + 2))^(kappa_0 + obs)
  bar_a <- (alpha * L_G0 + sum(w)) / (alpha + K_n)
  list(bar_a  = bar_a,
       hat_f0 = K_n * bar_a / (1 - bar_a))
}


## Naive PY plug-in estimator
## Pitman-Yor generalization
##
naive_py <- function(obs, theta = 5, sigma = 0,
                      kappa_0 = 1, beta_0 = 1) {
  K_n  <- length(obs)
  L_G0 <- (beta_0 / (beta_0 + 1))^kappa_0
  w    <- ((beta_0 + 1) / (beta_0 + 2))^(kappa_0 + obs)
  bar_a <- (theta + K_n * sigma) / (theta + K_n) * L_G0 +
           (1 - sigma) / (theta + K_n) * sum(w)
  list(bar_a  = bar_a,
       hat_f0 = K_n * bar_a / (1 - bar_a))
}



## Correct specification (Proposition 1, Corollary 2).
asymp_correct <- function(kappa_star, beta_star) {
  a_star    <- (beta_star / (beta_star + 1))^kappa_star
  L2_star   <- (beta_star / (beta_star + 2))^kappa_star
  bar_a_inf <- (a_star - L2_star) / (1 - a_star)
  rho       <- (1 - a_star) * (a_star - L2_star) /
               (a_star * (1 - 2 * a_star + L2_star))
  list(a_star = a_star, L2_star = L2_star,
       bar_a_inf = bar_a_inf, rho = rho)
}

## Misspecified G_0 (Proposition 3).
asymp_misspec <- function(kappa_star, beta_star, kappa_0, beta_0) {
  a_star  <- (beta_star / (beta_star + 1))^kappa_star
  s       <- 1 / (beta_0 + 2)
  r       <- (beta_0 + 1) / (beta_0 + 2)
  L_s     <- (beta_star / (beta_star + s))^kappa_star
  bar_a_inf <- (r^kappa_0) * (L_s - a_star) / (1 - a_star)
  if (bar_a_inf >= 1 || bar_a_inf <= 0) {
    rho <- NA
  } else {
    rho <- (1 - a_star) * bar_a_inf / (1 - bar_a_inf) / a_star
  }
  list(bar_a_inf = bar_a_inf, rho = rho)
}

## Pitman--Yor (Proposition 4).
asymp_py <- function(kappa_star, beta_star, sigma) {
  cs        <- asymp_correct(kappa_star, beta_star)
  bar_a_inf <- sigma * cs$a_star + (1 - sigma) * cs$bar_a_inf
  if (bar_a_inf >= 1) {
    rho <- NA
  } else {
    rho <- (1 - cs$a_star) * bar_a_inf /
           (1 - bar_a_inf) / cs$a_star
  }
  list(bar_a_inf = bar_a_inf, rho = rho)
}


## Compound Poisson simulation
simulate_cp <- function(kappa_star, beta_star, Lambda) {
  N_G  <- rpois(1, Lambda)
  if (N_G == 0) return(list(obs = numeric(0), f_0 = 0))
  lams <- rgamma(N_G, shape = kappa_star, rate = beta_star)
  ns   <- rpois(N_G, lams)
  list(obs = ns[ns >= 1], f_0 = sum(ns == 0))
}


## Table 1: asymptotic ratios across kappa
table_1 <- function() {
  kappas <- c(0.5, 1.0, 2.0, 5.0, 10.0)
  betas  <- kappas / 2
  tab <- t(mapply(function(k, b) {
    p <- asymp_correct(k, b)
    c(kappa = k, a_star = p$a_star,
      bar_a_inf = p$bar_a_inf, rho = p$rho)
  }, kappas, betas))
  rownames(tab) <- NULL
  cat("Table 1: Asymptotic ratios for G_* = Gamma(kappa, kappa/2)\n")
  print(round(tab, 3))
  invisible(tab)
}


## Table 2: dispersion grid at Lambda = 500 
table_2 <- function(n_reps = 1000, Lambda = 500, alpha = 5,
                     seed_base = 4000) {
  kappas <- c(0.5, 1.0, 2.0, 5.0)
  out <- data.frame(
    kappa     = kappas,
    K_n_mean  = NA_real_,
    f0_mean   = NA_real_,
    bar_a     = NA_real_,
    rho_emp   = NA_real_,
    rho_pred  = NA_real_,
    rel_bias  = NA_real_
  )
  for (i in seq_along(kappas)) {
    k <- kappas[i]; b <- k / 2
    set.seed(seed_base + i)
    Kn_v <- f0_v <- bar_v <- hat_v <- numeric(0)
    while (length(Kn_v) < n_reps) {
      d <- simulate_cp(k, b, Lambda)
      if (length(d$obs) < 5) next
      r <- naive_dp(d$obs, alpha = alpha,
                     kappa_0 = k, beta_0 = b)
      Kn_v  <- c(Kn_v, length(d$obs))
      f0_v  <- c(f0_v, d$f_0)
      bar_v <- c(bar_v, r$bar_a)
      hat_v <- c(hat_v, r$hat_f0)
    }
    p <- asymp_correct(k, b)
    out$K_n_mean[i] <- mean(Kn_v)
    out$f0_mean[i]  <- mean(f0_v)
    out$bar_a[i]    <- mean(bar_v)
    out$rho_emp[i]  <- mean(hat_v) / mean(f0_v)
    out$rho_pred[i] <- p$rho
    out$rel_bias[i] <- mean((hat_v - f0_v) / f0_v)
  }
  cat("Table 2: Dispersion grid at Lambda =", Lambda, "\n")
  cat("  Verifies analytical predictions of Table 1.\n\n")
  print(out, row.names = FALSE)
  invisible(out)
}


## Table 3: convergence in Lambda 
table_3 <- function(Lambdas = c(50, 200, 500, 2000, 10000),
                     n_reps = 1000, alpha = 5,
                     kappa_star = 1, beta_star = 0.5,
                     seed_base = 1000) {
  out <- data.frame(
    Lambda = Lambdas, f0_true = NA_real_, bar_a = NA_real_,
    hat_f0 = NA_real_, rel_bias = NA_real_
  )
  for (i in seq_along(Lambdas)) {
    set.seed(seed_base + i)
    f_true <- bar_a <- hat_f <- numeric(0)
    while (length(f_true) < n_reps) {
      d <- simulate_cp(kappa_star, beta_star, Lambdas[i])
      if (length(d$obs) < 5) next
      r <- naive_dp(d$obs, alpha = alpha,
                     kappa_0 = kappa_star, beta_0 = beta_star)
      f_true <- c(f_true, d$f_0)
      bar_a  <- c(bar_a, r$bar_a)
      hat_f  <- c(hat_f, r$hat_f0)
    }
    out$f0_true[i]  <- mean(f_true)
    out$bar_a[i]    <- mean(bar_a)
    out$hat_f0[i]   <- mean(hat_f)
    out$rel_bias[i] <- mean((hat_f - f_true) / f_true)
  }
  cat("Table 3: Convergence in Lambda for G_* = Gamma(1, 0.5)\n")
  cat("  Asymptotic predictions: bar_a_inf = 0.2, rel.bias -> -0.5\n\n")
  print(out, row.names = FALSE)
  invisible(out)
}


## Table 4: misspecified G_0 
table_4 <- function() {
  kappa_star <- 1.0; beta_star <- 0.5
  cases <- list(
    list(label = "Gamma(1, 0.5) (correct)",   k0 = 1.0, b0 = 0.5),
    list(label = "Gamma(1, 1)",               k0 = 1.0, b0 = 1.0),
    list(label = "Gamma(0.5, 0.5)",           k0 = 0.5, b0 = 0.5),
    list(label = "Gamma(0.5, 1)",             k0 = 0.5, b0 = 1.0),
    list(label = "Gamma(2, 1)",               k0 = 2.0, b0 = 1.0),
    list(label = "Gamma(2, 0.5)",             k0 = 2.0, b0 = 0.5),
    list(label = "Gamma(1, 2)",               k0 = 1.0, b0 = 2.0)
  )
  cat("Table 4: Misspecified G_0 with G_* = Gamma(1, 0.5)\n\n")
  cat(sprintf("%-30s %10s %10s\n", "G_0", "bar_a_inf", "rho"))
  for (c in cases) {
    p <- asymp_misspec(kappa_star, beta_star, c$k0, c$b0)
    cat(sprintf("%-30s %10.3f %10.3f\n", c$label,
                p$bar_a_inf, p$rho))
  }
}


## Verify Proposition 4
verify_py <- function(n_reps = 1000, Lambda = 5000,
                       theta = 5, kappa_star = 1, beta_star = 0.5,
                       seed_base = 7000) {
  sigmas <- c(0, 0.25, 0.5, 0.75, 0.9)
  out <- data.frame(
    sigma            = sigmas,
    bar_a_inf_pred   = NA_real_,
    bar_a_emp        = NA_real_,
    rho_pred         = NA_real_,
    rho_emp          = NA_real_
  )
  for (i in seq_along(sigmas)) {
    sg <- sigmas[i]
    p  <- asymp_py(kappa_star, beta_star, sg)
    set.seed(seed_base + i)
    bar_v <- f0_v <- hat_v <- numeric(0)
    while (length(bar_v) < n_reps) {
      d <- simulate_cp(kappa_star, beta_star, Lambda)
      if (length(d$obs) < 5) next
      r <- naive_py(d$obs, theta = theta, sigma = sg,
                     kappa_0 = kappa_star, beta_0 = beta_star)
      bar_v <- c(bar_v, r$bar_a)
      f0_v  <- c(f0_v, d$f_0)
      hat_v <- c(hat_v, r$hat_f0)
    }
    out$bar_a_inf_pred[i] <- p$bar_a_inf
    out$bar_a_emp[i]      <- mean(bar_v)
    out$rho_pred[i]       <- p$rho
    out$rho_emp[i]        <- mean(hat_v) / mean(f0_v)
  }
  cat("Proposition 4 verification (Pitman-Yor)\n")
  cat("  G_* = Gamma(1, 0.5), Lambda =", Lambda,
      ", theta =", theta, "\n")
  cat("  Predicted: bar_a_inf -> sigma * a_star + (1-sigma) * 0.2\n\n")
  print(out, row.names = FALSE)
  invisible(out)
}


## 
if (sys.nframe() == 0) {
  cat("================================================\n\n")
  table_1()
  cat("\n================================================\n\n")
  table_2()
  cat("\n================================================\n\n")
  table_3()
  cat("\n================================================\n\n")
  table_4()
  cat("\n================================================\n\n")
  verify_py()
  cat("\nDone.\n")
}
