# FIGARCH Analysis for Financial Returns
# This script processes financial data, computes log returns, and fits a FIGARCH(1,d,1) model
# using multiple optimization methods (BFGS, CG, L-BFGS-B, EM, and Adam).

# Install and load required packages
# ==============================================================================
# Install necessary packages if not already installed
install.packages(c("readxl", "dplyr", "arfima", "fracdiff"), quietly = TRUE)
library(readxl)
library(dplyr)
library(arfima)
library(fracdiff)

# Load and clean raw Excel data
# ==============================================================================
# Prompt user to select Excel file and load data
load_and_clean_data <- function() {
  file_path <- file.choose()
  returns_raw <- read_excel(file_path)
  cat("First 10 rows of raw data:\n")
  print(head(returns_raw, 10))
  
  # Remove first row if it contains incorrect headers, select first 5 columns
  returns <- returns_raw[-1, 1:5]
  
  # Rename columns for clarity
  colnames(returns) <- c("Date", "XAU", "RHOD_LON", "XPT", "XPD")
  
  # Convert 'Date' column to Date format
  if (is.numeric(returns$Date)) {
    returns$Date <- as.Date(as.numeric(returns$Date), origin = "1899-12-30")
  } else {
    returns$Date <- as.Date(returns$Date, format = "%d/%m/%Y")
  }
  
  # Remove rows with missing values and sort by date
  returns <- returns %>%
    filter(!is.na(Date) & !is.na(RHOD_LON) & !is.na(XAU) & !is.na(XPT) & !is.na(XPD)) %>%
    arrange(Date)
  row.names(returns) <- NULL
  
  cat("Number of observations after cleaning:", nrow(returns), "\n")
  return(returns)
}

# Compute log returns for selected assets
# ==============================================================================
compute_log_returns <- function(returns) {
  cols_to_logdiff <- c("XAU", "RHOD_LON", "XPT", "XPD")
  log_returns <- data.frame(Date = returns$Date[-1])  # One less row due to differencing
  
  for (col in cols_to_logdiff) {
    returns[[col]] <- as.numeric(returns[[col]])
    log_returns[[col]] <- diff(log(returns[[col]]), lag = 1)
  }
  
  cat("Number of log returns computed:", nrow(log_returns), "\n")
  cat("First 10 rows of log returns:\n")
  print(head(log_returns, 10))
  
  return(log_returns)
}

# Compute pi coefficients for FIGARCH model
# ==============================================================================
compute_pi <- function(d, J) {
  pi_j <- numeric(J)
  pi_j[1] <- 1
  for (j in 2:J) {
    pi_j[j] <- pi_j[j - 1] * (j - 1 - d) / j
  }
  return(pi_j)
}

# Compute FIGARCH variance
# ==============================================================================
figarch_variance <- function(r, theta, J) {
  mu <- theta[1]
  omega <- theta[2]
  phi1 <- theta[3]
  beta1 <- theta[4]
  d <- theta[5]
  
  T <- length(r)
  e <- r - mu
  e2 <- e^2
  
  sigma2 <- numeric(T)
  sigma2[1] <- var(r)
  
  pi_j <- compute_pi(d, J)
  pi_jm1 <- c(0, pi_j[-length(pi_j)])
  lambda_j <- pi_j - (phi1 + beta1) * pi_jm1
  
  omega_star <- omega / (1 - beta1)
  
  for (t in 2:T) {
    start_idx <- max(1, t - J)
    e2_history <- rev(e2[(start_idx):(t - 1)])
    pad_size <- J - length(e2_history)
    
    if (pad_size > 0) {
      arch_term <- sum(lambda_j[(pad_size + 1):J] * e2_history)
      avg_e2 <- mean(e2[1:min(100, T)])
      arch_term <- arch_term + sum(lambda_j[1:pad_size]) * avg_e2
    } else {
      arch_term <- sum(lambda_j[1:J] * rev(e2[(t - J):(t - 1)]))
    }
    
    sigma2[t] <- omega_star + arch_term + beta1 * (sigma2[t - 1] - omega_star)
    
    if (is.na(sigma2[t]) || sigma2[t] <= 1e-6) {
      sigma2[t] <- 1e-6
    }
  }
  
  return(sigma2)
}

# Compute FIGARCH log-likelihood
# ==============================================================================
loglik_figarch <- function(theta, r, J) {
  mu <- theta[1]
  omega <- theta[2]
  phi1 <- theta[3]
  beta1 <- theta[4]
  d <- theta[5]
  
  # Parameter constraints
  if (omega <= 0 || phi1 < 0 || beta1 < 0 || d <= 0 || d >= 1) {
    return(1e10)
  }
  
  tryCatch({
    sigma2 <- figarch_variance(r, theta, J)
    
    z <- (r - mu) / sqrt(sigma2)
    valid_idx <- !is.na(z) & is.finite(z)
    loglik <- -0.5 * log(2 * pi) - 0.5 * log(sigma2[valid_idx]) - 0.5 * z[valid_idx]^2
    
    negLL <- -sum(loglik)
    
    if (is.na(negLL) || is.infinite(negLL)) {
      negLL <- 1e10
    }
    return(negLL)
  }, error = function(e) {
    return(1e10)
  })
}

# Compute single time step log-likelihood
# ==============================================================================
single_loglik <- function(theta, r, J, t) {
  sigma2 <- figarch_variance(r, theta, J)
  mu <- theta[1]
  
  if (t > length(r) || is.na(sigma2[t]) || sigma2[t] <= 1e-10) {
    ll <- 0
  } else {
    z <- (r[t] - mu) / sqrt(sigma2[t])
    ll <- -0.5 * (log(2 * pi * sigma2[t]) + z^2)
  }
  return(ll)
}

# Transformed FIGARCH log-likelihood
# ==============================================================================
loglik_figarch_trans <- function(theta_trans, r, J) {
  mu <- theta_trans[1]
  omega <- exp(theta_trans[2])  # Transform log(omega)
  phi1 <- theta_trans[3]
  beta1 <- theta_trans[4]
  d <- theta_trans[5]
  theta <- c(mu, omega, phi1, beta1, d)
  return(loglik_figarch(theta, r, J))
}

# Initialize parameters for FIGARCH model
# ==============================================================================
initialize_parameters <- function(r) {
  mu_init <- mean(r)
  eps2 <- (r - mu_init)^2
  y <- eps2[3:length(eps2)]
  x1 <- eps2[2:(length(eps2) - 1)]
  x2 <- eps2[1:(length(eps2) - 2)]
  X <- cbind(1, x1, x2)
  ols_model <- lm(y ~ x1 + x2)
  omega_init <- max(abs(coef(ols_model)[1]), 1e-6)
  phi1_init <- max(coef(ols_model)[2], 0)
  beta1_init <- max(abs(coef(ols_model)[3]), 0)
  d_init <- 0.3
  d_init <- min(max(d_init, 0.01), 0.99)
  return(c(mu_init, log(omega_init), phi1_init, beta1_init, d_init))
}

# Fit FIGARCH model using multiple optimization methods
# ==============================================================================
fit_figarch_models <- function(r, theta_init_trans, J = 100) {
  # Test initial log-likelihood
  theta_init <- c(theta_init_trans[1], exp(theta_init_trans[2]), 
                  theta_init_trans[3], theta_init_trans[4], theta_init_trans[5])
  negLL <- loglik_figarch(theta_init, r, J)
  cat("Initial negative log-likelihood:", negLL, "\n")
  
  # BFGS Optimization
  fit_bfgs <- optim(
    par = theta_init_trans,
    fn = loglik_figarch_trans,
    r = r,
    J = J,
    method = "BFGS",
    control = list(maxit = 1000, trace = 1)
  )
  
  theta_hat <- c(
    mu = fit_bfgs$par[1],
    omega = exp(fit_bfgs$par[2]),
    phi1 = fit_bfgs$par[3],
    beta1 = fit_bfgs$par[4],
    d = fit_bfgs$par[5]
  )
  
  loglik_bfgs <- -fit_bfgs$value
  n <- length(r)
  k <- length(theta_hat)
  AIC_bfgs <- -2 * loglik_bfgs + 2 * k
  BIC_bfgs <- -2 * loglik_bfgs + log(n) * k
  
  cat("\n--- BFGS (Transformed) Results ---\n")
  cat("Estimated parameters:\n")
  print(theta_hat)
  cat("Log-likelihood:", loglik_bfgs, "\n")
  cat("AIC:", AIC_bfgs, "\n")
  cat("BIC:", BIC_bfgs, "\n")
  cat("Number of iterations:", fit_bfgs$counts[1], "\n")
  
  # Conjugate Gradient (CG) Optimization
  fit_cg <- optim(
    par = theta_init_trans,
    fn = loglik_figarch_trans,
    r = r,
    J = J,
    method = "CG",
    control = list(maxit = 50000, trace = 1)
  )
  
  theta_hat <- c(
    mu = fit_cg$par[1],
    omega = exp(fit_cg$par[2]),
    phi1 = fit_cg$par[3],
    beta1 = fit_cg$par[4],
    d = fit_cg$par[5]
  )
  
  loglik_cg <- -fit_cg$value
  AIC_cg <- -2 * loglik_cg + 2 * k
  BIC_cg <- -2 * loglik_cg + log(n) * k
  
  cat("\n--- CG (Transformed) Results ---\n")
  cat("Estimated parameters:\n")
  print(theta_hat)
  cat("Log-likelihood:", loglik_cg, "\n")
  cat("AIC:", AIC_cg, "\n")
  cat("BIC:", BIC_cg, "\n")
  cat("Number of iterations:", fit_cg$counts[1], "\n")
  
  # L-BFGS-B Optimization
  fit_lbfgs <- optim(
    par = theta_init_trans,
    fn = loglik_figarch_trans,
    r = r,
    J = J,
    method = "L-BFGS-B",
    control = list(maxit = 500, trace = 1)
  )
  
  theta_hat <- c(
    mu = fit_lbfgs$par[1],
    omega = exp(fit_lbfgs$par[2]),
    phi1 = fit_lbfgs$par[3],
    beta1 = fit_lbfgs$par[4],
    d = fit_lbfgs$par[5]
  )
  
  loglik_lbfgs <- -fit_lbfgs$value
  AIC_lbfgs <- -2 * loglik_lbfgs + 2 * k
  BIC_lbfgs <- -2 * loglik_lbfgs + log(n) * k
  
  cat("\n--- L-BFGS-B (Transformed) Results ---\n")
  cat("Estimated parameters:\n")
  print(theta_hat)
  cat("Log-likelihood:", loglik_lbfgs, "\n")
  cat("AIC:", AIC_lbfgs, "\n")
  cat("BIC:", BIC_lbfgs, "\n")
  cat("Number of iterations:", fit_lbfgs$counts[1], "\n")
}

# EM Algorithm for FIGARCH
# ==============================================================================
em_figarch_trans <- function(r, J = 100, theta0_trans, max_iter = 100, tol = 1e-8) {
  theta <- theta0_trans  # [mu, log(omega), phi1, beta1, d]
  iter <- 0
  loglik_old <- -Inf
  
  repeat {
    loglik <- -loglik_figarch_trans(theta, r, J)
    if (abs(loglik - loglik_old) < tol || iter >= max_iter) break
    loglik_old <- loglik
    iter <- iter + 1
    
    opt_result <- optim(
      par = theta,
      fn = function(params) loglik_figarch_trans(params, r, J),
      method = "BFGS",
      control = list(maxit = 500)
    )
    
    theta <- opt_result$par
  }
  
  return(list(
    par = theta,
    value = -loglik,
    iterations = iter
  ))
}

# Adam Optimization for FIGARCH
# ==============================================================================
adam_figarch_trans <- function(r, J = 100, theta0, max_iter = 50000, alpha = 0.01, 
                               beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8, tol = 1e-4, trace = 1) {
  theta <- theta0  # [mu, log(omega), phi1, beta1, d]
  m <- rep(0, length(theta))
  v <- rep(0, length(theta))
  t_iter <- 0
  loglik_old <- -Inf
  
  compute_grad <- function(theta, r, J) {
    grad <- numeric(length(theta))
    eps <- 1e-6
    for (i in 1:length(theta)) {
      theta_plus <- theta_minus <- theta
      theta_plus[i] <- theta[i] + eps
      theta_minus[i] <- theta[i] - eps
      grad[i] <- (loglik_figarch_trans(theta_plus, r, J) - 
                    loglik_figarch_trans(theta_minus, r, J)) / (2 * eps)
    }
    return(grad)
  }
  
  repeat {
    t_iter <- t_iter + 1
    loglik <- -loglik_figarch_trans(theta, r, J)
    grad <- compute_grad(theta, r, J)
    
    m <- beta1 * m + (1 - beta1) * grad
    v <- beta2 * v + (1 - beta2) * (grad^2)
    
    m_hat <- m / (1 - beta1^t_iter)
    v_hat <- v / (1 - beta2^t_iter)
    
    theta <- theta - alpha * m_hat / (sqrt(v_hat) + epsilon)
    
    # Apply constraints
    theta[2] <- max(log(1e-6), min(log(10), theta[2]))  # log(omega)
    theta[3] <- max(0, min(0.999, theta[3]))           # phi1
    theta[4] <- max(0, min(0.999, theta[4]))           # beta1
    theta[5] <- max(1e-6, min(0.999, theta[5]))        # d
    
    if (abs(loglik - loglik_old) < tol || t_iter >= max_iter) break
    loglik_old <- loglik
  }
  
  return(list(par = theta, value = loglik, iterations = t_iter))
}

# Main execution
# ==============================================================================
main <- function() {
  # Load and clean data
  returns <- load_and_clean_data()
  
  # Compute log returns
  log_returns <- compute_log_returns(returns)
  
  # Use XAU for modeling
  r <- log_returns$XAU
  
  # Initialize parameters
  theta_init_trans <- initialize_parameters(r)
  cat("Initial transformed parameters:", theta_init_trans, "\n")
  
  # Fit FIGARCH models
  J <- 100
  fit_figarch_models(r, theta_init_trans, J)
  
  # EM Optimization
  fit_em <- em_figarch_trans(r = r, J = J, theta0_trans = theta_init_trans)
  
  theta_hat <- c(
    mu = fit_em$par[1],
    omega = exp(fit_em$par[2]),
    phi1 = fit_em$par[3],
    beta1 = fit_em$par[4],
    d = fit_em$par[5]
  )
  
  loglik_em <- -fit_em$value
  n <- length(r)
  k <- length(theta_hat)
  AIC_em <- -2 * loglik_em + 2 * k
  BIC_em <- -2 * loglik_em + log(n) * k
  
  cat("\n--- EM (Transformed) Results ---\n")
  cat("Estimated parameters:\n")
  print(theta_hat)
  cat("Log-likelihood:", loglik_em, "\n")
  cat("AIC:", AIC_em, "\n")
  cat("BIC:", BIC_em, "\n")
  cat("Number of iterations:", fit_em$iterations, "\n")
  
  # Adam Optimization
  fit_adam <- adam_figarch_trans(r = r, J = J, theta0 = theta_init_trans)
  
  theta_hat <- c(
    mu = fit_adam$par[1],
    omega = exp(fit_adam$par[2]),
    phi1 = fit_adam$par[3],
    beta1 = fit_adam$par[4],
    d = fit_adam$par[5]
  )
  
  loglik_adam <- fit_adam$value
  AIC_adam <- -2 * loglik_adam + 2 * k
  BIC_adam <- -2 * loglik_adam + log(n) * k
  
  cat("\n--- Adam (Transformed) Results ---\n")
  cat("Estimated parameters:\n")
  print(theta_hat)
  cat("Log-likelihood:", loglik_adam, "\n")
  cat("AIC:", AIC_adam, "\n")
  cat("BIC:", BIC_adam, "\n")
  cat("Number of iterations:", fit_adam$iterations, "\n")
}

# Run the analysis
main()