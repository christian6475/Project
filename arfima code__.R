# ARFIMA(1,d,1) Analysis for Financial Returns
# This script processes financial data, computes log returns, and fits an ARFIMA(1,d,1) model
# using multiple optimization methods (BFGS, CG, L-BFGS-B, EM, and Adam).

# Install and load required packages
# ==============================================================================
install.packages(c("readxl", "dplyr", "arfima", "fracdiff"), quietly = TRUE)
library(readxl)
library(dplyr)
library(arfima)
library(fracdiff)

# Load and clean raw Excel data
# ==============================================================================
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
  log_returns <- data.frame(Date = returns$Date[-1])  # One less row due to diff()
  
  for (col in cols_to_logdiff) {
    returns[[col]] <- as.numeric(returns[[col]])
    log_returns[[col]] <- diff(log(returns[[col]]), lag = 1)
  }
  
  cat("Number of log returns computed:", nrow(log_returns), "\n")
  cat("First 10 rows of log returns:\n")
  print(head(log_returns, 10))
  
  # Plot original prices and log returns for XAU
  plot(returns$Date, returns$XAU, type = "l", xlab = "Date", ylab = "XAU", main = "Gold Price (XAU)")
  plot(log_returns$Date, log_returns$XAU, type = "l", xlab = "Date", ylab = "Log Return", main = "Log Return of XAU")
  
  return(log_returns)
}

# Define log-likelihood function for ARFIMA(1,d,1)
# ==============================================================================
log_likelihood_arfima <- function(params, data) {
  phi <- params[1]
  theta <- params[2]
  d <- params[3]
  
  # Enforce parameter constraints
  if (abs(phi) >= 1 || abs(theta) >= 1 || d <= -1 || d >= 1) {
    return(-Inf)  # Penalize invalid parameters
  }
  
  tryCatch({
    fit <- arfima(data, order = c(1, 0, 1), fixed = list(phi = phi, theta = theta, frac = d), 
                  lmodel = "d", quiet = TRUE)
    return(-logLik(fit)[1])
  }, error = function(e) {
    return(-Inf)  # Return -Inf on failure
  })
}

# Initialize parameters for ARFIMA(1,d,1)
# ==============================================================================
initialize_parameters <- function(data) {
  # Estimate AR(1) parameter using OLS
  Y <- data[2:length(data)]
  X <- data[1:(length(data) - 1)]
  ols_fit <- lm(Y ~ X)
  phi_init <- coef(ols_fit)[2]
  
  # Initialize MA(1) parameter
  theta_init <- 0
  
  # Estimate fractional differencing parameter
  d_init <- fracdiff(data, nar = 0, nma = 0)$d
  
  return(c(phi = phi_init, theta = theta_init, d = d_init))
}

# Fit ARFIMA model using multiple optimization methods
# ==============================================================================
fit_arfima_models <- function(data, initial_params) {
  # Default ARFIMA fit
  model <- arfima(data, order = c(1, 0, 1))
  cat("Default ARFIMA fit:\n")
  print(summary(model))
  
  # BFGS Optimization
  temps_bfgs <- system.time({
    optim_bfgs <- optim(par = initial_params, 
                        fn = log_likelihood_arfima, 
                        data = data, 
                        method = "BFGS",
                        control = list(maxit = 500, trace = 1))
  })
  
  cat("\nBFGS optimization results:\n",
      "phi:", optim_bfgs$par[1], "\n",
      "theta:", optim_bfgs$par[2], "\n",
      "d:", optim_bfgs$par[3], "\n",
      "Log-likelihood:", -optim_bfgs$value, "\n")
  cat("BFGS optimization time (seconds):\n")
  print(temps_bfgs)
  print(optim_bfgs$counts)
  
  # Calculate AIC and BIC
  loglik <- -optim_bfgs$value
  n <- length(data)
  k <- length(optim_bfgs$par)
  AIC <- -2 * loglik + 2 * k
  BIC <- -2 * loglik + log(n) * k
  cat("AIC:", AIC, "\nBIC:", BIC, "\n")
  
  # Conjugate Gradient (CG) Optimization
  counter <<- 0
  log_likelihood_arfima_counter <- function(params, data) {
    counter <<- counter + 1
    log_likelihood_arfima(params, data)
  }
  
  temps_cg <- system.time({
    counter <<- 0
    optim_cg <- optim(par = initial_params, 
                      fn = log_likelihood_arfima_counter, 
                      data = data, 
                      method = "CG",
                      control = list(maxit = 800, trace = 1))
  })
  
  cat("\nCG optimization results:\n",
      "phi:", optim_cg$par[1], "\n",
      "theta:", optim_cg$par[2], "\n",
      "d:", optim_cg$par[3], "\n",
      "Log-likelihood:", -optim_cg$value, "\n")
  cat("CG optimization time (seconds):\n")
  print(temps_cg)
  print(optim_cg$counts)
  
  loglik <- -optim_cg$value
  AIC <- -2 * loglik + 2 * k
  BIC <- -2 * loglik + log(n) * k
  cat("AIC:", AIC, "\nBIC:", BIC, "\n")
  
  # L-BFGS-B Optimization
  temps_lbfgs <- system.time({
    optim_lbfgs <- optim(par = initial_params, 
                         fn = log_likelihood_arfima, 
                         data = data, 
                         method = "L-BFGS-B",
                         control = list(maxit = 500, trace = 1))
  })
  
  cat("\nL-BFGS optimization results:\n",
      "phi:", optim_lbfgs$par[1], "\n",
      "theta:", optim_lbfgs$par[2], "\n",
      "d:", optim_lbfgs$par[3], "\n",
      "Log-likelihood:", -optim_lbfgs$value, "\n")
  cat("L-BFGS optimization time (seconds):\n")
  print(temps_lbfgs)
  print(optim_lbfgs$counts)
  
  loglik <- -optim_lbfgs$value
  AIC <- -2 * loglik + 2 * k
  BIC <- -2 * loglik + log(n) * k
  cat("AIC:", AIC, "\nBIC:", BIC, "\n")
}

# EM Algorithm for ARFIMA(1,d,1)
# ==============================================================================
em_arfima <- function(data, initial_params, max_iter = 100, tol = 1e-6) {
  params <- initial_params
  phi <- params[1]
  theta <- params[2]
  d <- params[3]
  
  loglik_old <- -Inf
  converged <- FALSE
  
  for (iter in 1:max_iter) {
    # E-step: Calculate current log-likelihood
    current_loglik <- -log_likelihood_arfima(params, data)
    
    # Check convergence
    if (abs(current_loglik - loglik_old) < tol && iter > 1) {
      converged <- TRUE
      break
    }
    
    loglik_old <- current_loglik
    
    # M-step: Maximize log-likelihood
    optim_result <- optim(par = params, 
                          fn = log_likelihood_arfima, 
                          data = data, 
                          method = "BFGS",
                          control = list(maxit = 100))
    
    params <- optim_result$par
    phi <- params[1]
    theta <- params[2]
    d <- params[3]
    
    cat(sprintf("Iteration %d: phi=%.4f, theta=%.4f, d=%.4f, log-likelihood=%.4f\n",
                iter, phi, theta, d, -optim_result$value))
  }
  
  if (converged) {
    cat("Convergence reached after", iter, "iterations.\n")
  } else {
    cat("No convergence after", max_iter, "iterations.\n")
  }
  
  return(list(phi = phi, theta = theta, d = d, 
              log_likelihood = -optim_result$value, 
              iterations = iter, 
              converged = converged))
}

# Adam Optimization for ARFIMA(1,d,1)
# ==============================================================================
adam_arfima <- function(data, initial_params, max_iter = 1000, alpha = 0.1, 
                        beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8, tol = 1e-6) {
  theta <- initial_params
  m <- rep(0, length(theta))  # First moment
  v <- rep(0, length(theta))  # Second moment
  t <- 0  # Iteration counter
  loglik_old <- -Inf
  
  compute_grad <- function(theta, data) {
    grad <- numeric(length(theta))
    eps <- 1e-6
    for (i in 1:length(theta)) {
      theta_plus <- theta_minus <- theta
      theta_plus[i] <- theta[i] + eps
      theta_minus[i] <- theta[i] - eps
      grad[i] <- (log_likelihood_arfima(theta_plus, data) - 
                    log_likelihood_arfima(theta_minus, data)) / (2 * eps)
    }
    return(grad)
  }
  
  for (t in 1:max_iter) {
    loglik <- log_likelihood_arfima(theta, data)
    if (is.na(loglik) || is.infinite(loglik)) {
      cat("Warning: log-likelihood is NA or infinite at iteration", t, ". Stopping.\n")
      break
    }
    
    grad <- compute_grad(theta, data)
    if (any(is.na(grad) | is.infinite(grad))) {
      cat("Warning: Gradient contains NA or infinite values at iteration", t, ". Stopping.\n")
      break
    }
    
    # Update moments
    m <- beta1 * m + (1 - beta1) * grad
    v <- beta2 * v + (1 - beta2) * (grad^2)
    
    # Bias correction
    m_hat <- m / (1 - beta1^t)
    v_hat <- v / (1 - beta2^t)
    
    # Update parameters
    theta_new <- theta - alpha * m_hat / (sqrt(v_hat) + epsilon)
    
    # Apply constraints
    theta_new[1] <- max(-0.999, min(0.999, theta_new[1]))  # phi
    theta_new[2] <- max(-0.999, min(0.999, theta_new[2]))  # theta
    theta_new[3] <- max(-0.999, min(0.999, theta_new[3]))  # d
    
    theta <- theta_new
    
    # Check convergence
    if (is.finite(loglik) && is.finite(loglik_old) && abs(loglik - loglik_old) < tol) break
    loglik_old <- loglik
  }
  
  return(list(theta = theta, loglik = -loglik, iterations = t))
}

# Main execution
# ==============================================================================
main <- function() {
  # Load and clean data
  returns <- load_and_clean_data()
  
  # Compute log returns
  log_returns <- compute_log_returns(returns)
  
  # Use RHOD_LON for modeling
  r <- log_returns$RHOD_LON
  
  # Initialize parameters
  initial_params <- initialize_parameters(r)
  cat("Initial parameters:", initial_params, "\n")
  
  # Fit ARFIMA models
  fit_arfima_models(r, initial_params)
  
  # EM algorithm
  temps_em <- system.time({
    em_result <- em_arfima(r, initial_params)
  })
  cat("\nEM algorithm results:\n",
      "phi:", em_result$phi, "\n",
      "theta:", em_result$theta, "\n",
      "d:", em_result$d, "\n",
      "Log-likelihood:", em_result$log_likelihood, "\n",
      "Converged:", em_result$converged, "\n",
      "Number of iterations:", em_result$iterations, "\n")
  cat("EM optimization time (seconds):\n")
  print(temps_em)
  
  # Calculate AIC and BIC for EM
  loglik <- em_result$log_likelihood
  n <- length(r)
  k <- 3
  AIC <- -2 * loglik + 2 * k
  BIC <- -2 * loglik + log(n) * k
  cat("AIC:", AIC, "\nBIC:", BIC, "\n")
  
  # Adam optimization
  temps_adam <- system.time({
    adam_result <- adam_arfima(r, initial_params)
  })
  cat("\nAdam optimization results:\n",
      "phi:", adam_result$theta[1], "\n",
      "theta:", adam_result$theta[2], "\n",
      "d:", adam_result$theta[3], "\n",
      "Log-likelihood:", adam_result$loglik, "\n",
      "Iterations:", adam_result$iterations, "\n")
  cat("Adam optimization time (seconds):\n")
  print(temps_adam)
}

main()