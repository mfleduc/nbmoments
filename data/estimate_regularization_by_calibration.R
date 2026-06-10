
library(pracma)

# ---------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------

n_sims <- 40

sample_sizes <- c(10, 20, 40)

ratio_means <- c(0.5, 1, 5)

dispersion_vals <- c(1/10)

moment_sets <- list(
  c(1, 2),
  c(-1, 1),
  c(-2, -1, 1),
  c(-1, 1, 2),
  c(-2, -1, 1, 2)
)

moment_names <- c(
  "12",
  "m11",
  "m211",
  "m112",
  "m2112"
)

#
# log10 spacing by 0.5
#
ridge_grid <- 10^c(seq(-6, -2, by = 0.33333),-2)

coverage_probs <- c(0.05,seq(0.1, 0.9, by=0.2), 0.95)

save_file <- "calibration_results.rds"

# ---------------------------------------------------------
# CvM statistic
# ---------------------------------------------------------

cvm_calibration <- function(pit_values) {
  
  pit_values <- pit_values[is.finite(pit_values)]
  
  n <- length(pit_values)
  
  if(n < 2)
    return(NA_real_)
  
  eps <- 1e-12
  
  u <- sort(
    pmin(
      pmax(pit_values, eps),
      1 - eps
    )
  )
  
  i <- seq_len(n)
  
  1 / (12 * n^2) +
    1/n*sum(
      (
        u - (2 * i - 1) / (2 * n)
      )^2
    )
}

# ---------------------------------------------------------
# PIT helper
# ---------------------------------------------------------

pit_value <- function(fit, theta_true) {
  
  out <- cdf.maxent_fit(
    fit,
    grid = theta_true
  )
  
  as.numeric(out$cdf[1])
}

# ---------------------------------------------------------
# Safe fitting wrapper
# ---------------------------------------------------------

safe_fit <- function(...) {
  
  tryCatch(
    estimate_dist(...),
    error = function(e) NULL,
    warning = function(w) invokeRestart("muffleWarning")
  )
}

# ---------------------------------------------------------
# Main storage
# ---------------------------------------------------------

results <- list()

counter <- 1

# =========================================================
# MAIN LOOP
# =========================================================

for(N in sample_sizes) {
  
  for(mu_ratio in ratio_means) {
    
    for(phi_num in dispersion_vals) {
      
      for(phi_den in dispersion_vals) {
        
        for(mm in seq_along(moment_sets)) {
          
          moments <- moment_sets[[mm]]
          moment_name <- moment_names[mm]
          
          for(ridge_pen in ridge_grid) {
            
            cat(
              "\n====================================\n",
              sprintf(
                "N=%d  mu=%.2f  phiN=%g  phiD=%g  moments=%s  ridge=%g\n",
                N,
                mu_ratio,
                phi_num,
                phi_den,
                moment_name,
                ridge_pen
              ),
              "====================================\n"
            )
            
            pit_vals <- numeric(n_sims)
            
            coverage_mat <- matrix(
              0,
              nrow = n_sims,
              ncol = length(coverage_probs)
            )
            
            widths <- matrix(
              NA_real_,
              nrow = n_sims,
              ncol = length(coverage_probs)
            )
            
            converged <- logical(n_sims)
            
            failures <- 0
            moment_error_mat <- matrix(
              NA_real_,
              nrow = n_sims,
              ncol = length(moments)
            )
            max_moment_error <- numeric(n_sims)
            for(sim in seq_len(n_sims)) {
              
              # -----------------------------------
              # Simulate data
              # -----------------------------------
              
              #
              # denominator mean fixed at 1
              #
              mu_den <- 10
              
              #
              # numerator mean determines ratio
              #
              mu_num <- mu_ratio*mu_den
              
              #
              # Negative binomial parameterization:
              # Var(X)=mu + phi*mu^2
              #
              size_num <- 1 / phi_num
              size_den <- 1 / phi_den
              
              x_num <- rnbinom(
                N,
                mu = mu_num,
                size = size_num
              )
              
              x_den <- rnbinom(
                N,
                mu = mu_den,
                size = size_den
              )
              
              # -----------------------------------
              # Fit MaxEnt distribution
              # -----------------------------------
              
              fit <- estimate_dist(
                x_num,
                x_den,
                moments = moments,
                ridge_pen = ridge_pen,
                ellmax = 120,
                n_nodes = 150,
                n_outer = 150
              )
              
              if(is.null(fit)) {
                
                failures <- failures + 1
                next
              }
              moment_error_mat[sim, ] <- fit$moment_errors
              
              max_moment_error[sim] <- max(
                fit$moment_errors,
                na.rm = TRUE
              )
             
              
              max_moment_error <- numeric(n_sims)
              
              converged[sim] <- isTRUE(fit$converged)
              
              # -----------------------------------
              # PIT
              # -----------------------------------
              
              pit_vals[sim] <- pit_value(
                fit,
                mu_ratio
              )
              
              # -----------------------------------
              # Coverage
              # -----------------------------------
              
              for(pp in seq_along(coverage_probs)) {
                
                ci = credible_interval.maxent_fit(
                  fit,
                  prob = coverage_probs[pp]
                )
                
                coverage_mat[sim, pp] <-
                  as.integer(
                    mu_ratio >= ci$lower &&
                      mu_ratio <= ci$upper
                  )
                
                widths[sim, pp] <-
                  ci$upper - ci$lower
              }
            }
            
            # -----------------------------------
            # Aggregate metrics
            # -----------------------------------
            mean_moment_error <- colMeans(
              moment_error_mat,
              na.rm = TRUE
            )
            
            median_max_moment_error <- median(
              max_moment_error,
              na.rm = TRUE
            )
            
            mean_max_moment_error <- mean(
              max_moment_error,
              na.rm = TRUE
            )
            
            empirical_coverage <-
              colMeans(
                coverage_mat,
                na.rm = TRUE
              )
            
            avg_width <-
              colMeans(
                widths,
                na.rm = TRUE
              )
            
            cvm <- cvm_calibration(pit_vals)
            
            pit_mean <- mean(
              pit_vals,
              na.rm = TRUE
            )
            
            pit_var <- var(
              pit_vals,
              na.rm = TRUE
            )
            
            # -----------------------------------
            # Store result
            # -----------------------------------
            
            results[[counter]] <- list(
              N = N,
              mu_ratio = mu_ratio,
              phi_num = phi_num,
              phi_den = phi_den,
              moments = moments,
              moment_name = moment_name,
              ridge_pen = ridge_pen,
              
              empirical_coverage = empirical_coverage,
              avg_width = avg_width,
              
              cvm = cvm,
              pit_mean = pit_mean,
              pit_var = pit_var,
              
              converged_rate = mean(converged),
              failures = failures,
              
              #
              # NEW
              #
              mean_moment_error = mean_moment_error,
              mean_max_moment_error = mean_max_moment_error,
              median_max_moment_error = median_max_moment_error,
              
              pit_vals = pit_vals
            )
            counter <- counter + 1
            
            # -----------------------------------
            # Save intermediate results
            # -----------------------------------
            
            saveRDS(
              results,
              file = save_file
            )
            
            cat(
              sprintf(
                "DONE | CvM=%.5f | Conv=%.2f | Fail=%d\n",
                cvm,
                mean(converged),
                failures
              )
            )
          }
        }
      }
    }
  }
}

# =========================================================
# Convert to data.frame summary
# =========================================================

summary_df <- do.call(
  rbind,
  lapply(results, function(x) {
    
    data.frame(
      N = x$N,
      mu_ratio = x$mu_ratio,
      phi_num = x$phi_num,
      phi_den = x$phi_den,
      moment_name = x$moment_name,
      ridge_pen = x$ridge_pen,
      
      cvm = x$cvm,
      pit_mean = x$pit_mean,
      pit_var = x$pit_var,
      
      converged_rate = x$converged_rate,
      failures = x$failures,
      
      #
      # NEW
      #
      mean_max_moment_error =
        x$mean_max_moment_error,
      
      median_max_moment_error =
        x$median_max_moment_error,
      
      #
      # Per-moment reconstruction errors
      #
      t(as.data.frame(
        setNames(
          as.list(x$mean_moment_error),
          paste0(
            "moment_err_",
            seq_along(x$mean_moment_error)
          )
        )
      )),
      
      cov_05 = x$empirical_coverage[1],
      cov_25 = x$empirical_coverage[2],
      cov_50 = x$empirical_coverage[3],
      cov_75 = x$empirical_coverage[4],
      cov_95 = x$empirical_coverage[5],
      
      width_05 = x$avg_width[1],
      width_25 = x$avg_width[2],
      width_50 = x$avg_width[3],
      width_75 = x$avg_width[4],
      width_95 = x$avg_width[5]
    )
  })
)
saveRDS(summary_df, "calibration_summary_df.rds")
write.csv(
  summary_df,
  "calibration_summary_df.csv",
  row.names = FALSE
)
cat("\nFinished.\n")