dual_objective_hermite =function(lambda,moments,powers, quad, ridge_pen = 1e-5) {
 lambda = transform_lambda(lambda, powers)
 
 pos_idx <- which.max(powers)
 neg_idx <- which.min(powers)
 # max_log_w <- max(quad$log_w)
 # disp(quad)
 # w <- quad$weights# * exp(quad$log_w - max_log_w)
 # Z <- sum(w)
 max_lp <- max(quad$log_integrand)
 
 w_unnorm <-
   quad$weights *
   exp(quad$log_integrand - max_lp)
 
 Z <- sum(w_unnorm)
 # browser()
 p = w_unnorm / Z
 
 # E_q[phi * exp(eta)] / Z
 # print(dim(p))
 grad <- colSums(quad$ExpBasis * p)
 value <- (max_lp + log(Z)) -
    sum(lambda * moments) +
    ridge_pen * sum(lambda^2)
 # value <- log(Z) - sum(lambda * moments) + ridge_pen * sum(lambda^2)
 
 gradient <- grad - moments + 2 * ridge_pen * lambda
 
 # list(value = value, gradient = gradient)
  gradient[pos_idx] <-
    gradient[pos_idx] * lambda[pos_idx]
  
  if(powers[neg_idx] < 0) {
    gradient[neg_idx] <-
      gradient[neg_idx] * lambda[neg_idx]
  }
  # gradient=gradient+2 * ridge_pen * lambda
  list(
    value = value,
    gradient = gradient
  )
}
maxent_distribution_logspace <- function(
    moments,
    values,
    n_nodes = 100,
    n_outer = 100,
    optim_tol = 1e-5,
    ridge_pen = 1e-6,
    verbose = FALSE
) {
  powers=moments
  #
  # Basic checks
  #
  lambda_tol =optim_tol
  stopifnot(length(moments) == length(values))
  stopifnot(all(is.finite(values)))
  stopifnot(all(values > 0))
  scale_x <- values[which(moments == 1)]
  
  values = values / (scale_x ^ moments)
  #
  # Sort moments
  #
  
  ord <- order(moments)
  
  moments <- moments[ord]
  values  <- values[ord]
  
  m <- length(moments)
  
  #
  # Initial values
  #
  # Use coercive tails automatically
  #
  
  theta <- rep(0, m)
  
  # pos_idx <- which.max(moments)
  # neg_idx <- which.min(moments)
  
  # if(moments[pos_idx] > 0)
  #   lambda[pos_idx] <- -0.1
  # 
  # if(moments[neg_idx] < 0)
  #   lambda[neg_idx] <- -0.1
  
  #
  # Mean initialization
  #
  
  if(any(moments == 1)) {
    
    idx1 <- which(moments == 1)
    
    theta[idx1] <-
      theta[idx1] - 1 / values[idx1]
  }
  
  #
  # Main optimization loop
  #
  
  converged <- FALSE
  
  for(outer in seq_len(n_outer)) {
    
    #
    # Build adapted quadrature
    #
    lambda =transform_lambda(theta, powers)
    quad <- adapted_hermite_quad(
      lambda = lambda,
      powers = moments,
      n_nodes = n_nodes
    )
  
    # Optimize
    #
    
    opt <- optim(
      par = theta,
      
      fn = function(par) {
        lambda_local <- transform_lambda(par, powers)

        quad_local <- adapted_hermite_quad(
          lambda = lambda_local,
          powers = moments,
          n_nodes = n_nodes
        )

        dual_objective_hermite(
          lambda  = par,
          moments = values,
          powers  = moments,
          quad    = quad_local  ,
          ridge_pen = ridge_pen
        )$value
      },
      
      gr = function(par) {
        lambda_local <- transform_lambda(par, powers)

        quad_local <- adapted_hermite_quad(
          lambda = lambda_local,
          powers = moments,
          n_nodes = n_nodes
        )
        
        dual_objective_hermite(
          lambda  = par,
          moments = values,
          powers  = moments,
          quad    = quad_local ,
          ridge_pen = ridge_pen
        )$gradient
      },
      
      method = "BFGS",
      
      control = list(
        maxit = 1000,
        reltol = optim_tol
      )
    )
    theta_new <- opt$par
    max_step <- 50
    step <- theta_new - theta
    if(max(abs(step)) > max_step) {
      step <- step * max_step / max(abs(step))
    }
    
    theta_new <- theta + step
    #
    # Convergence check
    #
    delta <- max(abs(theta_new - theta))
   theta = theta_new
    
    #
    # Compute achieved moments
    #
    lambda = transform_lambda(theta, powers)
    quad <- adapted_hermite_quad(
      lambda = lambda,
      powers = moments,
      n_nodes = n_nodes
    )
    
    max_lp <- max(quad$log_integrand)
    
    w_unnorm <-
      quad$weights *
      exp(quad$log_integrand - max_lp)
    
    Z <- sum(w_unnorm)
    
    pw <- w_unnorm / Z
    
    moments_achieved <-
      colSums(
        quad$ExpBasis * pw
      )
    moment_errors <-
      abs(moments_achieved - values) / values
    max_moment_error <- max(moment_errors)
    if(verbose) {
      cat(
        sprintf(
          "iter=%d  delta=%.3e  max_err=%.3e\n",
          outer,
          delta,
          max_moment_error
        )
      )
    }
    #
    # Convergence criterion
    #
    if(delta < lambda_tol &&
       max_moment_error < lambda_tol) {
      
      converged <- TRUE
      break
    }
  }
  #
  # Final quadrature
  #
  values = values * (scale_x ^ moments)
  lambda = lambda/(scale_x^moments)
  quad_final <- adapted_hermite_quad(
    lambda = lambda,
    powers = moments,
    n_nodes = 2 * n_nodes
  )
  max_lp = max(quad_final$log_integrand)
  
  w_unnorm =
    quad_final$weights *
    exp(quad_final$log_integrand - max_lp)
  
  Z = sum(w_unnorm)
  
  pw = w_unnorm / Z
  
  log_Z <- max_lp + log(Z)
  
  moments_achieved <-
    colSums(
      quad_final$ExpBasis * pw
    )
  
  moment_errors =
    abs(moments_achieved - values) / values
  
  #
  # Return object
  #
   
  structure(
    list(
      lambda = lambda,
      powers = moments,
      moments_target = values,
      moments_achieved = moments_achieved,
      moment_errors = moment_errors,
      converged = converged,
      #
      # log-space quadrature
      #
      t_nodes = quad_final$t_nodes,
      x_nodes = exp(quad_final$t_nodes),
      weights = quad_final$weights,
      prob_weights = pw,
      ExpBasis = quad_final$ExpBasis,
      log_Z = log_Z,
      mode_t = quad_final$mode,
      sigma_t = quad_final$sigma,
      quad_final=quad_final
    ),
    
    class = "maxent_fit"
  )
}
transform_lambda <- function(theta, powers) {
  
  lambda <- theta
  
  pos_idx <- which.max(powers)
  
  lambda[pos_idx] <- -exp(theta[pos_idx])
  
  neg_idx <- which.min(powers)
  
  if(powers[neg_idx] < 0) {
    lambda[neg_idx] <- -exp(theta[neg_idx])
  }
  
  lambda
}
cdf.maxent_fit <- function(object, grid = NULL, n = 400) {

  # if (is.null(grid)) {
  #   t <- get_hermite_quad(n)$nodes
  #   x <- exp(t)
  # } else {
  #   x <- grid[grid > 0]
  #   t <- log(x)
  # }
  #   
  # cdf.maxent_fit <- function(object, grid = NULL) {
  #   
    x <- object$x_nodes
    p <- object$prob_weights
    
    #
    # keep finite values
    #
    keep <- is.finite(x) &
      is.finite(p)
    
    x <- x[keep]
    p <- p[keep]
    
    #
    # sort
    #
    ord <- order(x)
    
    x <- x[ord]
    p <- p[ord]
    
    #
    # merge duplicate nodes
    #
    ux <- unique(x)
    
    up <- numeric(length(ux))
    
    for(i in seq_along(ux)) {
      up[i] <- sum(p[x == ux[i]])
    }
    
    #
    # normalize
    #
    up <- up / sum(up)
    
    #
    # cumulative distribution
    #
    cdf <- cumsum(up)
    
    #
    # native grid
    #
    if(is.null(grid)) {
      
      return(list(
        x = ux,
        cdf = cdf
      ))
    }
    
    #
    # interpolate
    #
    vals <- approx(
      x = ux,
      y = cdf,
      xout = grid,
      yleft = 0,
      yright = 1,
      ties = "ordered"
    )$y
    
    list(
      x = grid,
      cdf = vals
    )
  }
adaptive_hermite_rule <- function(lambda, powers, n_nodes = 100) {
  
  base <- get_hermite_quad(n_nodes)
  
  z <- base$nodes
  w <- base$weights
  
  # ----- model -----
  eta <- function(t) {
    as.vector(exp(outer(t, powers, "*")) %*% lambda)
  }
  
  d_eta <- function(t) {
    as.vector(exp(outer(t, powers, function(t, k) k * t)) %*% lambda)
  }
  
  d2_eta <- function(t) {
    as.vector(exp(outer(t, powers, function(t, k) k^2 * t)) %*% lambda)
  }
  
  # ----- adaptive proposal (ONLY for sampling) -----
  mode_res <- optimize(eta, c(-20, 20), maximum = TRUE)
  t_mode <- mode_res$maximum
  
  curv <- -d2_eta(t_mode)
  curv <- pmax(curv, 1e-8)
  sigma <- 1 / sqrt(curv)
  
  t_nodes <- t_mode + sigma * z
  
  # ----- FIXED reference measure -----
  log_q <- dnorm(t_nodes, log = TRUE)
  
  # ----- model log-density -----
  log_f <- eta(t_nodes)
  
  # ----- invariant importance weights -----
  log_w <- log_f - log_q
  
  list(
    t_nodes = t_nodes,
    weights = w,
    log_w = log_w,
    ExpBasis = exp(outer(t_nodes, powers, "*")),
    mode = t_mode,
    sigma = sigma
  )
}