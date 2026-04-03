# Maximum Entropy Measure on (0, infinity)
# Given finitely many raw integer moments c_k = integral x^k p(x) dx
# Finds the density p(x) maximizing entropy S[p] = -integral p(x) log p(x) dx
# subject to moment constraints.
#
# The maximum entropy density has the exponential family form:
#   p(x) = exp(lambda_0 + lambda_1 * x + lambda_2 * x^2 + ... + lambda_m * x^m)
# where the lambda_i are Lagrange multipliers determined by the moment constraints.
#
# For the strong Stieltjes moment problem on (0, infinity) we also support
# negative integer moments, i.e. moments of the form E[X^{-k}] for k > 0.
# In that case the exponent includes terms lambda_{-k} * x^{-k}.
#
# References:
#   Mead & Papanicolaou (1984), J. Mathematical Physics 25, 2404-2417
#   Dowson & Wragg (1973), IEEE Trans. Information Theory 19, 689-693
#   Bandyopadhyay et al. (2005), Phys. Rev. E 71, 057701

# ==============================================================================
# UTILITIES
# ==============================================================================
library("statmod")
#' Log-sum-exp (no external dependencies)
logSumExp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#' Evaluate log unnormalized density at quadrature nodes
#'
#' @param lambda numeric vector of Lagrange multipliers
#' @param powers integer vector of moment orders
#' @param nodes numeric vector of quadrature nodes
#' @return numeric vector of log density values at nodes
log_density_nodes <- function(lambda, powers, nodes) {
  log_p <- numeric(length(nodes))
  for (i in seq_along(powers)) {
    log_p <- log_p + lambda[i] * nodes^powers[i]
  }
  log_p
}

#' Compute normalized weights at quadrature nodes
#'
#' @param lambda numeric vector of Lagrange multipliers
#' @param powers integer vector of moment orders
#' @param nodes numeric vector of quadrature nodes
#' @param weights numeric vector of quadrature weights
#' @return list with log_Z and normalized probability weights at nodes
normalized_weights <- function(lambda, powers, nodes, weights) {
  log_p <- log_density_nodes(lambda, powers, nodes)
  max_lp <- max(log_p)
  unnorm <- weights * exp(log_p - max_lp)
  Z <- sum(unnorm)
  list(
    log_Z = max_lp + log(Z),
    prob_weights = unnorm / Z
  )
}

#' #' Generate quadrature nodes and weights on (0, inf) via log transformation
#' #' x = exp(t), t uniformly spaced, trapezoidal rule
#' #' Jacobian: integral f(x) dx = integral f(e^t) e^t dt
#' #'
#' #' @param n number of nodes
#' #' @param t_min lower bound in log space
#' #' @param t_max upper bound in log space
#' #' @return list with nodes (x values) and weights
log_transform_quadrature <- function(n, t_min = -10, t_max = 10) {
  t_nodes <- seq(t_min, t_max, length.out = n)
  h <- t_nodes[2] - t_nodes[1]
  x_nodes <- exp(t_nodes)
  # Trapezoidal weights with Jacobian x = e^t
  weights <- x_nodes * h
  weights[1] <- weights[1] / 2
  weights[n] <- weights[n] / 2
  list(nodes = x_nodes, weights = weights)
}
#'
#' Gauss-Hermite quadrature varying with lambda
#' @param lambda The Lagrange multipliers for the moments
#' @param powers Integer vector of moment orders
#' @param n_nodes Number of quadrature nodes
#' @return The approximate values of the moments?
#'
laplace_gauss_hermite_quad = function(lambda, moments, n_nodes=100){
  logp = function(x){sum(lambda*x^moments)}
  dlogp = function(x){sum(moments*lambda*x^(moments-1))}
  hesslogp = function(x){sum(moments*(moments-1)*lambda*x^(moments-2))}
  optim_bds = c(-10,10)
  optval = optim(1,logp, gr=dlogp, method="L-BFGS-B",lower=exp(optim_bds[1]),
                 upper=exp(optim_bds[2]),control=list(fnscale=-1 ))
  xhat = optval$par #Mode of the distribution
  lhat = logp(xhat) # Value of log(p) at the mode
  stopifnot(hesslogp(xhat)<0)
  sigma = sqrt(-(1/hesslogp(xhat))) # "width" at the mode
  umin = -xhat/sigma #Lower bound where the function is nonzero
  quad_rule = statmod::gauss.quad( n_nodes, kind="hermite")
  # quad_rule$weights=quad_rule$weights
  # quad_rule$nodes=quad_rule$nodes*sqrt(2)
  quad_rule$weights = quad_rule$weights[quad_rule$nodes>umin]*sqrt(2)
  quad_rule$nodes = quad_rule$nodes[quad_rule$nodes>umin]
  x_nodes = xhat+sigma*quad_rule$nodes
  g_vals = exp(sapply(x_nodes, logp)-lhat+quad_rule$nodes^2)
  total_weights = sigma*quad_rule$weights*g_vals
  return(
    list(
      nodes=x_nodes,
      weights=total_weights,
      sigma=sigma,
      lhat=lhat,
      g_nodes=g_vals
    )
    )
  # block_u = function(u)(max(umin,u))
  #This function evaluates g(u)
  # eval_g = function(u){
  #   if(u<umin){
  #     return( 0 )
  #   }else{
  #     return( exp(logp(xhat+sigma*u) +(u^2)-lhat) )
  #   }
  # }
  # g_vals = sapply(quad_rule$nodes, eval_g)
  # #Now: Do the moment calculations
  # power_portion = outer( xhat+sigma*quad_rule$nodes, powers, "^" )
  # these_moments = colSums( g_vals*power_portion*quad_rule$weights  )*sigma*exp(lhat)
  # return(these_moments)
}
#' Set up quadrature adapted to the moment information
#' Centers the log-space grid near the expected location of the mass
#'
#' @param moments numeric vector of moment values
#' @param powers integer vector of moment orders
#' @param n_nodes number of quadrature nodes
#' @return list with nodes and weights
setup_quadrature <- function(moments, powers, n_nodes = 4000) {
  pos <- powers >= 0
  if (sum(pos) >= 2 && 1 %in% powers && 2 %in% powers) {
    mu  <- moments[powers == 1]
    m2  <- moments[powers == 2]
    cv2 <- max(m2 / mu^2 - 1, 0.01)
    t_center <- log(mu) - 0.5 * log(1 + cv2)
    t_spread  <- sqrt(log(1 + cv2)) * 4 + 3
  } else if (1 %in% powers) {
    mu <- moments[powers == 1]
    t_center <- log(max(mu, 1e-6))
    t_spread  <- 20
  } else {
    t_center <- 0
    t_spread  <- 20
  }
  log_transform_quadrature(n_nodes,#t_min=-8,t_max=8)
                           t_min = t_center - t_spread,
                           t_max = t_center + t_spread)
}

# ==============================================================================
# DUAL OBJECTIVE
# ==============================================================================

#' Dual objective and gradient for the MaxEnt problem
#'
#' The primal problem maximizes entropy subject to moment constraints.
#' The dual is: minimize D(lambda) = log Z(lambda) - sum_k lambda_k * c_k
#' where log Z(lambda) = log integral exp(sum_k lambda_k x^k) dx
#'
#' @param lambda numeric vector of Lagrange multipliers
#' @param moments numeric vector of target moments c_k
#' @param powers integer vector of moment orders k
#' @param nodes quadrature nodes
#' @param weights quadrature weights
#' @return list with value (scalar dual objective) and gradient
dual_objective <- function(lambda, moments, powers, nodes, weights) {
  nw <- normalized_weights(lambda, powers, nodes, weights)
  log_Z <- nw$log_Z
  pw    <- nw$prob_weights

  # Gradient: d log Z / d lambda_k = E_p[X^k]
  grad_log_Z <- numeric(length(powers))
  for (i in seq_along(powers)) {
    grad_log_Z[i] <- sum(pw * nodes^powers[i])
  }

  value    <- log_Z - sum(lambda * moments)
  gradient <- grad_log_Z - moments

  list(value = value, gradient = gradient)
}

#' Hessian of the dual objective (for diagnostics / Newton steps)
#'
#' H_{ij} = Cov_p[X^{k_i}, X^{k_j}]
#'
#' @param lambda numeric vector of Lagrange multipliers
#' @param powers integer vector of moment orders
#' @param nodes quadrature nodes
#' @param weights quadrature weights
#' @return matrix of second derivatives
dual_hessian <- function(lambda, powers, nodes, weights) {
  nw <- normalized_weights(lambda, powers, nodes, weights)
  pw <- nw$prob_weights
  m  <- length(powers)
  H  <- matrix(0, m, m)
  for (i in seq_len(m)) {
    xi <- nodes^powers[i]
    ei <- sum(pw * xi)
    for (j in i:m) {
      xj <- nodes^powers[j]
      ej <- sum(pw * xj)
      H[i, j] <- sum(pw * xi * xj) - ei * ej
      H[j, i] <- H[i, j]
    }
  }
  H
}

# ==============================================================================
# MAIN SOLVER
# ==============================================================================

#' Solve for the maximum entropy measure on (0, inf)
#'
#' Given target moments c_k = integral_0^inf x^k p(x) dx for k in `powers`,
#' finds the density p(x) = exp(sum_k lambda_k x^k) that maximizes entropy
#' S[p] = -integral p(x) log p(x) dx subject to the moment constraints.
#'
#' Supports positive and negative integer moment orders.
#'
#' @param moments numeric vector of moment values c_k
#' @param powers integer vector of moment orders k (can include negatives)
#'   Example: c(-1L, 1L, 2L) for E[X^{-1}], E[X], E[X^2]
#' @param n_nodes number of quadrature nodes (default 1000)
#' @param tol convergence tolerance on gradient norm (default 1e-8)
#' @param max_iter maximum L-BFGS iterations (default 2000)
#' @param verbose print optimization progress (default FALSE)
#' @return list with components:
#'   $lambda          Lagrange multipliers (same order as powers)
#'   $powers          moment orders
#'   $nodes           quadrature nodes on (0, inf)
#'   $weights         quadrature weights
#'   $p_nodes         normalized density values at nodes
#'   $converged       logical convergence flag
#'   $moments_target  input moment targets
#'   $moments_achieved moments of the fitted density
#'   $moment_errors   absolute errors in moment constraints
maxent_measure <- function(moments, powers, n_nodes = 4000,
                           tol = 1e-8, max_iter = 10000,
                           verbose = FALSE) {

  stopifnot(length(moments) == length(powers))
  stopifnot(all(is.finite(moments)))
  stopifnot(all(moments > 0 | powers < 0))  # positive moments for k >= 0

  powers <- as.integer(powers)
  quad   <- setup_quadrature(moments, powers, n_nodes)
  nodes  <- quad$nodes
  weights <- quad$weights

  # Initial lambda: zeros except nudge lambda for highest positive power
  # to encourage decay at infinity
  lambda_init <- rep(0, length(powers))
  max_pos <- which(powers == max(powers[powers > 0]))
  if (length(max_pos) > 0) lambda_init[max_pos] <- -0.1

  # Optimization via L-BFGS
  opt <- optim(
    par     = lambda_init,
    fn      = function(lv) {
      res <- dual_objective(lv, moments, powers, nodes, weights)
      if (!is.finite(res$value)) return(1e15)
      res$value
    },
    gr      = function(lv) {
      res <- dual_objective(lv, moments, powers, nodes, weights)
      if (!all(is.finite(res$gradient))) return(rep(0, length(powers)))
      res$gradient
    },
    method  = "L-BFGS-B",
    control = list(
      maxit  = max_iter,
      factr  = (tol / .Machine$double.eps) * 1e7,
      pgtol  = tol,
      trace  = if (verbose) 1L else 0L
    )
  )
  lambda_opt <- opt$par
  converged  <- opt$convergence == 0

  if (!converged && verbose) {
    message("L-BFGS-B: ", opt$message)
  }

  # Compute normalized density at nodes
  nw      <- normalized_weights(lambda_opt, powers, nodes, weights)
  p_nodes <- nw$prob_weights / weights  # density = prob_weight / quad_weight

  # Verify moment constraints
  moments_achieved <- numeric(length(powers))
  for (i in seq_along(powers)) {
    moments_achieved[i] <- sum(weights * nodes^powers[i] *
                                 nw$prob_weights / weights * weights)
    # simplifies to:
    moments_achieved[i] <- sum(nw$prob_weights * nodes^powers[i])
  }
  names(moments_achieved) <- as.character(powers)
  moment_errors <- abs(moments_achieved - moments)
  names(moment_errors) <- as.character(powers)

  if (verbose) {
    cat("\nMoment verification:\n")
    for (i in seq_along(powers)) {
      cat(sprintf("  k=%3d: target=%.8f  achieved=%.8f  error=%.2e\n",
                  powers[i], moments[i],
                  moments_achieved[i], moment_errors[i]))
    }
  }

  structure(
    list(
      lambda           = lambda_opt,
      powers           = powers,
      nodes            = nodes,
      weights          = weights,
      p_nodes          = nw$prob_weights / weights,
      prob_weights     = nw$prob_weights,
      log_Z            = nw$log_Z,
      converged        = converged,
      optim_result     = opt,
      moments_target   = moments,
      moments_achieved = moments_achieved,
      moment_errors    = moment_errors
    ),
    class = "maxent_fit"
  )
}
##################
#' Solve for the maximum entropy measure on (0, inf)
#'
#' Given target moments c_k = integral_0^inf x^k p(x) dx for k in `powers`,
#' finds the density p(x) = exp(sum_k lambda_k x^k) that maximizes entropy
#' S[p] = -integral p(x) log p(x) dx subject to the moment constraints.
#'
#' Supports positive and negative integer moment orders.
#'
#' @param moments numeric vector of moment values c_k
#' @param powers integer vector of moment orders k (can include negatives)
#'   Example: c(-1L, 1L, 2L) for E[X^{-1}], E[X], E[X^2]
#' @param n_nodes number of quadrature nodes (default 1000)
#' @param tol convergence tolerance on gradient norm (default 1e-8)
#' @param max_iter maximum L-BFGS iterations (default 2000)
#' @param verbose print optimization progress (default FALSE)
#' @return list with components:
#'   $lambda          Lagrange multipliers (same order as powers)
#'   $powers          moment orders
#'   $nodes           quadrature nodes on (0, inf)
#'   $weights         quadrature weights
#'   $p_nodes         normalized density values at nodes
#'   $converged       logical convergence flag
#'   $moments_target  input moment targets
#'   $moments_achieved moments of the fitted density
#'   $moment_errors   absolute errors in moment constraints
maxent_distribution <- function(moments, values, n_nodes = 150,
                           n_outer=2500,lambda_tol=0.5e-5,optim_tol=1e-7,lambda_edges = -15,
                           verbose = FALSE) {
  stopifnot(length(moments) == length(values))
  stopifnot(any(moments>0))
  stopifnot(all(is.finite(values)))
  stopifnot(all(values > 0 ))  # positive moments since integrand is non-negative
  # if(!any(moments==0)){
  #   moments <- append(moments, 0)
  #   values <- append(values, 1)
  # }
  moments_sort <- sort.int(moments, index.return=TRUE)
  moments <- moments_sort$x
  values <- values[moments_sort$ix]
  #
  # special cases to start things off
  #
  # if(all(moments==c(-1,0,1))){
  #   result = fit_gig( values[3], values[1] )
  #   l0 = log(sqrt(result$a/result$b)/(2*besselK( sqrt(result$a*result$b),1  )))
  #   # Ensures the distribution normalizes correctly
  #   return(structure(
  #     list(
  #       lambda=c( lambda_m1, l0,lambda_1 ),
  #       moments_target=values,
  #       moments_achieved=c(result$moments_achieved[2],1,result$moments_achieved[1]),
  #       converged=result$converged,
  #       moment_error = abs( ( values-moments_achieved ) )/values,
  #       powers=moments,
  #     ),class="maxent_fit"
  #   )
  #   )
  # }else if(all(moments==c(0,1))){
  #   l0 = log( values[2] )
  #   l1 = -values[2]
  #   return(structure(
  #     list(
  #       lambda=c( l0,l1 ),
  #       moments_target=values,
  #       moments_achieved=values,
  #       converged=TRUE,
  #       moment_error = c(0,0),
  #       powers=moments,
  #     ),class="maxent_fit"
  #   )
  #   )
  # }
  # Initially: See what moments we do have
  # If we have -1,1 then start with the GIG guess
  # If we have  just 1 then start with the exponential guess
  # In the future, if we have 1 and 2 then we start with the truncated normal guess
  # Then just nudge stuff to ensure proper decay
  lambda_init <- rep(0, length(moments))
  if(1%in%moments&-1%in%moments){
    ndx1 = which(moments==1)
    ndxn1 = which(moments==-1)
    ndx0=ndxn1+1
    gig_guess = fit_gig( values[ndx1], values[ndxn1] )
    lambda_init[ndx1] = gig_guess$lambda_1
    lambda_init[ndxn1] = gig_guess$lambda_m1
    # lambda_init[ndxn1+1] = log(sqrt(gig_guess$a/gig_guess$b)/(2*besselK( sqrt(gig_guess$a*gig_guess$b),1  )))
    if(any(moments > 1)){  lambda_init[length(moments)]=lambda_edges }
    if(any(moments< -1)){  lambda_init[1]=lambda_edges }
  }else if(1%in%moments){
    ndx1 = which(moments==1)
    lambda_init[ndx1] = -1/values[ndx1]
    # lambda_init[ndx1-1] = log(values[ndx1])
    if(any(moments>1)){  lambda_init[length(moments)]=lambda_edges }
    if(any(moments<0)){  lambda_init[1]=lambda_edges }
  }else{
    max_pos <- which(moments == max(moments[moments > 0]))
    if (length(max_pos) > 0){ lambda_init[max_pos] = lambda_edges }
    lambda_init[1] = lambda_edges
  }
  lambda = lambda_init
  bounds_lower = array(-Inf, dim=c(1, length(lambda)))
  bounds_upper = array(Inf,dim=c(1, length(lambda)) )
  bounds_upper[length(lambda)] = -1e-8
  if(any(moments<0)){bounds_upper[1] = -1e-8}
  for(outer in seq_len(n_outer)){
    # if(lambda[1]==-1e-8){lambda[1] = lambda_edges}
    # if(lambda[length(lambda)]==-1e-8){lambda[length(lambda)] = lambda_edges}
    quadrule = laplace_gauss_hermite_quad(lambda, moments, n_nodes=n_nodes  )
    if(any(!is.finite(quadrule$weights))){
      # disp(outer)
      quadrule = laplace_gauss_hermite_quad(lambda_prev, moments, n_nodes=n_nodes  )
    }
    optval = optim(
      par=lambda,
      fn = function(X){
        res = dual_objective(X, values, moments, quadrule$nodes,quadrule$weights )
        return(res$value)
      },
      gr = function(X){
        res = dual_objective(X, values, moments, quadrule$nodes,quadrule$weights )
        return(res$gradient)
      },
      lower=bounds_lower,
      upper=bounds_upper,
      method="L-BFGS-B",
      control=list(
        maxit=1000,
        factr= 1e6*(optim_tol/.Machine$double.eps),
        pgtol=optim_tol
      )
    )

    # moment_est =
    lambda_new=optval$par
    # disp(lambda_new)
    # disp(outer)
    #What are the moments?
    nw = normalized_weights( lambda, moments, quadrule$nodes, quadrule$weights   )
    moments_achieved = sapply( moments, function(X){sum(nw$prob_weights*quadrule$nodes^X)} )
    #
    disp(abs(moments_achieved-values) /abs(values))
    moment_delta = max(abs(moments_achieved-values)/abs(values))
    delta = max(abs(lambda-lambda_new))
    # disp(delta)
    if((delta<lambda_tol)|(moment_delta<lambda_tol)){break}
    lambda_prev=lambda
    lambda=lambda_new

  }
  lambda=lambda_new
  quad_final = laplace_gauss_hermite_quad(lambda, moments, n_nodes*2)
  if(any(is.nan(quad_final$weights))){
    quad_final = laplace_gauss_hermite_quad(lambda, moments, n_nodes)
  }
  nw = normalized_weights( lambda, moments, quad_final$nodes, quad_final$weights*exp(quad_final$lhat)  )

  moments_achieved = numeric(length(moments))
  names(moments_achieved) = as.character(moments)
  moments_achieved = sapply( moments, function(X){sum(nw$prob_weights*quad_final$nodes^X)} )
  moment_errors=abs(moments_achieved-values)
  if (verbose) {
    cat("\nMoment verification:\n")
    for (i in seq_along(moments)) {
      cat(sprintf("  k=%3d: target=%.8f  achieved=%.8f  error=%.2e\n",
                  moments[i], values[i],
                  moments_achieved[i], moment_errors[i]))
    }
  }

  structure(
    list(
      lambda           = lambda,
      powers           = moments,
      nodes            = quad_final$nodes,
      weights          = quad_final$weights,
      p_nodes          = nw$prob_weights / quad_final$weights,
      prob_weights     = nw$prob_weights,
      log_Z            = nw$log_Z,
      converged        = (outer<n_outer),
      optim_result     = optval,
      moments_target   = values,
      moments_achieved = moments_achieved,
      moment_errors    = abs(moments_achieved-values)/abs(values)
    ),
    class = "maxent_fit"
  )
}
# ==============================================================================
# S3 METHODS
# ==============================================================================

#' Print method for maxent_fit
print.maxent_fit <- function(x, ...) {
  cat("Maximum Entropy density on (0, inf)\n")
  cat("Moment orders:", x$powers, "\n")
  cat("Converged:", x$converged, "\n")
  cat("Max moment error:", max(x$moment_errors), "\n")
  invisible(x)
}

#' Summary method for maxent_fit
summary.maxent_fit <- function(object, probs = c(0.80, 0.95), ...) {
  cat("=== Maximum Entropy Density on (0, inf) ===\n\n")
  cat("Moment orders:", object$powers, "\n")
  cat("Converged:    ", object$converged, "\n\n")

  cat("Lagrange multipliers:\n")
  for (i in seq_along(object$powers)) {
    cat(sprintf("  lambda[k=%3d] = %12.6f\n",
                object$powers[i], object$lambda[i]))
  }

  cat("\nMoment constraints:\n")
  for (i in seq_along(object$powers)) {
    cat(sprintf("  E[X^%3d]: target=%12.6f  achieved=%12.6f  rel. error=%.2e\n",
                object$powers[i], object$moments_target[i],
                object$moments_achieved[i], object$moment_errors[i]))
  }

  # Summary statistics from density
  pw <- object$prob_weights
  nd <- object$nodes

  mean_val <- sum(pw * nd)
  m2_val   <- sum(pw * nd^2)
  var_val  <- m2_val - mean_val^2
  sd_val   <- sqrt(max(var_val, 0))

  cat(sprintf("\nMean:     %.6f\n", mean_val))
  cat(sprintf("Std Dev:  %.6f\n", sd_val))
  cat(sprintf("CV:       %.4f\n", sd_val / mean_val))

  cat("\nCredible intervals:\n")
  for (p in probs) {
    ci <- credible_interval(object, p)
    cat(sprintf("  %3.0f%%: (%.6f, %.6f)  width=%.6f\n",
                100 * p, ci[1], ci[2], diff(ci)))
  }
  invisible(object)
}

#' Evaluate MaxEnt density at new points
#'
#' @param object maxent_fit object
#' @param newdata numeric vector of x values
#' @return numeric vector of normalized density values
predict.maxent_fit <- function(object, newdata, ...) {
  log_p <- numeric(length(newdata))
  for (i in seq_along(object$powers)) {
    log_p <- log_p + object$lambda[i] * newdata^object$powers[i]
  }
  # Normalize using partition function estimated from quadrature
  log_p_nodes <- log_density_nodes(object$lambda, object$powers, object$nodes)
  max_lp <- max(log_p_nodes)
  Z <- sum(object$weights * exp(log_p_nodes - max_lp))
  exp(log_p - max_lp) / Z
}

#' Compute credible interval from MaxEnt density
#'
#' @param object maxent_fit object
#' @param prob probability mass (default 0.95)
#' @return named numeric vector with lower and upper bounds
credible_interval <- function(object, prob = 0.95) {
  ord    <- order(object$nodes)
  x_ord  <- object$nodes[ord]
  pw_ord <- object$prob_weights[ord]

  cdf    <- cumsum(pw_ord)
  cdf    <- cdf / max(cdf)

  alpha  <- (1 - prob) / 2
  lower  <- x_ord[which(cdf >= alpha)[1]]
  upper  <- x_ord[which(cdf >= 1 - alpha)[1]]
  c(lower = lower, upper = upper)
}

#' Plot the MaxEnt density
#'
#' @param x maxent_fit object
#' @param n_plot number of points for smooth curve
#' @param log_scale plot on log x scale (default TRUE)
#' @param add add to existing plot (default FALSE)
#' @param ... additional arguments passed to plot/lines
plot.maxent_fit <- function(x, n_plot = 500, log_scale = TRUE,
                            add = FALSE, ...) {
  # Use a smooth grid for plotting
  node_range <- range(x$nodes[x$prob_weights > max(x$prob_weights) * 1e-6])
  x_plot <- exp(seq(log(node_range[1]), log(node_range[2]), length.out = n_plot))
  p_plot <- predict(x, x_plot)

  xlab <- if (log_scale) "x (log scale)" else "x"
  if (add) {
    lines(x_plot, p_plot, ...)
  } else {
    plot(x_plot, p_plot, type = "l",
         xlab = xlab, ylab = "Density",
         main = "Maximum Entropy Density",
         log  = if (log_scale) "x" else "",
         ...)
  }
  invisible(x)
}
#===============================================================================
#Truncated normal special case
#===============================================================================
tnorm_obj <- function(params, moments){
  #Calculate the parameters of the truncated normal distribution given the first two moments
  # Rewrite to get this like fig_gig
  snr <- params[1]/params[2]
  dnsnr <- dnorm(-snr)
  Z = 1-pnorm(-snr)
  eX <- params[1]+dnsnr*params[2]/Z
  varX <- params[2]^2*(1-snr*dnsnr/Z-(dnsnr/Z)^2)
  err <- (eX-moments[1])^2/moments[1]^2 + ( varX+eX^2-moments[2] )^2/moments[2]^2
}
# ==============================================================================
# GIG SPECIAL CASE
# ==============================================================================

#' Fit GIG distribution by matching E[X] and E[X^{-1}]
#' The GIG is the MaxEnt distribution for these two moment constraints.
#'
#' @param mean_val  E[X]
#' @param inv_mean  E[X^{-1}]
#' @return list with GIG parameters p, a, b and achieved moments
#'   Parameterization: p(x) propto x^{p-1} exp(-0.5*(a*x + b/x))
fit_gig <- function(mean_val, inv_mean) {
  # MaxEnt GIG has fixed p = 1
  # p(x) propto exp(lambda_1 * x + lambda_{-1} / x)
  # = exp(-a/2 * x - b/2 / x)  with a = -2*lambda_1, b = -2*lambda_{-1}
  # Moments for GIG(p=1, a, b):
  # E[X]    = sqrt(b/a) * K_2(sqrt(ab)) / K_1(sqrt(ab))
  # E[1/X]  = sqrt(a/b) * K_0(sqrt(ab)) / K_1(sqrt(ab))

  gig1_moments <- function(log_a, log_b) {
    a <- exp(log_a)
    b <- exp(log_b)
    sqrt_ab <- sqrt(a * b)
    K0 <- besselK(sqrt_ab, 0)
    K1 <- besselK(sqrt_ab, 1)
    K2 <- besselK(sqrt_ab, 2)
    c(
      sqrt(b / a) * K2 / K1,   # E[X]
      sqrt(a / b) * K0 / K1    # E[1/X]
    )
  }

  obj <- function(par) {
    m <- tryCatch(gig1_moments(par[1], par[2]),
                  error = function(e) c(NA, NA))
    if (any(!is.finite(m))) return(1e10)
    (m[1] - mean_val)^2 / mean_val^2 +
      (m[2] - inv_mean)^2  / inv_mean^2
  }

  # Initial guess: a ~ 1/mean, b ~ 1/inv_mean
  init <- c(log(1 / mean_val), log(1 / inv_mean))
  opt  <- optim(init, obj, method = "Nelder-Mead",
                control = list(maxit = 10000, reltol = 1e-14))

  a_fit <- exp(opt$par[1])
  b_fit <- exp(opt$par[2])
  sqrt_ab <- sqrt(a_fit * b_fit)
  K0 <- besselK(sqrt_ab, 0)
  K1 <- besselK(sqrt_ab, 1)
  K2 <- besselK(sqrt_ab, 2)

  list(
    p = 1,
    a = a_fit,
    b = b_fit,
    lambda_1  = -a_fit / 2,
    lambda_m1 = -b_fit / 2,
    converged = opt$convergence == 0,
    moments_achieved = c(
      "E[X]"    = sqrt(b_fit / a_fit) * K2 / K1,
      "E[X^-1]" = sqrt(a_fit / b_fit) * K0 / K1
    )
  )
}
#' Evaluate GIG density
#'
#' @param x numeric vector of evaluation points
#' @param p GIG parameter p
#' @param a GIG parameter a (> 0)
#' @param b GIG parameter b (> 0)
#' @return numeric vector of density values
dgig <- function(x, p, a, b) {
  sqrt_ab <- sqrt(a * b)
  log_const <- p / 2 * log(a / b) - log(2 * besselK(sqrt_ab, p))
  exp(log_const + (p - 1) * log(x) - 0.5 * (a * x + b / x))
}

