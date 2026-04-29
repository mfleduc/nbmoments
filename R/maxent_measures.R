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
# negative integer moments, i.e. moments of the form \eqn{\mathbb{E}[X^k]} for k > 0.
# In that case the exponent includes terms \eqn{lambda_{-k} * x^{-k}}.
#
# References:
#   Mead & Papanicolaou (1984), J. Mathematical Physics 25, 2404-2417
#   Dowson & Wragg (1973), IEEE Trans. Information Theory 19, 689-693
#   Bandyopadhyay et al. (2005), Phys. Rev. E 71, 057701
library(pracma)
library(statmod)
# ==============================================================================
# UTILITIES
# ==============================================================================
#' @title Unnormalized log density
#' @description Evaluate log unnormalized density at quadrature nodes
#' @param lambda numeric vector of Lagrange multipliers
#' @param powers integer vector of moment orders
#' @param nodes numeric vector of quadrature nodes
#' @returns numeric vector of log density values at nodes
#' @export
log_density_nodes <- function(lambda, powers, nodes) {
  log_p <- numeric(length(nodes))
  for (i in seq_along(powers)) {
    log_p <- log_p + lambda[i] * nodes^powers[i]
  }
  log_p
}
#' @title Calculates normalized weights so that \eqn{\int fdx=1}
#' @description Compute normalized weights at quadrature nodes
#' @param lambda numeric vector of Lagrange multipliers
#' @param powers integer vector of moment orders
#' @param nodes numeric vector of quadrature nodes
#' @param weights numeric vector of quadrature weights
#' @returns list with log_Z and normalized probability weights at nodes
#' @export
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
#' @title Gauss-Laguerre based quadrature rule
#' @description Quadrature rule based on Gauss-Laguerre quadrature for estimating the moments
#' @param lambda The Lagrange multipliers for the moments
#' @param powers Integer vector of moment orders
#' @param n_nodes Number of quadrature nodes
#' @return A list describing the quadrature rule
#' @export
gauss_laguerre_quad = function(lambda, moments, n_nodes=100){
  logp = function(x){sum(lambda*x^moments)}
  dlogp = function(x){sum(moments*lambda*x^(moments-1))}
  hesslogp = function(x){sum(moments*(moments-1)*lambda*x^(moments-2))}
  optim_bds = c(-10,10)
  optval = optim(1,logp, gr=dlogp, method="L-BFGS-B",lower=exp(optim_bds[1]),
                 upper=exp(optim_bds[2]),control=list(fnscale=-1 ))
  xhat = optval$par #Mode of the distribution
  lhat = logp(xhat) # Value of log(p) at the mode
  sigma=sqrt(2)#+0.013514+0.003921
  # stopifnot(hesslogp(xhat)<0)
  # sigma = sqrt(-(1/hesslogp(xhat))) # "width" at the mode
  # umin = -xhat/sigma #Lower bound where the function is nonzero
  quad_rule = statmod::gauss.quad( n_nodes, kind="laguerre")
  # quad_rule$weights=quad_rule$weights
  # quad_rule$nodes=quad_rule$nodes*sqrt(2)
  # quad_rule$weights = quad_rule$weights#[quad_rule$nodes>umin]#*sqrt(2)
  # quad_rule$nodes = quad_rule$nodes#[quad_rule$nodes>umin]
  # x_nodes = xhat+sigma*quad_rule$nodes
  x_nodes = quad_rule$nodes
  g_vals = exp(sapply(x_nodes, logp)-lhat+quad_rule$nodes)
  total_weights = sigma*quad_rule$weights*g_vals
  return(
    list(
      nodes=x_nodes,
      weights=total_weights,
      g_nodes=g_vals
    )
    )
}

# ==============================================================================
# DUAL OBJECTIVE
# ==============================================================================

#' @title Dual objective and gradient for the MaxEnt problem
#'
#' @description The primal problem maximizes entropy subject to moment constraints.
#' The dual is: minimize \eqn{D(\lambda) = log Z(\lambda) - \sum_k \lambda_k * c_k}
#' where \eqn{log Z(\lambda) = log \int exp(\sum_k \lambda_k x^k) dx}
#'
#' @param lambda numeric vector of Lagrange multipliers
#' @param moments numeric vector of target moments \eqn{c_k}
#' @param powers integer vector of moment orders k
#' @param nodes quadrature nodes
#' @param weights quadrature weights
#' @returns list with value (scalar dual objective) and gradient
#' @export
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
  if( tail(lambda,1) >0){value=Inf;gradient=Inf}
  list(value = value, gradient = gradient)
}

#' @title Hessian of the dual problem
#' @description Hessian of the dual objective (for diagnostics / Newton steps)
#' \eqn{H_{ij} \approx Cov_p[X^{k_i}, X^{k_j}]} in the Laplace approximation
#'
#' @param lambda numeric vector of Lagrange multipliers
#' @param powers integer vector of moment orders
#' @param nodes quadrature nodes
#' @param weights quadrature weights
#' @returns matrix of second derivatives
#' @export
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
##################
#' @title  maximum entropy measure on \eqn{(0, \infty)}
#'
#' @description Given target moments \eqn{c_k = \int_0^\infty x^k p(x) dx} for k in `powers`,
#' finds the density \eqn{p(x) = \exp(\sum_k \lambda_k x^k)} that maximizes entropy
#' \eqn{S[p] = -\int p(x) log p(x) dx} subject to the moment constraints.
#' The code works by first calculating the maximum entropy distribution given
#' the mean (an exponential) and then one by one incorporating the other moments
#' based on how poorly the current information determines that moment. The moments
#' that match current information least are included first to maximize the additional
#' "information" added to the calculation at each step. See Ximing Wu (2003)
#' Calculation of maximum entropy densities with application to income distribution,
#' Journal of Econometrics, which inspired the idea.
#' Supports positive and negative integer moment orders.
#'
#' @param moments numeric vector of moment values \eqn{c_k}
#' @param powers integer vector of moment orders k (can include negatives)
#'   Example: c(-1L, 1L, 2L) for \eqn{E[X^{-1}], E[X], E[X^2]}
#' @param n_nodes number of quadrature nodes (default 150)
#' @param n_outer maximum L-BFGS iterations (default 2000)
#' @param lambda_tol convergence tolerance for Lagrange multipliers, outer
#' iterations converge when sup-norm of difference beeween consecutive iterations
#' is less than lambda_tol
#' @param optim_tol convergence tolerance on gradient norm (default 1e-7)
#' @param lambda_edges Adjustable parameter: upper bound on Lagrange multipliers
#' for the largest positive and negative moments. Should keep this at 0, things
#' may behave poorly if adjusted. Default 0
#' @param verbose print optimization progress (default FALSE)
#' @returns list with components:
#'   $lambda          Lagrange multipliers (same order as powers)
#'   $powers          moment orders
#'   $nodes           quadrature nodes on (0, inf)
#'   $weights         quadrature weights
#'   $p_nodes         normalized density values at nodes
#'   $converged       logical convergence flag
#'   $moments_target  input moment targets
#'   $moments_achieved moments of the fitted density
#'   $moment_errors   absolute errors in moment constraints
#' @export
maxent_distribution <- function(moments, values, n_nodes = 150,
                           n_outer=2000,lambda_tol=0.5e-5,optim_tol=1e-7,lambda_edges = 0,
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
  # Initially: See what moments we do have
  # If we have -1,1 then start with the GIG guess
  # If we have  just 1 then start with the exponential guess
  # In the future, if we have 1 and 2 then we start with the truncated normal guess
  # Then just nudge stuff to ensure proper decay
  lambda_init = rep(0, length(moments))
  # if(1%in%moments&2%in%moments){
  #   ndx1 = which(moments==1)
  #   ndx2 = which(moments==2)
  #   tnorm_guess = fit_tnorm(values[ndx1],values[ndx2])
  #   lambda_init[ndx1] = tnorm_guess[1]
  #   lambda_init[ndx2] = tnorm_guess[2]
  #   disp(paste("Incorporating moments 1 and 2"))
  #   # lambda_init[ndxn1+1] = log(sqrt(gig_guess$a/gig_guess$b)/(2*besselK( sqrt(gig_guess$a*gig_guess$b),1  )))
  #   # if(any(moments > 1)){  lambda_init[length(moments)]=lambda_edges }
  #   # if(any(moments< -1)){  lambda_init[1]=lambda_edges }
  # }else if(1%in%moments&-1%in%moments){
  #   disp(paste("Incorporating moments 1 and -1"))
  #   ndx1 = which(moments==1)
  #   ndxn1 = which(moments==-1)
  #   ndx0=ndxn1+1
  #   gig_guess = fit_gig( values[ndx1], values[ndxn1] )
  #   lambda_init[ndx1] = gig_guess$lambda_1
  #   lambda_init[ndxn1] = gig_guess$lambda_m1
  #   # if(any(moments > 2)){  lambda_init[length(moments)]=lambda_edges }
  #   # if(any(moments< 0)){  lambda_init[1]=lambda_edges }
  if(1%in%moments){
    ndx1 = which(moments==1)
    lambda_init[ndx1] = -1/values[ndx1]
    disp(paste("Incorporating moment   1"))
    # lambda_init[ndx1-1] = log(values[ndx1])
    # if(any(moments>1)){  lambda_init[length(moments)]=lambda_edges }
    # if(any(moments<0)){  lambda_init[1]=lambda_edges }
  }else{
    max_pos <- which(moments == max(moments[moments > 0]))
    if (length(max_pos) > 0){ lambda_init[max_pos] = lambda_edges }
    lambda_init[1] = lambda_edges
  }
  moments_incorp = which(lambda_init!=0)
  ## Now: Figure out how to incorporate the other moments
  mom_left_out = length(moments)-length(moments_incorp)
  if(mom_left_out==0){
    lambda=lambda_init
    outer=0
  }else{
    #Idea: Whichever ones match the given information least? IDK
    # if(!all(moments==1)){
      quad_init = gauss_laguerre_quad(lambda_init[moments_incorp], moments[moments_incorp], n_nodes=n_nodes  )
      nw <- normalized_weights(lambda_init[moments_incorp], moments[moments_incorp], quad_init$nodes, quad_init$weights)
      log_Z <- nw$log_Z
      pw    <- nw$prob_weights
      est_moments = 0*values
      for(mm in 1:length(moments)){
        est_moments[mm] =  sum(pw * quad_init$nodes^moments[mm])
      }
      est_moments[moments_incorp] = values[moments_incorp]
      next_index = which.max(  abs(est_moments-values)/values )
   #}#else{
    #
    # }
    ##
    while(mom_left_out>0){
      disp(paste("Incorporating moment ",moments[next_index]))
      moments_incorp = unique(sort(c(next_index, moments_incorp)))
      mom_left_out = length(moments)-length(moments_incorp)
      if(mom_left_out==0){n_outer=10*n_outer}
      lambda = lambda_init[moments_incorp]-0.01
      # lambda[lambda==0] = -0.01
      bounds_lower = array(-Inf, dim=c(1, length(lambda)))
      bounds_upper = array(Inf,dim=c(1, length(lambda)) )
      bounds_upper[length(lambda)] = 0#-1e-5
      if(any(moments<0)){bounds_upper[1] =0}# -1e-5}
      for(outer in seq_len(n_outer)){
        # if(lambda[1]==-1e-8){lambda[1] = lambda_edges}
        # if(lambda[length(lambda)]==-1e-8){lambda[length(lambda)] = lambda_edges}
        quadrule = gauss_laguerre_quad(lambda, moments[moments_incorp], n_nodes=n_nodes  )
        if(any(!is.finite(quadrule$weights))){
          # disp(outer)
          quadrule = gauss_laguerre_quad(lambda_prev, moments[moments_incorp], n_nodes=n_nodes  )
        }
        optval = optim(
          par=lambda,
          fn = function(X){
            res = dual_objective(X, values[moments_incorp], moments[moments_incorp],
                                 quadrule$nodes,quadrule$weights )
            return(res$value)
          },
          gr = function(X){
            res =dual_objective(X, values[moments_incorp], moments[moments_incorp],
                                quadrule$nodes,quadrule$weights )
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
        nw = normalized_weights( lambda, moments[moments_incorp], quadrule$nodes, quadrule$weights   )
        moments_achieved = sapply(moments[moments_incorp], function(X){sum(nw$prob_weights*quadrule$nodes^X)} )
        #
        # disp(abs(moments_achieved-values[moments_incorp]) /abs(values[moments_incorp]))
        moment_delta = max(abs(moments_achieved-values[moments_incorp])/abs(values[moments_incorp]))
        delta = max(abs(lambda-lambda_new))
        # disp(delta)
        if((delta<lambda_tol)|(moment_delta<lambda_tol)){break}
        lambda_prev=lambda
        lambda=lambda_new
      }
      lambda_init[moments_incorp]=lambda
      nw <- normalized_weights(lambda_init[moments_incorp], moments[moments_incorp], quadrule$nodes, quadrule$weights)
      log_Z <- nw$log_Z
      pw    <- nw$prob_weights
      est_moments = 0*values
      for(mm in 1:length(moments)){
        est_moments[mm] =  sum(pw * quadrule$nodes^moments[mm])
      }
      est_moments[moments_incorp] = values[moments_incorp]
      next_index = which.max(  abs(est_moments-values)/values )
    }
    lambda=lambda_new
  }

  quad_final = gauss_laguerre_quad(lambda, moments, n_nodes*2)
  if(any(is.nan(quad_final$weights))){
    quad_final = gauss_laguerre_quad(lambda, moments, n_nodes)
  }
  nw = normalized_weights( lambda, moments, quad_final$nodes, quad_final$weights*exp(quad_final$lhat)  )

  moments_achieved = numeric(length(moments))
  names(moments_achieved) = as.character(moments)
  moments_achieved = sapply( moments, function(X){sum(nw$prob_weights*quad_final$nodes^X)} )
  moment_errors=abs(moments_achieved-values)/values
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
#' @export
print.maxent_fit <- function(x, ...) {
  cat("Maximum Entropy density on (0, inf)\n")
  cat("Moment orders:", x$powers, "\n")
  cat("Converged:", x$converged, "\n")
  cat("Max moment error:", max(x$moment_errors), "\n")
  invisible(x)
}

#' Summary method for maxent_fit
#'@export
summary.maxent_fit <- function(object, probs = c(0.80, 0.95),grid=NA, ...) {
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
  ci <- credible_interval(object, probs,grid=grid)
  cat("\nCredible intervals:\n")
  for (p in 1:length(probs)) {

    cat(sprintf("  %3.0f%%: (%.6f, %.6f)  width=%.6f\n",
                100 * probs[p], ci$lower[p], ci$upper[p], ci$upper[p]-ci$lower[p]))
  }
  invisible(object)
}
#' @title Predict method for maxent_fit
#' @description Evaluate MaxEnt density at new points
#'
#' @param object maxent_fit object
#' @param newdata numeric vector of x values
#' @return numeric vector of normalized density values
#' @export
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

#' @title Calculate equal-tail credible interval
#' @description Compute equal-tail credible interval from MaxEnt density
#'
#' @param object maxent_fit object
#' @param prob probability mass (default 0.95)
#' @return numeric vectors with lower and upper bounds
#' @export
credible_interval <- function(object, prob = 0.95,grid=NA) {
  if(any(is.na(grid))){
    quadrule = gauss_laguerre_quad( object$lambda,object$powers, n_nodes = 10*length(object$nodes) )
    ord    <- order(quadrule$nodes)
    x_ord  <- quadrule$nodes[ord]
    pw_ord = predict(object, x_ord)

    cdf    <- cumsum(pw_ord*quadrule$weights[ord])
    cdf    <- cdf / max(cdf)
  }else{
    pred_pdf = predict(object, grid)
    cdf = cumtrapz( grid, pred_pdf )
    cdf    <- cdf / max(cdf)
    x_ord=grid
  }
  lower=upper=array(0,dim=c(1, length(prob)))
  for(ii in 1:length(prob)){
  alpha  <- (1 - prob[ii]) / 2
  lower[ii]  <- x_ord[which(cdf >= alpha)[1]]
  upper[ii]  <- x_ord[which(cdf >= 1 - alpha)[1]]
  }
  return(list(prob=as.array(prob),lower = as.array(lower), upper = as.array(upper)))
}

#' @title Plot the MaxEnt density
#' @description S3 method for plotting the maxent density
#' @param x maxent_fit object
#' @param n_plot number of points for smooth curve
#' @param log_scale plot on log x scale (default TRUE)
#' @param add add to existing plot (default FALSE)
#' @param ... additional arguments passed to plot/lines
plot.maxent_fit <- function(x, n_plot = 500, log_scale = FALSE,
                            add = FALSE, ...) {
  # Use a smooth grid for plotting
  node_range <- range(x$nodes[x$prob_weights > max(x$prob_weights) * 1e-50])
  # node_range <- range(x$nodes)
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

##################
#' @title Fit the truncated normal distribution
#' @description Solve for the maximum entropy measure on (0, inf) given a mean and second moment
#'
#' Given a target mean and second moment, calculates the truncated normal distribution that maximizes the entropy given those constraints.
#'
#' @param mn Desired mean of the distribution
#' @param sqmn Desired value of \eqn{E[X^2]}
#' @returns vector of lagrange multipliers such that the posterior distribution is a truncated normal.
#' @export
fit_tnorm = function(mn, sqmn){
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
  opt = optim( c(mn,sqrt(sqmn-mn^2)) , function(X)tnorm_obj(X,c(mn,sqmn)),
               method = "CG",control = list(maxit = 50000, reltol = 1e-8))
  #This gives us the parameters, now we need to write it in terms of lagrange multipliers
  lambda1 = opt$par[1]/opt$par[2]^2
  lambda2 = -1/2/opt$par[2]^2
  return(c(lambda1,lambda2))
}
# ==============================================================================
# GIG SPECIAL CASE
# ==============================================================================
#' @title Fit a generalized inverse Gaussian with p=1
#' @description Fit GIG distribution by matching \eqn{E[X]} and \eqn{E[X^{-1}]}
#' The GIG is the MaxEnt distribution for these two moment constraints when p=1.
#' @param mean_val  \eqn{E[X]}
#' @param inv_mean  \eqn{E[X^{-1}]}
#' @return list with GIG parameters p, a, b and achieved moments
#'   Parameterization: \eqn{p(x) \propto x^{p-1} \exp(-(a*x + b/x)/2)}
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
#' @title GIG PDF
#' @description Evaluate GIG density
#'
#' @param x numeric vector of evaluation points
#' @param p GIG parameter p
#' @param a GIG parameter a (> 0)
#' @param b GIG parameter b (> 0)
#' @return numeric vector of density values
#' @export
dgig <- function(x, p, a, b) {
  sqrt_ab <- sqrt(a * b)
  log_const <- p / 2 * log(a / b) - log(2 * besselK(sqrt_ab, p))
  exp(log_const + (p - 1) * log(x) - 0.5 * (a * x + b / x))
}

