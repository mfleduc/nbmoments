#################################################################
# Functions for moment calculation via polynomial approximation #
#################################################################
#
#'@title Moment computation for the negative binomial posteror
#'@description Uses the polynomial approximation of Bradlow et al as described in LeDuc and Kissler 2026
#' to estimate moments of the negative binomial posterior distribution. Can compute both negative and positive integer moments
#' with conditions on the moments described in the paper.
#' @param k Integer: The moment order, \eqn{\mathbb{E}[X^k]}
#' @param data Array: The count data with which to estimate the moment
#' @param ell_max Integer, the maximum degree of the polynomial approximation to the Gamma function
#' @param a,b the parameters for the Pearson-VI prior for the parameter m (odds ratio). Default
#' a=0, b=max(2, k+1). The prior has the form \eqn{m^a/(1+m)^b}.
#' @param c,d,z1,z2 The parameters for the PEarson-VI prior for the parameter r
#' (dispersion parameter). Defaults c=0,d=1,z1=0,z2=-1. The prior takes the form
#' \eqn{(r-z1)^c/(r-z2)^d}. Note that the distribution is defined only for 0>=z1>z2>-sum(data), and
#' is supported on \eqn{(-z1, \infty)}. A minimally informative prior must set z1=0, and both z1 and z2 should be integers.
#' @returns The approximated moment, using the method in LeDuc and Kissler 2026.
#' @export
compute_moment <- function( k, data,ell_max=50,a=1, b=max(2, k+1), c=0, d=1 , z1=0,z2=-1 ){
  #First: A lof of stuff needs to be precalculated
  stopifnot( abs(k-round(k))<10^-12   )
  k=round(k)
  min_0k = min(0,k)
  max_0k = max(0,k)
  #
  n = length(data)
  xstar = max(data)
  x_sum=sum(data)
  if(b<=k){b=k+1}
  K1 <- (a+b+x_sum+1)/n #K parameter
  #
  nk1 <- array(0, dim=c(1, xstar+2))
  for(kk in 0:(xstar+1)){
    nk1[kk+1] <- sum(data==kk)
  }
  # sk1 <- 0*nk1
  # for(kk in 0:(1+xstar)){
  #   sk1[kk+1] <- sum(nk1[(kk+1):(xstar+2)])
  # }
  sk1 <- rev(cumsum(rev(nk1)))
  vt1 <- wt1 <- 0*seq(1, 1+xstar)
  vt1[1] <- 1
  wt1[abs(z2)+1] <- -1
  # cs1 <- a+x_sum+1+ c(min_0k,max_0k)#C_k values
  Umin = compute_U_polynomial( min_0k, x_sum, a, ell_max )
  Umax = compute_U_polynomial( max_0k, x_sum, a, ell_max )
  ## Compute a_j "stably" (or at least moreso)
  hi1 <- ( sk1[2:(xstar+2)]+c*vt1+d*wt1 );
  h1 <- sum(hi1)
  acoeffs1 <- 0*seq(0,hi1[1])
  acoeffs1[1] <- 1
  for(ii  in 1:(xstar)){
    tmpcoeffs <- (ii)^seq(hi1[ii+1],0,by=-1)*choose(hi1[[ii+1]],seq(hi1[ii+1],0,by=-1))
    acoeffs1 <- polymul(tmpcoeffs/max(tmpcoeffs), acoeffs1)
    # acoeffs1 = acoeffs1/max(acoeffs1)
  }
  acoeffs1 <- acoeffs1[length(acoeffs1):1]
  # Q terms now. The real test
  nz <- max(4,sum( abs(acoeffs1)==0 ))
  # Q_min <- Q_max <- 0*seq(0,h1)
  Q_min <- numeric(h1 + 1)
  Q_max <- numeric(h1 + 1)
  for(jj in (nz-3):(h1)){
    betafnterm <- lbeta(a+x_sum+1-jj+(0:ell_max)-2, jj+1+k)
    K1term <- log(1/K1)*(a+x_sum+1-jj+(0:ell_max)-2)
    Uterm <- log(Umin$coeffs)-log(n)*(0:ell_max)
    Q_min[jj+1] <- sum(exp(betafnterm+(K1term)+(Uterm)))
  }
  for(jj in (nz-3):(h1)){
    betafnterm <- lbeta(a+x_sum+1-jj+(0:ell_max)-2, jj+1+k)
    K1term <- log(1/K1)*(a+x_sum+1-jj+(0:ell_max)-2)
    Uterm <- log(Umax$coeffs)-log(n)*(0:ell_max)
    Q_max[jj+1] <- sum(exp(betafnterm+(K1term)+(Uterm)))
  }
  if(k>0){
    # moment = sum(acoeffs1*Q_max)/sum(acoeffs1*Q_min)*(1/n)^k*exp(sum(log(a+x_sum+(1:k)-1)))
    moment = sum(acoeffs1*Q_max)/sum(acoeffs1*Q_min)*(1/n)^k*exp(lgamma(a + x_sum + k) - lgamma(a + x_sum))
  }else if(k<0){
    moment = sum(acoeffs1*Q_min)/sum(acoeffs1*Q_max)*(1/n)^k/exp(sum(log(a+x_sum+(k:-1))))
  }else if(k==0){
    moment=1
  }
  return(moment)
}
#'@title Computation of the polynomial approximant to the Gamma function.
#'@description Computes the polynomial approximation to the Gamma function used in the
#'moment estimation. Does so in a relatively stable manner by normalizing the approximation
#'at each step.
#'@param k The moment to be approximated, \eqn{\mathbb{E}[X^k]}
#'@param x_sum the sum of the data
#'@param a The a parameter for the prior distribution of m
#'@param ell_max the maximum degree of the polynomial approximation
#'@returns The coefficients of the polynomial U as well as the log of the scaling
#'@export
compute_U_polynomial <- function(k, x_sum, a, ell_max) {
  # C_k <- a + x_sum + k + 1
  # poly <- c(1, rep(0, ell_max))
  # total_log_scale <- 0
  # 
  # for (j in seq_len(C_k)) {
  #   factor <- j^(0:ell_max)
  #   new_poly <- numeric(ell_max + 1)
  #   for (p in 0:ell_max) {
  #     for (q in 0:(ell_max - p)) {
  #       new_poly[p + q + 1] <- new_poly[p + q + 1] +
  #                              exp(sum(log(poly[p + 1])+log(factor[q + 1])))
  #     }
  #   }
  #   # poly = new_poly
  #   # Normalize at each step
  #   max_val <- max(abs(new_poly[new_poly != 0]))
  #   if (max_val > 0 && is.finite(max_val)) {
  #     total_log_scale <- total_log_scale + log(max_val)
  #     poly <- new_poly / max_val
  #   } else {
  #     poly <- new_poly
  #   }
  # }
  # list(coeffs = poly[(1+ell_max):1], logscale = total_log_scale)
  # compute_U_polynomial <- function(k, x_sum, a, ell_max) {
  C_k <- a + x_sum + k + 1
  poly <- numeric(ell_max + 1)
  poly[1] <- 1
  total_log_scale <- 0
  factor <- numeric(ell_max + 1)
  new_poly <- numeric(ell_max + 1)
  for (j in seq_len(C_k)) {
    
    
    factor[1] <- 1
    
    for(q in 1:ell_max){
      factor[q + 1] <- factor[q] * j
    }
    
    # new_poly <- numeric(ell_max + 1)
    
    for (p in 0:ell_max) {
      
      poly_p <- poly[p + 1]
      
      if(poly_p == 0)
        next
      
      max_q <- ell_max - p
      
      for (q in 0:max_q) {
        
        new_poly[p + q + 1] <-
          new_poly[p + q + 1] +
          poly_p * factor[q + 1]
      }
    }
    
    max_val <- max(abs(new_poly))
    
    if (max_val > 0 && is.finite(max_val)) {
      
      total_log_scale <- total_log_scale + log(max_val)
      
      poly <- new_poly / max_val
      
      poly[abs(poly) < 1e-15] <- 0
      
    } else {
      poly <- new_poly
    }
  }
  list(
    coeffs = rev(poly),
    logscale = total_log_scale
  )
}
#'@title Moment computation for the negative binomial posterior via the Tricomi expansion
#'@description Uses the Tricomi expansion
#' to estimate moments of the negative binomial posterior distribution. Still a
#' work in progress, but hopefully will allow generalization to arbitrary real moments
#' @param k Integer: The moment order,\eqn{\mathbb{E}[X^k]}
#' @param data Array: The count data with which to estimate the moment
#' @param a,b the parameters for the Pearson-VI prior for the parameter m (odds ratio). Default
#' a=0, b=max(2, k+1). The prior has the form \eqn{m^a/(1+m)^b}.
#' @param c,d,z1,z2 The parameters for the PEarson-VI prior for the parameter r
#' (dispersion parameter). Defaults c=0,d=1,z1=0,z2=-1. The prior takes the form
#' \eqn{(r-z1)^c/(r-z2)^d}. Note that the distribution is defined only for 0>=z1>z2>-sum(data), and
#' is supported on \eqn{(-z1, \infty)}. A minimally informative prior must set z1=0, and both z1 and z2 should be integers.
#' @returns Right now, an error. Eventually, the approximated moment using a method based on the Tricomi expansion.
#' @export
compute_moment_tricomi <- function( k, data,a=1, b=max(2, k+1), c=0, d=1 , z1=0,z2=-1 ){
  #First: A lof of stuff needs to be precalculated
  stop("This function is under development and should not be used yet")
  stopifnot( abs(k-round(k))<10^-12   )
  k=round(k)
  min_0k = min(0,k)
  max_0k = max(0,k)
  #
  n = length(data)
  xstar = max(data)
  x_sum=sum(data)
  K1 <- (a+b+x_sum+1)/n #K parameter

}
#'@title Convolution of two vectors via FFT
#'@description Uses the FFT to convolve two vectors
#' @param a,b the vectors to convolve
#' @returns The convolution
#' @export
fft_convolve <- function(a, b) {
  n <- length(a) + length(b) - 1
  N <- 2^ceiling(log2(n))
  
  fa <- fft(c(a, rep(0, N - length(a))))
  fb <- fft(c(b, rep(0, N - length(b))))
  
  Re(fft(fa * fb, inverse = TRUE))[1:n] / N
}





