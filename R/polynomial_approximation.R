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
  K1 <- (a+b+x_sum+1)/n #K parameter
  #
  nk1 <- array(0, dim=c(1, xstar+2))
  for(kk in 0:(xstar+1)){
    nk1[kk+1] <- sum(data==kk)
  }
  sk1 <- 0*nk1
  for(kk in 0:(1+xstar)){
    sk1[kk+1] <- sum(nk1[(kk+1):(xstar+2)])
  }
  vt1 <- wt1 <- 0*seq(1, 1+xstar)
  vt1[1] <- 1
  wt1[abs(z2)+1] <- -1
  cs1 <- a+x_sum+1+ c(min_0k,max_0k)#C_k values
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
  Q_min <- Q_max <- 0*seq(0,h1)
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
    moment = sum(acoeffs1*Q_max)/sum(acoeffs1*Q_min)*(1/n)^k*exp(sum(log(a+x_sum+(1:k)-1)))
  }else if(k<0){
    moment = sum(acoeffs1*Q_min)/sum(acoeffs1*Q_max)*(1/n)^k/exp(sum(log(a+x_sum+(k:-1))))
  }else if(k==0){
    moment=1
  }
  return(moment)
}
compute_U_polynomial <- function(k, x_sum, a, ell_max) {
  C_k <- a + x_sum + k + 1
  poly <- c(1, rep(0, ell_max))
  total_log_scale <- 0
  
  for (j in seq_len(C_k)) {
    factor <- j^(0:ell_max)
    new_poly <- numeric(ell_max + 1)
    for (p in 0:ell_max) {
      for (q in 0:(ell_max - p)) {
        new_poly[p + q + 1] <- new_poly[p + q + 1] + 
                               exp(sum(log(poly[p + 1])+log(factor[q + 1])))
      }
    }
    # Normalize at each step
    max_val <- max(abs(new_poly[new_poly != 0]))
    if (max_val > 0 && is.finite(max_val)) {
      total_log_scale <- total_log_scale + log(max_val)
      poly <- new_poly / max_val
    } else {
      poly <- new_poly
    }
  }
  list(coeffs = poly[(1+ell_max):1], logscale = total_log_scale)
}








