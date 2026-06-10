# Calculate the MLEs and the MLE of theta at mu = mu0
MLE.mu = function(x, mu0=1){
  n = length(x); xb = mean(x); # MLE of mu
  ss = function(x, a){
    s = 0; ss = 0; k= max(x)
    for(i in 1:(k-1)){
      ni = sum(x >= (i+1)); ai = i/(1+a*i)
      s = s + ni*ai
      ss = ss + ni*ai^2
    }
    return(c(s,ss))}
  fa = function(a){
    sc = ss(x,a); sums = sc[1]
    res = log(1+a*xb)-a*xb+a^2*sums/n
    return(res)}
  k = 1; a0 =(var(x)/xb-1)/xb
  j = 0; lr = .0001; ur = a0
  while(fa(lr)*fa(ur) > 0 & j <=30){
    ur = ur + .5
    j = j + 1}
  if(j > 30){
    an = a0}
  else{
    an = uniroot(fa, c(lr, ur))[[1]]} # an = MLE of alpha
  th = 1/an # th = MLE of theta
  fc = function(a){
    sc = ss(x,a); sums = sc[1]
    res = log(1+a*mu0)-a*mu0*(1+a*xb)/(1+a*mu0)+a^2*sums/n
    return(res)}
  k = 1; a0 =an
  lt = .0001; ut = an; j = 1
  while(fc(lt)*fc(ut) >0 & j <= 30){
    ut = ut + .5
    j = j + 1}
  if(j > 30){ml0 = a0}
  else{
    ml0 = uniroot(fc, c(lt, ut))[[1]]}
  return(c(xb, 1/an, 1/ml0))
}
# Computes the SLRT stat at mu = mu0
LRTstat = function(x, muh, th, mu0,th0){
  n = length(x);
  xb = muh
  LF = function(x,mu, th){
    
    # JOURNAL OF STATISTICAL COMPUTATION AND SIMULATION 21
    
    n = length(x); xb = mean(x)
    al = 1/th; lv = numeric(0)
    slv = 0; k= max(x)
    for(i in 1:(k-1)){
      ni = sum(x >= (i+1)); ai =log(1+al*i)
      slv = slv + ni*ai}
    y = xb*log(mu)-(xb+th)*log(1+al*mu)+slv/n
    return(y)
  }
  st = -2*n*(LF(x,mu0,th0)-LF(x,muh,th))
  if(st < 0){st = 0}
  stat = sign(muh-mu0)*sqrt(st)
  return(c(stat))
}
# One-Sample CI for the mean: Score and Likelihood CIs
ci.mu = function(x, cl){
  n = length(x)
  al = (1-cl)/2; zu = qnorm(1-al); zl = -zu
  mls = MLE.mu(x); xb = muh = mls[1]; th = mls[2]
  calsq = zu^2;
  cent = (xb+calsq/2/n)/(1-calsq/n/th)
  me = (zu/sqrt(n))*sqrt(calsq/4/n+xb+xb^2/th)/(1-calsq/n/th)
  Lscr = cent- me; Uscr = cent + me
  fn = function(mu){
    ml0 = MLE.mu(x, mu); muh = ml0[1]; th = ml0[2]; th0 = ml0[3]
    return(LRTstat(x, muh, th, mu,th0)-zu)} # lower limit
  fu = function(mu){
    ml0 = MLE.mu(x, mu); muh = ml0[1]; th = ml0[2]; th0 = ml0[3]
    return(LRTstat(x, muh, th, mu,th0)-zl)} #upper limit
  lc = max(.001, Lscr-.1); uc = muh
  while(fn(lc)*fn(uc) > 0){
    lc = max(0.0001, lc - .1)}
  Low = uniroot(fn, c(lc, uc))[[1]]
  lc = muh; uc = Uscr
  while(fu(lc)*fu(uc) > 0){
    uc = uc + .1}
  Upp = uniroot(fu, c(muh, Uscr+1))[[1]]
  return(c(Lscr, Uscr, Low, Upp))
}
# Confidence Interval for the Ratio of Means
MOV.Ratio.CI = function(x, y, cl){
  al = (1-cl)/2; zcrt = qnorm(1-al)
  cix = ci.mu(x, cl); lx = cix[1]; ux = cix[2]; mux = (lx+ux)/2
  ciy = ci.mu(y, cl); ly = ciy[1]; uy = ciy[2]; muy = (ly+uy)/2
  pr = mux*muy
  ltr2 = muy^2-(muy-uy)^2; ltr1 = mux^2-(mux-lx)^2
  Lscr = (pr - sqrt(pr^2-ltr1*ltr2))/ltr2
  ltr2 = muy^2-(muy-ly)^2; ltr1 = mux^2-(mux-ux)^2
  Uscr = ( pr + sqrt(pr^2-ltr1*ltr2))/ltr2
  #
  lx = cix[3]; ux = cix[4]; mux = (lx+ux)/2
  ly = ciy[3]; uy = ciy[4]; muy = (ly+uy)/2
  
  # 22 K. KRISHNAMOORTHY AND M. M. HASAN
  pr = mux*muy
  ltr2 = muy^2-(muy-uy)^2; ltr1 = mux^2-(mux-lx)^2
  Llrt = (pr - sqrt(pr^2-ltr1*ltr2))/ltr2
  ltr2 = muy^2-(muy-ly)^2; ltr1 = mux^2-(mux-ux)^2
  Ulrt = ( pr + sqrt(pr^2-ltr1*ltr2))/ltr2
  return(c(Lscr, Uscr, Llrt, Ulrt))
}
# Confidence Interval for the Difference between Means
MOV.Diff.CI = function(x, y, cl){
  al = (1-cl)/2; zcrt = qnorm(1-al)
  cix = ci.mu(x, cl);
  ciy = ci.mu(y, cl);
  lx = cix[1]; ux = cix[2]; mux = (lx+ux)/2
  ly = ciy[1]; uy = ciy[2]; muy = (ly+uy)/2
  Lscr = mux-muy - sqrt((mux-lx)^2+(muy-uy)^2)
  Uscr = mux-muy + sqrt((mux-ux)^2+(muy-ly)^2)
  lx = cix[3]; ux = cix[4]; mux = (lx+ux)/2
  ly = ciy[3]; uy = ciy[4]; muy =(ly+uy)/2
  Llrt = mux-muy - sqrt((mux-lx)^2+(muy-uy)^2)
  Ulrt = mux-muy + sqrt((mux-ux)^2+(muy-ly)^2)
  return(c(Lscr, Uscr, Llrt, Ulrt))
}
#Example
# x1 = c(3, 3, 0, 1, 9, 8, 0, 23, 6, 6, 6,8, 6, 12, 10, 0, 3, 28,
#        
#        2, 6, 3,+ 3, 3, 2, 76, 2, 4, 13)
# 
# x2 = c(9, 9, 3, 1, 7, 1, 19, 7, 0, 7, 2,1, 4, 0, 25, 3, 6, 2, 8,
#        
#        0, 72,+ 2, 5, 1, 28, 4, 4, 19, 0, 0, 3)