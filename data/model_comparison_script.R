## TESTING: Compare frequentist behavior to that of Krishnamoorthy and Hassan
set.seed(16846)
p_test=c(0.05,0.25, 0.5, 0.75, 0.95)
library("pracma")
library(PoissonRatioUQ)
n <- 10
n2<- 10
r1 <- 10
m1 <- 0.25
p1 <- 1/(1+m1)
r2 <- 2
m2 <- 5
p2 <- 1/(1+m2)
ellmax <- 120 #
true_ratio = (r1*m1)/(r2*m2)
###
counter=0
moments = list(
 c(1,2), c(-1,1),c(-2,-1,1),c(-1,1,2),c(-2,-1,1,2)
)
ntrials= 100
nfail=nfail2=0
in_int = array(0, dim=c(length(p_test)+1,ntrials,length(moments)+2))
for(rth in 1:ntrials){
  data1 <- rnbinom(n , r1, p1)
  data2 <- rnbinom(n2 , r2, p2)
  # data1=c( 7 ,15 ,13,  6,  8,  7, 14,  4, 15, 21,  5,  8, 11,  7, 13, 10,  5,
  #          8, 18,  2, 11,  7,  7,  5, 13)
  # data2=c( 6,  0, 12, 16,  7, 28,  0, 12, 22, 16,  4, 37,  7,  7, 17,  4, 12,  3,  8,  3,  8, 18,  8, 22,  0)
  stepflag=T
  p=1
  # disp(credible_ints)
  # kh_conf_ints = MOV.Ratio.CI(data1, data2, p_test[4])[1:2]
  # disp(kh_conf_ints)
  while(stepflag==T){
     
    kh_conf_ints = try(MOV.Ratio.CI(data1, data2, p_test[p])[1:2])
    if(is.numeric(kh_conf_ints)){
      if((true_ratio>kh_conf_ints[1])&(true_ratio<kh_conf_ints[2])){
        in_int[p:length(p_test),rth,1]=1
        stepflag=F
      }else if(p==length(p_test)){
        in_int[p+1,rth,1]=1
        stepflag=F
      }else{
        p=p+1
      }
    }else{
      nfail2=nfail2+1
      stepflag=F
    }
      
  }
  stepflag=T
  p=1
  # disp(credible_ints)
  while(stepflag==T){
    kh_conf_ints = try(MOV.Ratio.CI(data1, data2, p_test[p])[3:4])
    # disp(kh_conf_ints)
    if(is.numeric(kh_conf_ints)){
      if((true_ratio>kh_conf_ints[1])&(true_ratio<kh_conf_ints[2])){
        in_int[p:length(p_test),rth,2]=1
        stepflag=F
      }else if(p==length(p_test)){
        in_int[p+1,rth,2]=1
        stepflag=F
      }else{
        p=p+1
      }
    }else{
      nfail=nfail+1
      stepflag=F
    }
    }
  
  for(mm in 1:length(moments)){
    # disp(mm)
    medist = estimate_dist(data1, data2, moments = moments[[mm]],  ridge_pen=0.0032, 
                           ellmax=150, n_nodes = 150,n_outer=150)
    # disp(medist$moment_errors)
    credible_ints = credible_interval.maxent_fit(medist, prob = p_test,
                                                 grid=seq(1e-4, 10, by=1e-4))
    p=1
    stepflag=T
    if(mm==1){
      disp(c(credible_ints$lower[4] ,credible_ints$upper[4]))
    }
    while(stepflag==T){
      if((true_ratio>credible_ints$lower[p])&(true_ratio<credible_ints$upper[p])){
        in_int[p:length(p_test),rth,mm+2]=1
        stepflag=F
      }else if(p==length(p_test)){
        in_int[p+1,rth,mm+2]=1
        stepflag=F
      }else{
        p=p+1
      }
    }
  }
}

AA=apply(in_int, MARGIN = c(1,3),sum)