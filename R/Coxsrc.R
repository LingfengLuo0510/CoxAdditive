AR1 <- function(tau, m) {
  if(m==1) {R <- 1}
  if(m > 1) {
    R <- diag(1, m)
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        R[i,j] <- R[j,i] <- tau^(abs(i-j))
      }
    }
  }
  return(R)
}
#' function for standardize covariate
normalize = function(x){
  y = sqrt(length(x))*(x-mean(x))/sqrt(sum((x-mean(x))^2))
  return(y)
}

#' The evaluation function 
#' @return a list of the following values
#' FP:  false positive 
#' FN:  false negative
#' Se:  Senstitivity
#' Sp:  Specificity
#' FDP: False positive portion
#' 
#' @export
FPFNSeSpLik=function(TrueBeta=TrueBeta,beta=beta){
  FP <- length(which(TrueBeta==0 & beta!=0))
  FN <- length(which(TrueBeta!=0 & beta==0))
  Se <- length(which(TrueBeta!=0 & beta!=0))/length(which(TrueBeta!=0))
  Sp <- length(which(TrueBeta==0 & beta==0))/length(which(TrueBeta==0))
  FDP=FP/max(length(which(beta!=0)),1)
  output=c(FP, FN, Se, Sp, FDP)
  return(output)
}


#' Simulation function
#' the function for simulating data
#' @param N sample size
#' @param p number of parameters
#' @param p_true number of true parameters
#' @param m1
#' @param knot
#' @import mvtnorm
#' @return a list of the following values
#' 
#' delta: event indicator
#' z: Covariate matrix
#' time: the death time
#' p: number of variables
#' N: the sample size
#' knot: the number of basis functions for non-linear predictors
#' 
#' 
#' @export
simul  <- function(N = 1000, p = 500, p_true = 5, m1 = 100, knot = 10){
  Sigma_z1    <- diag(p)
  Corr1       <- AR1(0.3,m1)
  diag(Corr1) <- 1
  z           <- NULL
  j           <- 0
  
  while(j<(p/m1)){
    j    <- j+1
    z    <- cbind(z,rmvnorm(N, mean=rep(0,m1), sigma=Corr1))
  }
  
  z      <- apply(z,2,normalize)
  z2     <- z**2
  
  ####True beta
  TrueBeta       <- rep(0, p)
  TrueBeta_index <- sample(1:p,p_true,replace=FALSE)
  signbeta       <- sample(c(-1,1),p_true,replace=T)
  #signbeta=rep(1,p_true)
  mag            <- runif(p_true, 1,2)
  TrueBeta[TrueBeta_index]  <- mag*signbeta
  xbeta          <- z2%*%TrueBeta
  U              <- runif(N, 0, 1)
  pre_time       <- -log(U)/(1*exp(xbeta))
  pre_censoring  <- runif(N,1,30)
  pre_censoring  <- pre_censoring*(pre_censoring<3)+3*(pre_censoring>=3)
  tcens          <- (pre_censoring<pre_time) # censoring indicator
  delta          <- 1-tcens
  time           <- pre_time*(delta==1)+pre_censoring*(delta==0)
  delta          <- delta[order(time)]
  z              <- z[order(time),]
  time           <- time[order(time)]
  
  return( list(delta=delta,  z =z, time = time, p=p, N = N, knot = knot, TrueBeta = TrueBeta) )
  
}






