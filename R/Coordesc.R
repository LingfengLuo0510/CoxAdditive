#' This function implements the coordinate-wise gradient boosting algorithm
#' 
#' @param delta event indicator
#' @param z Covariate matrix
#' @param time the death time
#' @param p number of variables
#' @param N the sample size
#' @param knot the number of basis functions for non-linear predictors
#' @param maxit the number of maximum iterations
#' @param rate the learning rate
#' @param tol the convergence threshold
#' 
#' @return return a list of the following values
#' beta :  1 denotes selected variables and 0 denotes unselected variables.
#' time :  the used time
#' convergence : whether it's converge within maximum iteration steps
#' key  :  the number of iterations



CoorDesc <- function(delta, z, time, p ,N, knot , maxit = 300, rate = 0.05, tol = 1e-3){
  
  # set initial values
  time_old      <- proc.time()           # record the time used
  key           <- 0                     # the number of iterations 
  likelihood_all<- NULL                  # the likelihood
  j_star_all    <- NULL                  # the chosen optimal block
  index         <- rep(1:p, each=knot-1) # the index to record the chosen beta 
  converge      <- FALSE                 # record whether it converges
  z_spline      <- NULL                  
  bs8           <- list(rep(0,p))
  beta          <- rep(0, p*(knot-1))    # record beta
  
  #Using B-splines to span the g function
  for (j in 1:p){
    knot_set    <- quantile(z[,j],prob=seq(1:(knot-4))/(knot-3))
    bs7         <- bs(z[,j],df=knot, knot=knot_set, intercept=FALSE, degree=3)
    bs8[[j]]    <- matrix(bs7, nrow=N) # for build beta_t
    z_spline    <- cbind(z_spline, bs8[[j]])
  }
  
  #start the loop
  repeat{
    key       <- key + 1
    offset    <- z_spline%*%beta
    GD        <- NULL
    likelihood<- NULL
    update_all<- NULL
    #z_spline=NULL
    #L1_summary=NULL
    #L2_summary=NULL
    for (j in 1:p){
      beta_j     <- rep(0,knot-1)                                # initilize         
      result     <- ddloglik(N,delta,bs8[[j]],beta_j,offset)     # call the rcpp function to get first-order derivative, and partial likelihood                                  
      update     <- as.vector(result$L1)                         # record the first derivate
      update_all <- rbind(update_all,update)                   
      GD[j]      <- sum(result$L1*as.numeric(update))            
      likelihood[j] <- sum(result$partial_likelihood)/N
      # L1_summary=cbind(L1_summary,result$L1)
      #L2_summary=cbind(L2_summary,result$L2)
      #z_spline=cbind(z_spline, bs8)
    }
    
    j_star              <- which(GD==max(GD))                    # choose the maximum as the optimal block.
    j_star_all          <- c(j_star_all, j_star)
    beta[index==j_star] <- beta[index==j_star] + rate*update_all[j_star,]
    likelihood_all      <- c(likelihood_all, likelihood[j_star])
    
    # check the convergence
    # use relative likelihood ratio change here. 
    if (key>=10) {
      llk.diff = abs((likelihood_all[key]-likelihood_all[key-1])/(likelihood_all[key-1]))
      if(llk.diff < tol) {
        converge <- TRUE    # relative likelihood ratio change is less than the threshold, ,mark the result's convergence true.
        break
      }
    }
    # if the loop gets the maximum iteration, we break the loop.
    if (key==maxit){
      break
    }
  }#end loop
  
  # get the chosen beta 
  beta_key    <- (abs(beta)>0)
  group_key   <- index[beta_key]
  select_group_key  <- unique(group_key)
  select_group_key2 <- rep(0,p)
  select_group_key2[select_group_key] <- 1
  time_used   <- proc.time()-time_old
  
  return( list(beta = select_group_key2, time = time_used[3] ,convergence = converge,key = key ) )  
}



