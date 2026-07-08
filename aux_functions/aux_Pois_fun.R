## Online Supplement to manuscript
## Submitted to The American Statistician
##
## 07 / 2026

## Aux functions for Poisson regression, from scracth 




## -- Poisson Linking function                         ####
## P01_Pois_Ch01
linkingFun <- function(beta_iter, X){             
  exp(X %*% beta_iter)
}

## -- Poisson Likelihood functions                     ####
## P01_Pois_Ch02

PoisLik <- function(lambda_iter, y){
  n_prod <- length(y)
  prod_arg <- sapply(1:n_prod, function(i){
    (lambda_iter[i]^y[i]/factorial(y[i]))*exp(-lambda_iter[i]) })
  prod(prod_arg) }


## -- Poisson LogLikelihood                            ####
loglik_pois <- function(y, lambda_iter){         
  Loglik_fun<-(t(y)%*%log(lambda_iter))-sum(lambda_iter)-sum(lgamma(y+1)) 
}


## -- Poisson gradient for Newton-Raphson              ####
gradient_pois <- function(lambda_iter, X, y){
  residue_iter <- y - lambda_iter  
  stacked_g <- sapply(1:nrow(X), function(i){                                       
    residue_iter[i]*X[i,]})
  if(!is.matrix(stacked_g)){ g_gradient <- sum(stacked_g)
  } else { g_gradient <- as.matrix(apply(stacked_g,1,sum))}
  return(list(gradient = g_gradient, stacked_g = stacked_g,
              pois_residue =  residue_iter)) 
}






## -- Poisson Hessian for Newton-Raphson               ####
hessian_pois <- function(lambda_iter, X){
  n_regressors <- ncol(X)
  n_observ <- nrow(X)
  Outer_Prod_X_sum<-matrix(0L,nrow=n_regressors, ncol=n_regressors)
  for(i in 1:n_observ){
    x <- X[i,]
    x <- as.numeric(x)
    outprod_X <- (x %o% x)
    Outer_Prod_X_sum<-Outer_Prod_X_sum-lambda_iter[i]*outprod_X
  }
  Outer_Prod_X_sum
}


## -- "half stepping" for iterative optimisation       ####
myHalfStepping <- function(Xnr, beta_iter, Hessian_iter, g_gradient){
  accept_candidate <- FALSE
  max_iter_halfstep <- 1000
  lambda_accepted <- linkingFun(beta_iter, Xnr)
  iter_halfstep <- 0
  stepsize <- 1
  target_iter <- qr.solve(Hessian_iter, g_gradient, tol = 1e-18)            
  while(!accept_candidate & (iter_halfstep < max_iter_halfstep)){
    temp_delta <- beta_iter - stepsize*target_iter
    lambda_iter <- linkingFun(temp_delta, Xnr)  
    iter_halfstep  <- iter_halfstep + 1
    criterion_half_stepping <- as.numeric(loglik_pois(y, lambda_iter)) >  
      as.numeric(loglik_pois(y, lambda_accepted))
    if(!criterion_half_stepping){ stepsize <- stepsize/2 } else { 
      accept_candidate <- TRUE 
      target_iter <- stepsize*target_iter                                        ## Fixed May 2026, previously omitted
    }
  } 
  list(target_iter = target_iter, temp_delta = temp_delta, stepsize = stepsize)
}      


## -- Poisson Newton-Raphson iterative optimisation    ####
NR_MLE <- function(Xnr, ynr){
  if(!is.matrix(Xnr)){Xnr <- as.matrix(Xnr)}
  if(!is.matrix(ynr)){ynr <- as.matrix(ynr)}
  n_regressors <- ncol(Xnr)
  n_observ <- nrow(Xnr)
  beta_iter <- matrix(0, nrow=n_regressors, ncol = 1)
  rownames(beta_iter) <- colnames(Xnr)
  colnames(beta_iter) <- "estimate"
  ln_target <- as.matrix(log(abs(ynr)))
  ln_target[which(!is.finite(ln_target))]<-0                       
  beta_iter[1] <- mean(ln_target)                                
  iter_count <- 0
  iter_max <- 200               
  tol <- rep(1e-6, n_regressors)
  iter_stop <- FALSE
  hessian_is_singular <- FALSE
  while(iter_count <= iter_max & !iter_stop){
    iter_count <- iter_count + 1
    lambda_iter <- linkingFun(beta_iter, Xnr)          
    grad <- gradient_pois(lambda_iter, Xnr, ynr)           
    g_gradient <- grad$gradient
    Hessian_iter <- hessian_pois(lambda_iter, Xnr)  
    inv_test <- try(qr.solve(Hessian_iter, g_gradient), silent = TRUE)
    if(inherits(inv_test, "try-error")){
      hessian_is_singular <- TRUE
      iter_stop <- TRUE } 
    if(!hessian_is_singular){
      HalfStep_proc <- myHalfStepping(Xnr, beta_iter, Hessian_iter, g_gradient)
      target_iter <- HalfStep_proc$target_iter
      temp_delta <-  HalfStep_proc$temp_delta                                  
    } else { warning("Hessian is singular")
      target_iter <- 0 
      stepsize  <- 1
    }
    if(length(which(sqrt(target_iter^2) < tol)) < n_regressors){
      beta_iter <- temp_delta } else {iter_stop <- TRUE} }
  return(list(coef_est = beta_iter, predict_y = lambda_iter, 
              Hessian = Hessian_iter, iter = iter_count  ))
}


