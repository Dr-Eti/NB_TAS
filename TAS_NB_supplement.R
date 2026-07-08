## Online Supplement to manuscript
## Submitted to The American Statistician
##
## 07 / 2026

## Input prep
## -- load functions for auxiliary Poisson regression        ####
source("aux_functions/aux_Pois_fun.R")   


## -- read data                                              ####
mypath <- "./data"
myData <- read.csv(paste0(mypath, "/TAS_NB_example.csv", collapse = " "), header = TRUE, row.names=1)


target_feat <- "g_DISRUPTIONS"
regressors_lablels <- colnames(myData)[-which(colnames(myData) %in% target_feat)]

y <- myData[,target_feat, drop = FALSE]                     # target
y <- as.matrix(y)
n <- nrow(myData)                                           # observations
p <- length(regressors_lablels)                             # features excluding intercept

X <-  cbind(1, myData)                                      # add intercept 
X <- as.matrix(X[,-which(colnames(X) %in% target_feat)])
colnames(X) <- c("intercept", regressors_lablels)
rownames(X) <- rownames(myData)
n_regressors <- ncol(X)                                     # p + 1











##                                                           ####
## Neg Bin cross-section                                     ####
## -- own NB PMF                                             ####
## The pre-built function ?pnbinom uses a different parametrisation of the PMF. This one uses lambda and alpha

dnegbin_own <-   function(x, p_lambda, p_alpha){
  if(p_alpha == 0){
    theta_nb <- gamma_alpha <- 0
  } else {
    theta_nb <- 1/p_alpha
    gamma_alpha <-gamma(p_alpha)
  }
  sapply(1:length(p_lambda), function(j){
    mu <- p_lambda[j]                                                         # equivalent to p_alpha*a/b
    a <- mu/(mu + p_alpha)
    b <- 1-a                                                                  # equivalent to p_alpha/(p_alpha + mu)
    #bin_coeff <- (gamma(x + p_alpha)/(gamma(x + 1)*gamma_alpha))              
    bin_coeff <- choose((x + p_alpha - 1), x)                                 # should be the same
    pr_x <-  bin_coeff*(a^x)*(b^p_alpha)                                      # should be the same as dnbinom(x, size = p_alpha, p = b)
    return(lpdf = pr_x)
  })
}


## -- own CDF                                                ####

## pnbinom(x, size = p_alpha, p = b)

negbin_CDF <- function(x, p_lambda, p_alpha){
  sum(sapply(0:x, function(k){
    dnegbin_own(k, p_lambda = p_lambda, p_alpha = p_alpha)                     # pnbinom(x, size = p_alpha, p = b)
  }))
}   



## -- Log likelihood                                         ####
loglik_NB <- function(y,  p_alpha, p_lambda){
  loglambda <- log(p_lambda)
  loglambda[which(!is.finite(loglambda))] <- 0
  if(p_alpha == 0){
    theta_nb <- alpha_lgam <- alpha_log <- 0
  } else {
    theta_nb <- 1/p_alpha
    alpha_lgam <- lgamma(p_alpha)
    alpha_log <- log(p_alpha)
  }
  a <- lgamma(y+p_alpha)-alpha_lgam-lgamma(y+1)
  b <- y*loglambda+p_alpha*alpha_log
  c <- (p_alpha+y)*log(p_lambda+p_alpha)
  nb_loglik <- sum(a+b-c)
  return(nb_loglik)
}




## -- Gradient wrt beta                                      ####

## Assumes:
## -- p_lambda <- linkingFun(beta_iter = reg_coeff_iter, X = X)

gradient_NB <- function(p_lambda, p_alpha, X, y){
  p_theta <- 1/p_alpha
  n_observ <- nrow(X)
  if(!is.finite(p_theta)){p_theta <- 0} 
  residue_iter <- y - p_lambda       
  stacked_g <- sapply(1:n_observ, function(i){                                       
    i <- as.numeric(i)
    temp_g_i <- (residue_iter[i] / (1 + p_theta*p_lambda[i])) 
    temp_g_i*X[i,]
  })
  if(!is.matrix(stacked_g)){ g_gradient <- sum(stacked_g)
  } else { g_gradient <- as.matrix(apply(stacked_g,1,sum)) }
  return(list(gradient = g_gradient, stacked_g = stacked_g,
              pois_residue =  residue_iter)) 
}  





## -- Hessian wrt beta                                       ####

hessian_NB <- function(p_lambda, p_alpha, X, y){
  ## NOTE: unlike the poisson Hessian, y is an argument
  p_theta <- 1/p_alpha
  n_regressors <- ncol(X)
  n_observ <- nrow(X)
  Hessian_iter <- matrix(0L,nrow = n_regressors, ncol = n_regressors)
  for(i in 1:n_observ){
    x <- X[i,]
    x <- as.numeric(x)
    outprod_X <- (x %o% x)                                                              
    negbin_hess_num <- (1 + p_theta*y[i])*p_lambda[i]
    negbin_hess_den <- (1 + p_theta*p_lambda[i])^2
    Hessian_iter <- Hessian_iter - ((negbin_hess_num/negbin_hess_den))*outprod_X     
  }
  return(Hessian_iter)
}






## -- Gradient wrt alpha                                     ####
gradient_NB_alpha <- function(y, p_alpha, p_lambda){
  if(!is.matrix(y)) y <- as.matrix(y)
  n_observ <- nrow(y)
  test_1 <- which(!is.finite(p_alpha))
  if(length(test_1)>0){p_alpha[test_1] <- 0} 
  if(p_alpha == 0){
    alpha_log <-  0
    alpha_digam <- 0
  } else {
    alpha_log <-  log(p_alpha)
    alpha_digam <-  digamma(p_alpha)
  }  
  a <- digamma(p_alpha + y) - alpha_digam + (alpha_log + 1)
  b <- log(p_lambda + p_alpha)
  c <- (p_alpha + y)/(p_lambda + p_alpha)
  test_2 <- which(!is.finite(c))
  if(length(test_2)>0){c[test_2] <- 0} 
  d_dalpha <- sum(a - b - c)
  return(d_dalpha)
}




## -- Hessian wrt alpha                                      ####
Hessian_NB_alpha <- function(y, p_alpha, p_lambda){
  if(!is.finite(p_alpha)){p_alpha <- 0} 
  a <- trigamma(y + p_alpha) - trigamma(p_alpha)
  b <- (p_lambda^2 + p_alpha*y)
  c <- (p_alpha*(p_alpha + p_lambda)^2)
  dL2_dalpha2 <- sum( a + (b/c))
  return( dL2_dalpha2)
}


## -- Half-stepping                                          ####
myHalfStepping_NB <- function(Xnr, ynr, beta_iter, Hessian_iter, g_gradient, p_alpha_fix){
  accept_candidate <- FALSE
  converge_flag <- FALSE
  max_iter_halfstep <- 1000
  lambda_accepted <- linkingFun(beta_iter, Xnr)
  loglik_acceptd <- loglik_NB(ynr, p_lambda = lambda_accepted, p_alpha = p_alpha_fix)
  lambda_accepted <- linkingFun(beta_iter, Xnr)
  iter_halfstep <- 0
  stepsize <- 1
  target_iter <- qr.solve(Hessian_iter, g_gradient, tol = 1e-18)            
  while(!accept_candidate & (iter_halfstep < max_iter_halfstep)){
    temp_delta <- beta_iter - stepsize*target_iter
    lambda_iter <- linkingFun(temp_delta, Xnr)  
    iter_halfstep  <- iter_halfstep + 1
    loglik_iter <- loglik_NB(ynr, p_lambda = lambda_iter, p_alpha = p_alpha_fix)
    criterion_half_stepping <-  loglik_iter > loglik_acceptd
    if(!criterion_half_stepping){ stepsize <- stepsize/2 } else { 
      accept_candidate <- TRUE 
      target_iter <- stepsize*target_iter
    }
  } 
  if(iter_halfstep < max_iter_halfstep) converge_flag <- TRUE
  list(target_iter = target_iter, temp_delta = temp_delta, stepsize = stepsize,
       n_iter = iter_halfstep, converge_flag = converge_flag)
}
  
  
  
  
  
  


## -- Newton-Raphson, wrt Beta                               ####
NR_MLE_NB <- function(Xnr, ynr, p_alpha_fix, iter_max = 200, halfstep = TRUE){
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
  tol <- rep(1e-6, n_regressors)
  iter_stop <- FALSE
  hessian_is_singular <- FALSE
  while(iter_count <= iter_max & !iter_stop){
    iter_count <- iter_count + 1
    lambda_iter <- linkingFun(beta_iter, Xnr)   
    grad <- gradient_NB(p_lambda = lambda_iter, 
                        p_alpha = p_alpha_fix, 
                        X = Xnr, 
                        y = ynr)           
    g_gradient <- grad$gradient
    Hessian_iter <- hessian_NB(p_lambda = lambda_iter, 
                               p_alpha = p_alpha_fix, 
                               X = Xnr, 
                               y = ynr)        
    inv_test <- try(qr.solve(Hessian_iter, g_gradient), silent = TRUE)
    if(inherits(inv_test, "try-error")){
      hessian_is_singular <- TRUE
      iter_stop <- TRUE } 
    if(!hessian_is_singular){
      if(halfstep){
        HalfStep_proc <- myHalfStepping_NB(Xnr,
                                           ynr,
                                           beta_iter, 
                                           Hessian_iter, 
                                           g_gradient,
                                           p_alpha_fix)
        target_iter <- HalfStep_proc$target_iter
        temp_delta <-  HalfStep_proc$temp_delta
      } else {
        target_iter <- qr.solve(Hessian_iter, g_gradient)
        temp_delta <- beta_iter - target_iter
      }
    } else { warning("Hessian is singular")
      target_iter <- 0 
    }
    if(length(which(sqrt((target_iter)^2) < tol)) < n_regressors){
      beta_iter <- temp_delta } else {iter_stop <- TRUE}
  }
  return(list(coef_est = beta_iter, predict_y = lambda_iter, 
              Hessian = Hessian_iter, iter = iter_count  ))
}




## -- Newton-Raphson, wrt alpha                              ####
NR_MLE_NB_alpha <- function(Xnr, ynr, p_alpha_init, p_beta_fix, iter_max = 100){
  
  ## NOTE:
  ## For my initial data, which wouldn't converge I had to tweak max iter = 8 initally and then 23 when called recursively in the alternating procedure to replicate glm.nb
  
  if(!is.matrix(Xnr)){Xnr <- as.matrix(Xnr)}
  if(!is.matrix(ynr)){ynr <- as.matrix(ynr)}
  n_regressors <- ncol(Xnr)
  n_observ <- nrow(Xnr)
  p_lambda_0 <- linkingFun(p_beta_fix, Xnr)        # call own link function
  p_alpha_0 <- p_alpha_init                        # initial guess on alpha, to be provided externally (for flexi e.g. could be Hilbe's on first pass) 
  iter_count <- 0
  tol <- 1e-6
  iter_stop <- FALSE
  hessian_is_singular <- FALSE
  while(iter_count <= iter_max & !iter_stop){
    iter_count <- iter_count + 1
    if(iter_count == 1){
      p_alpha_iter <-  p_alpha_0
      p_lambda_iter <- p_lambda_0
    } 
    g_gradient <-  gradient_NB_alpha(y = ynr,
                                     p_alpha = p_alpha_iter, 
                                     p_lambda = p_lambda_iter
    )           
    Hessian_iter <- Hessian_NB_alpha(y = ynr,
                                     p_alpha = p_alpha_iter, 
                                     p_lambda = p_lambda_iter
    )                 
    inv_test <- !is.finite(g_gradient/Hessian_iter)                              # dealing with scalars now...
    if(inv_test){
      hessian_is_singular <- TRUE
      iter_stop <- TRUE } 
    if(!hessian_is_singular){
      target_iter <- g_gradient/Hessian_iter
    } else { warning("Hessian is singular")
      target_iter <- 0 
    }
    temp_delta <- p_alpha_iter - target_iter
    if(sqrt(target_iter^2) > tol){
      p_alpha_iter <- temp_delta 
      if(p_alpha_iter < 0) p_alpha_iter <- 0
    } else {iter_stop <- TRUE}
  }
  return(list(p_alpha = p_alpha_iter, predict_y = p_lambda_iter, 
              Hessian = Hessian_iter, iter = iter_count  ))
}



## -- Alternating fitting procedure                          ####
NegBReg_altern <- function(myData, target_feat =  target_feat){
  regressors_lablels <- colnames(myData)[-which(colnames(myData)==target_feat)] 
  y <- as.matrix(myData[,target_feat, drop = FALSE])        
  X <- cbind(1,as.matrix(myData[,regressors_lablels]))
  colnames(X) <- c("intercept",regressors_lablels)
  n_regressors <- ncol(X) 
  n_observ <- nrow(X)
  degrOfFreed <- n_observ - n_regressors                  
  my_Pois_NR <- NR_MLE(Xnr = X, ynr = y)
  temp_lambda <- my_Pois_NR$predict_y                 
  temp_beta <- my_Pois_NR$coef_est        
  temp_alpha <- n_observ/sum(((y/temp_lambda)-1)^2)
  alpha_1 <- NR_MLE_NB_alpha(Xnr = X, ynr = y, p_alpha_init = temp_alpha, p_beta_fix = temp_beta)
  if(alpha_1$p_alpha < 0){ chosen_alpha <- 0} else { chosen_alpha <- alpha_1$p_alpha }
  alpha_iter_lst  <- list()
  max_iter3 <- big_alpha <- 100 
  tol_gap3 <- 1e-6
  gap3 <- 1                                                               
  hessian_is_singular <- FALSE
  d1 <- sqrt(2*max(1,degrOfFreed))                                                                                       
  d2 <- 1
  Lm <- loglik_NB(y = y, p_alpha = chosen_alpha, p_lambda = temp_lambda)
  Lm0 <- Lm + 2*d1                                                                             
  stopping_criterion <- (abs(Lm0 - Lm)/d1 + abs(gap3)/d2)        # as in MASS::glm.nb                                    
  tol_stop <- 1e-4
  iter_count3 <- 1
  alpha_iter_lst[[iter_count3]] <- chosen_alpha 
  while ((iter_count3 < max_iter3) && (abs(gap3) > tol_gap3) && (stopping_criterion > tol_stop )) {
    alpha_iter <- as.numeric(alpha_iter_lst[[iter_count3]])                                        
    negbin_fit_iter <- NR_MLE_NB(Xnr = X, ynr = y, p_alpha_fix = alpha_iter)
    lambda_iter <- negbin_fit_iter$predict_y           
    beta_iter <- negbin_fit_iter$coef_est
    temp_alpha <- n_observ/sum(((y/lambda_iter)-1)^2)
    alpha_refined <- NR_MLE_NB_alpha(Xnr = X, ynr = y, p_alpha_init = temp_alpha, p_beta_fix = beta_iter)
    if(alpha_refined$p_alpha < 0){ chosen_alpha <- 0} else { chosen_alpha <- alpha_refined$p_alpha }
    gap3 <- chosen_alpha - alpha_iter
    Lm0 <- Lm
    Lm <-  loglik_NB(y = y, p_alpha = chosen_alpha, p_lambda = lambda_iter)
    stopping_criterion <- (abs(Lm0 - Lm)/d1 + abs(gap3)/d2) 
    iter_count3 <- iter_count3 + 1
    alpha_iter_lst[[iter_count3]] <- chosen_alpha 
  }
  return(list(p_alpha_MLE=chosen_alpha, p_lambda_MLE=lambda_iter, p_beta_MLE=beta_iter, log_likelihood=Lm, n_iter=iter_count3 ))
}


## -- Fisher Info matrix elements: wrt beta                  ####

FI_Beta <- function(p_lambda, p_alpha, X, y){
  p_theta <- 1/p_alpha
  n_regressors <- ncol(X)
  n_observ <- nrow(X)
  FI_Beta_iter <- matrix(0L,nrow = n_regressors, ncol = n_regressors)
  for(i in 1:n_observ){
    x <- X[i,]
    x <- as.numeric(x)
    outprod_X <- (x %o% x)                                                              
    FI_Beta_num <- p_lambda[i]
    FI_Beta_den <- 1 + p_theta*p_lambda[i]
    FI_Beta_iter <- FI_Beta_iter - ((FI_Beta_num/FI_Beta_den))*outprod_X     
  }
  FI_Beta_iter <- -1*FI_Beta_iter
  return(FI_Beta_iter)
}



## -- Scoring wrt Beta                                       ####
SCOR_MLE_NB <- function(Xnr, ynr, p_alpha_fix, iter_max = 200, halfstep = TRUE){
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
  tol <- rep(1e-6, n_regressors)
  iter_stop <- FALSE
  FI_is_singular <- FALSE
  while(iter_count <= iter_max & !iter_stop){
    iter_count <- iter_count + 1
    lambda_iter <- linkingFun(beta_iter, Xnr)   
    grad <- gradient_NB(p_lambda = lambda_iter, 
                        p_alpha = p_alpha_fix, 
                        X = Xnr, y = ynr)           
    g_gradient <- grad$gradient
    FI_iter <- -1*FI_Beta(p_lambda = lambda_iter, p_alpha = p_alpha_fix, X = Xnr, y = ynr)        
    inv_test <- try(qr.solve(FI_iter, g_gradient), silent = TRUE)
    if(inherits(inv_test, "try-error")){
      FI_is_singular <- TRUE
      iter_stop <- TRUE 
    } 
    if(!FI_is_singular){
      if(halfstep){
        HalfStep_proc <- myHalfStepping_NB(Xnr,  ynr, beta_iter, FI_iter, g_gradient, p_alpha_fix)
        target_iter <- HalfStep_proc$target_iter
        temp_delta <-  HalfStep_proc$temp_delta
      } else {
        target_iter <- qr.solve(FI_iter, g_gradient, tol = 1e-18)
        temp_delta <- beta_iter - target_iter
      }
    } else { warning("Hessian is singular")
      target_iter <- 0 
    }
    if(length(which(sqrt((target_iter)^2) < tol)) < n_regressors){
      beta_iter <- temp_delta } else {iter_stop <- TRUE}
  }
  return(list(coef_est = beta_iter, predict_y = lambda_iter, 
              Hessian = FI_iter, iter = iter_count  ))
}


##                                                           ####
## Fisher Info: wrt alpha                                    ####

## With the approximation for Expected Gamma from Yu and Mousavi 2024

FI_alpha <- function(p_lambda, p_alpha, X, M =  20){
  fi_a2 <- sapply(1:length(p_lambda), function(i){
    x_i <- X[i,]
    p_lambda_i <- p_lambda[i]
    fi_a1 <-  sapply(0:M, function(j){
      nb_CDF_j <- negbin_CDF(j, p_lambda = p_lambda_i, p_alpha)
      (1-nb_CDF_j )/(j+p_alpha)^2
    })
    sum(fi_a1) - (p_lambda_i/(p_alpha*(p_alpha + p_lambda_i)))
  })
  sum(fi_a2)
}


## -- full alternating procedure                             ####
test_NB <- NegBReg_altern(myData = myData, 
              target_feat =  target_feat)

## -- Variance for beta                                      ####
Info_NB <- FI_Beta(p_lambda = linkingFun(test_NB$p_beta_MLE, X),
                   p_alpha = test_NB$p_alpha_MLE, 
                   X = X, y = y)
vcov_NB <- qr.solve(Info_NB, tol = 1e-18)

## -- variance for alpha                                     ####
p_lambda <- test_NB$p_lambda_MLE
p_alpha <- test_NB$p_alpha_MLE
est_FI_alpha <- FI_alpha(p_lambda, p_alpha, X, M = 50)
var_alpha_NB <- 1/est_FI_alpha

## -- Observed info                                          ####
Info_NB_obs <- -1*hessian_NB(p_lambda = test_NB$p_lambda_MLE, p_alpha = test_NB$p_alpha_MLE, X = X, y = y)
vcov_NB_obs <-  qr.solve(Info_NB_obs, tol = 1e-18)
SE_model_NB_obs <- sqrt(diag(vcov_NB_obs))                                      

Info_NB_obs_alpha <- -1*Hessian_NB_alpha(y=y, p_alpha = test_NB$p_alpha_MLE, p_lambda = test_NB$p_lambda_MLE)
var_NB_obs_alpha <- 1/Info_NB_obs_alpha
SE_mode_NB_obs_alpha <- sqrt(var_NB_obs_alpha)

## -- SE expected vs observed info                           ####
SE_model_NB_expected <- sqrt(diag(vcov_NB))
SE_model_NB_expected_alpha <- sqrt(var_alpha_NB)





##                                             ####
## Comparison with pre-built                   ####
## --- MASS:glm.nb                             ####
obj_NB <- MASS::glm.nb(g_DISRUPTIONS ~.,
                       data = myData)

obj_NB$coefficients
obj_NB$theta
1/obj_NB$theta           # dispersion parameter, to be compared with GAMLSS                               



## --- GAMLSS                                  ####
obj_NB_gam <- gamlss::gamlss(g_DISRUPTIONS ~.,
                             family = "NBI",
                             data = myData) 
obj_NB_gam$mu.coefficients
exp(obj_NB_gam$sigma.coefficients)                # dispersion parameter
alpha_gamlss <- 1/exp(obj_NB_gam$sigma.coefficients)

##
summary(obj_NB_gam)[1:5,]

## extract standar errors
var_gamlss <- vcov(obj_NB_gam)
var_gamlss_coeff <- diag(var_gamlss[1:(nrow(var_gamlss)-1), 1:(ncol(var_gamlss)-1)])
var_gamlss_logTheta <- diag(var_gamlss)[nrow(var_gamlss)]
SE_gamlss <- sqrt(var_gamlss_coeff)


## --- Table comparative                       ####
tab_pre_built <-  round(summary(obj_NB)$coefficients[,1:2],5)
tab_pre_built_GAMLSS <- round(cbind.data.frame(obj_NB_gam$mu.coefficients, SE_gamlss),5)
tab_fromScratch <- round(cbind.data.frame(test_NB$p_beta_MLE, SE_model_NB_expected, SE_model_NB_obs),5)

tab_pre_built_alpha <- round(cbind.data.frame(summary(obj_NB)[17], summary(obj_NB)[18]),5)
tab_pre_built_alpha_GAMLSS <- cbind.data.frame(round(alpha_gamlss,5), "n.a.")
tab_fromScratch_alpha <- round(cbind.data.frame(test_NB$p_alpha_MLE, SE_model_NB_expected_alpha, SE_mode_NB_obs_alpha),5)

colnames(tab_fromScratch)  <- c("Estimate", "SE expected", "SE observed")
colnames(tab_pre_built) <- paste0(c("glmnb_"),colnames(tab_pre_built))
colnames(tab_pre_built_GAMLSS) <- paste0(c("gamlss_coeff", "gamlss_SE"))

a <- cbind.data.frame(tab_fromScratch, tab_pre_built, tab_pre_built_GAMLSS)
b <- cbind.data.frame(tab_fromScratch_alpha, tab_pre_built_alpha, tab_pre_built_alpha_GAMLSS)
colnames(b) <- colnames(a)
rownames(b) <- "alpha"
myTab4viz <- rbind.data.frame(a, b)




