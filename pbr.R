## matrix predictors: x_train; x_test
## scalar response: y_train; y_test
## prior_type:  SSL (default) or Uniform
## b0: hyperparameter of the mixture Laplace distribution (spike component in the SSL prior)
## b1: hyperparameter of the mixture Laplace distribution  (slab component in the SSL prior)
## rho: hyparameter of the prior of the basis coefficients 
## K: number of ridge functions 
## nwarmup: number of warmup iterations
## nsample: number of sampling iterations (including warmup and post warmup), nsample > nwarmup 
## init_params: list of initial parameter values 
## prior_params: list of parameters in the prior distributions 
## pred: make predictions or not 


pbr <- function(x_train, z_train = NULL, y_train, prior_type = "SSL",  x_test = NULL, z_test = NULL,  y_test = NULL, K, nknots, rho = 0, b0 = 0.1, b1 = 1, prior_params, init_params, nwarmup, nsample, pred = FALSE){
  
  set.seed(123)
  
  # parameters for the prior distribution 
  mu_0 <- prior_params$mu_0; # prior for mu (Normal)
  sigma_0 <- prior_params$sigma_0; # prior for mu (Normal)
  alpha_1 <- prior_params$alpha_1; # prior for sigma (IG)
  beta_1 <- prior_params$beta_1  # prior for sigma (IG)
  alpha_2 <- prior_params$alpha_2; # prior for "w" ("a" in the code) (Beta, mixing proportion, Spike-and-slab)
  beta_2 <- prior_params$beta_2; # prior for "w" ("a" in the code) (Beta, mixing proportion, Spike-and-slab)

  
  
  n_train <- dim(x_train)[1]
  n_test <- dim(x_test)[1]
  p <-  dim(x_train)[2]
  lambda_0 <- rep(10000,K)
  
  J <- nknots + 2  # number of basis functions
  acc <- matrix(NA, nsample, K) # acceptance rate 
  acc[1,] <- rep(1,K)
  
  # save samples and initialize
  mu_s <- rep(NA, nsample); 
  mu_s[1] <- init_params$mu; 
  sigma2_s <- rep(NA, nsample);  
  sigma2_s[1] <- init_params$sigma2;
  
  beta_s <- array(NA, dim = c(nsample, J, K)) # intercept = true; save the basis coefficients 
  a_s <- array(NA, dim = c(nsample, p-1, K)) # save the angles 
  
  if(!is.null(z_train)){
    beta_s_z <- array(NA, dim = c(nsample, J, ncol(z_train))) # might be some empty columns for characterized z
    coef_s_z <- matrix(NA, nsample, ncol(z_train)) 
  }
  
  for(k in 1:K){
    beta_s[1,,k] <- rep(0, J)
  }
  
  gamma_s <-  array(NA, dim = c(nsample, p, K)) # save the gammas 
  knots_s <-  array(NA, dim = c(nsample, nknots + 2, K)) # save the knots 
  
  if(prior_type == "SSL"){
    mix_a_s <- array(NA, dim = c(nsample,K)) # mixing proportions
    m_s <- array(NA, dim = c(nsample,p-1, K)) # save the allocation indicator of each angle
  }
  
  
  if(!is.null(init_params$gammas)){
    print("init gamma provided!")
    gamma_s[1,,] <- init_params$gammas
  }else{
    for(k in 1:K){
      mu_gamma_0 <- rnorm(p)
      mu_gamma_0 <- mu_gamma_0/sqrt(sum(mu_gamma_0^2))  # some randomly initialized directions
      gamma_s[1,,k] <- rvmf(1, mu = mu_gamma_0, k = lambda_0[k])
    }
  }
  
  if(prior_type == "SSL"){
    if(!is.null(init_params$mix_a)){ 
      mix_a_s[1,] <-  rep(init_params$mix_a, K)
    }else{
      mix_a_s[1,] <-  rep(0.5,K)
    }
  }
  
  for(k in 1:K){
    a_prop <- c_to_p(gamma_s[1,,k]) # proposed angles 
    a_s[1,,k] <-  a_prop 
  }
  
  # save the estimated components at the current iteration
  comp_save <- matrix(0, nrow = n_train, ncol = K)
  g_center_save <- matrix(0, nrow = nsample, ncol = K) # save the constant to center each g
  
  if(!is.null(z_train)){
    comp_save_z <- matrix(0, nrow = n_train, ncol = ncol(z_train))
    g_center_save_z <-  matrix(0, nrow = nsample, ncol = ncol(z_train))
    knots_s_z <-  array(NA, dim = c(nsample, nknots + 2, ncol(z_train)))
  }
  
  for(i in 2:nsample){
    
    if(i%%100 == 0){
      print(paste(i, "samples!"))
    }
    
    for(k in 1:K){  
      
      gamma_pre <- gamma_s[i-1,,k]
      a_pre <- a_s[i-1,,k]
      
      # Step 1.1: sample the allocation 
      if(prior_type == "SSL"){
        pp1 <- (mix_a_s[i-1,k]*exp(-abs(a_pre)/b0))/b0
        pp2 <- ((1- mix_a_s[i-1,k])* exp(-abs(a_pre)/b1))/b1
        m_s[i,,k] <- rbinom(p-1,1, pp1/(pp1+pp2)) 
        m_pre <- m_s[i,,k] 
      }
      
      # Step 1.2: sample the mixing proportion 
      if(prior_type == "SSL"){
        mix_a_s[i,k] <- rbeta(1, sum(m_s[i,,k])+alpha_2,  sum(1-m_s[i,,k])+beta_2)         
        mix_a_pre <-  mix_a_s[i,k] 
      }else{
        mix_a_pre <- m_pre <- NULL
      }
      
      # Step 1.3 
      # compute the residuals 
      if(!is.null(z_train)){
        rr <- switch(as.character(K), "1" = y_train  - apply(comp_save_z, 1, sum) -  mu_s[i-1],
                     "2" = y_train - comp_save[, -k] -  apply(comp_save_z, 1, sum) - mu_s[i-1], 
                     y_train - apply(comp_save[, -k], 1, sum) - apply(comp_save_z, 1, sum) -  mu_s[i-1])
      }else{
        rr <- switch(as.character(K), "1" = y_train  -  mu_s[i-1],
                     "2" = y_train - comp_save[, -k] -  mu_s[i-1], 
                     y_train - apply(comp_save[, -k], 1, sum) -  mu_s[i-1])
      }
      
      mh_tmp <- gamma.mh(x_train, prior_type, gamma_pre, mix_a_pre, m_pre, nknots, alpha_1, beta_1, alpha_2, beta_2, b0, b1, rho, lambda_0[k], rr)
      acc[i,k] <-  mh_tmp$acc # binary, accept or not 
      
      
      if(mh_tmp$acc){ # found an accepted sample
        gamma_s[i,,k] <-  mh_tmp$gamma
        a_s[i,,k] <-  mh_tmp$a
      }else{
        gamma_s[i,,k] <-  gamma_s[i-1,,k]
        a_s[i,,k] <-  a_s[i-1,,k] 
      }
      
      
      index_tmp <- apply(x_train, 1, function(u){t(gamma_s[i,,k])%*%u%*%gamma_s[i,,k]}) # compute the index 
      
      knots_tmp <- quantile(index_tmp, seq(0,1,length.out = nknots+2)) # include the boundary
      knots <-    knots_tmp

      B_gamma <- ns(index_tmp,  Boundary.knots =c(knots[1], knots[nknots +2]), knots = knots[2:(2+nknots-1)],  intercept = TRUE)
      B_gamma_T_B_gamma <- t(B_gamma)%*% B_gamma 
      
      beta_s[i,,k] <- solve(B_gamma_T_B_gamma + rho* diag(J))%*%t(B_gamma)%*%rr
      knots_s[i,,k] <- knots
      comp_save[,k] <-  B_gamma%*%beta_s[i,,k]
      g_center_save[i,k] <- mean(comp_save[,k])
      comp_save[,k] <-  comp_save[,k] - mean(comp_save[,k]) # center
      
        if(i %% 100 == 0){
          if(mean(acc[(i-100+1):i,k]) < 0.2){
            if(lambda_0[k] <100000){
              lambda_0[k] <- lambda_0[k]*1.1
            }
          }
          if(mean(acc[(i-100+1):i,k]) > 0.4){
            lambda_0[k] <- lambda_0[k]/1.1
          }
      }
    }
    
    # -- optional step to fit Z 
    if(!is.null(z_train)){
      pz <- ncol(z_train)
      for(jj in 1:pz){
        
        rr_z <- switch(as.character(pz), "1" = y_train  - apply(comp_save, 1, sum) -  mu_s[i-1],
                       "2" = y_train - comp_save_z[, -jj] -  apply(comp_save, 1, sum) - mu_s[i-1], 
                       y_train - apply(comp_save_z[, -jj], 1, sum) - apply(comp_save, 1, sum) -  mu_s[i-1])
        
          knots_tmp <- quantile(z_train[,jj], seq(0,1,length.out = nknots+2)) # include the boundary
          knots <-    knots_tmp
          B_gamma <- ns(z_train[,jj],  Boundary.knots =c(knots[1], knots[nknots +2]), knots = knots[2:(2+nknots-1)],  intercept = TRUE)
          B_gamma_T_B_gamma <- t(B_gamma)%*% B_gamma 
          beta_s_z[i,,jj] <- solve(B_gamma_T_B_gamma + rho* diag(J))%*%t(B_gamma)%*%rr_z
          knots_s_z[i,,jj] <- knots
          comp_save_z[,jj] <-  B_gamma%*%beta_s_z[i,,jj]
          g_center_save_z[i,jj] <- mean(comp_save_z[,jj])
          comp_save_z[,jj] <-  comp_save_z[,jj] - mean(comp_save_z[,jj]) # center 
        }
    }
    
    # -- Step 2
    ri <- y_train - apply(comp_save, 1, sum)
    alpha_sigma2_post <- n_train/2 + alpha_1
    beta_sigma2_post <-  sum((ri - mu_s[i-1])^2)/2 + beta_1
    sigma2_s[i] <- 1/rgamma(1, shape = alpha_sigma2_post,  beta_sigma2_post) 
    
    # -- Step 3
    sigma2_tilde <- 1/(1/(sigma_0^2) + n_train/sigma2_s[i])
    mu_tilde <- sigma2_tilde*(mu_0/(sigma_0^2) + sum(ri)/sigma2_s[i])
    mu_s[i] <- rnorm(1,  mu_tilde, sd =   sqrt(sigma2_tilde))
  }
  
  obj_to_return <- list( acc = acc, g_center_save =  g_center_save, mu_s = mu_s, beta_s = beta_s, gamma_s = gamma_s,  knots_s = knots_s, 
                         sigma2_s =  sigma2_s, nwarmup = nwarmup,  nsample = nsample)
  
  if(!is.null(z_train)){
    obj_to_return <- c(obj_to_return, list(knots_s_z  = knots_s_z, g_center_save_z = g_center_save_z, 
                                            comp_save_z = comp_save_z,  beta_s_z =  beta_s_z))
  }
  if(prior_type == "SSL"){
    obj_to_return <- c(obj_to_return, list(mix_a_s =  mix_a_s,  m_s =  m_s))
  }
  
  if(pred == TRUE){
    pred_train <- pbr.predict(obj_to_return, x_new = x_train, z_new = z_train, y_new = y_train)
    pred_test <- pbr.predict(obj_to_return, x_new = x_test, z_new = z_test, y_new = y_test)
    obj_to_return <- c(obj_to_return, list(pred_train = pred_train, pred_test = pred_test))
  }
  return(obj_to_return)
}


c_to_p <- function(beta){
  
  p <- length(beta)
  theta <- rep(NA, p-1)
  tmp <- beta[p]
  for(j in (p-1):1){
    theta[j] <- atan(beta[j]/tmp)
    tmp <- tmp/cos(theta[j])
  }
  return(theta)
}

p_to_c <- function(theta){
  
  p <- length(theta)+1
  beta <- rep(NA, p)
  beta[1] <-sin(theta[1])
  tmp1 <- cumprod(cos(theta))
  
  for(j in 2:(p-1)){
    beta[j] <- sin(theta[j])*tmp1[j-1]
  }
  
  beta[p] <- tail(tmp1,1)
  
  return(beta)
}


generate.prop.von <- function(gamma_pre, lambda_0){
  
  gamma_prop <-  c(rvmf(1, mu = gamma_pre, k = lambda_0))
  return(gamma_prop)
  
}




gamma.mh <- function(x_train, prior_type,  gamma_pre,  mix_a_pre, m_pre,  nknots, alpha_1, beta_1, alpha_2, beta_2, b_0, b_1, rho, lambda_0, rr){
  
  n_train <- dim(x_train)[1]
  
  p <- dim(x_train)[2]
  J <- nknots + 2
  a_pre <- c_to_p(gamma_pre) # angles from previous iteration 
  
  cal_log_prop <- function(x_train, prior_type, a_tmp,  nknots, rr){
    
    gamma_tmp <- p_to_c(a_tmp)
    
    index_tmp <- apply(x_train, 1, function(u){t(gamma_tmp)%*%u%*%gamma_tmp}) # compute the index 
    
    knots_tmp <- quantile(index_tmp, seq(0,1,length.out = nknots+2)) # include the boundary
    knots <-    knots_tmp
    
    B_gamma <- ns(index_tmp,  Boundary.knots =c(knots[1], knots[nknots +2]), knots = knots[2:(2+nknots-1)],  intercept = TRUE)
    
    B_gamma_T_B_gamma <- t(B_gamma)%*% B_gamma 
    Sigma_rho <- solve(B_gamma_T_B_gamma + rho* diag(J))
    Sigma_0 <-  solve(B_gamma_T_B_gamma)
    inv_Sigma_0 <- B_gamma_T_B_gamma
    
    S_theta <- t(rr)%*%rr - t(rr)%*%B_gamma %*% (0.5 * Sigma_0 + Sigma_rho  - (Sigma_rho%*%inv_Sigma_0%*% Sigma_rho)/2) %*% t(B_gamma)%*%rr
    
    if(prior_type == "Uniform"){
      log_D_theta <-  -(alpha_1 + n_train/2) * log(S_theta + 2*beta_1)  - sum((p-1-(1:(p-2)))*log(abs(cos(a_tmp[1:(p-2)]))))
    }
    
    if(prior_type == "SSL"){
      tmp_ss <- m_pre * (-abs(a_tmp)/b_0) + (1 - m_pre)*(-abs(a_tmp)/b_1)
      log_D_theta <-  -(alpha_1 + n_train/2) * log(S_theta + 2*beta_1) + sum(tmp_ss) -  sum((p-1-(1:(p-2)))*log(abs(cos(a_tmp[1:(p-2)]))))  
    }
    
    return(log_D_theta)
  }
  
  acc  <- 0  
  gamma_prop <-  generate.prop.von(gamma_pre, lambda_0)
  a_prop <- c_to_p(gamma_prop)
  acc_p <- min(1, exp(cal_log_prop(x_train, prior_type, a_prop,  nknots,  rr) -
                        cal_log_prop(x_train, prior_type, a_pre,  nknots, rr)))
  
  prob_tmp <- runif(1) 
  if(prob_tmp <= acc_p){
    acc <- 1
  }
  
  if(acc == 1){
    return(list(acc = 1, gamma = gamma_prop, a = a_prop, a_prop= a_prop))
  }else{
    return(list(acc = 0))
  }
}

# make predictions for the test data 
pbr.predict <- function(bppr_obj, x_new, z_new = NULL, y_new = NULL){
  
  nsamples <- length(bppr_obj$mu_s)
  nwarmup <- bppr_obj$nwarmup
  n_new <- dim(x_new)[1]
  K <- dim(bppr_obj$gamma_s)[3]
  g_pred_s <- array(NA, dim = c(nsamples -  nwarmup, K, n_new))
  mu_pred_s <- rep(NA, nsamples -  nwarmup)
  nknots_all <- dim(bppr_obj$knots_s)[2]
  
  if(!is.null(z_new)){
    g_pred_s_z <- array(NA, dim = c(nsamples -  nwarmup, ncol(z_new), n_new))
    knots_s_z  <- bppr_obj$knots_s_z
    g_center_save_z <- bppr_obj$g_center_save_z 
    beta_s_z <-  bppr_obj$beta_s_z
    coef_s_z <- bppr_obj$coef_s_z
  }
  
  for(i in (nwarmup+1):nsamples){
    
    mu_pred_s[i- nwarmup] <-  bppr_obj$mu_s[i]
    
    for(k in 1:K){
      gamma_pre <- bppr_obj$gamma_s[i,,k]
      index_tmp <-  apply(x_new, 1, function(u){t(gamma_pre)%*%u%*%gamma_pre})
      B_gamma <- ns(index_tmp,  Boundary.knots =c(bppr_obj$knots_s[i,1,k], bppr_obj$knots_s[i,nknots_all,k]), knots =  bppr_obj$knots_s[i,,k][2:(nknots_all-1)],  intercept = TRUE)
      g_pred_s[i - nwarmup,k, ]  <-   B_gamma%*% bppr_obj$beta_s[i,,k ] -  bppr_obj$g_center_save[i,k]
    }
    
    if(!is.null(z_new)){
      for(jj in 1:ncol(z_new)){
          B_gamma <- ns(z_new[,jj],  Boundary.knots =c(bppr_obj$knots_s_z[i,1,jj], bppr_obj$knots_s_z[i,nknots_all,k]), knots =  bppr_obj$knots_s_z[i,,k][2:(nknots_all-1)],  intercept = TRUE)
          g_pred_s_z[i - nwarmup,jj, ]  <-   B_gamma%*% bppr_obj$beta_s_z[i,,jj] -  bppr_obj$g_center_save_z[i,jj]
        }
    }
  }
  
  if(!is.null(y_new)){
    if(!is.null(z_new)){
      mspe <- mean((y_new - apply(apply(g_pred_s, c(1,3), sum), 2, mean) - apply(apply(g_pred_s_z, c(1,3), sum), 2, mean)  - mean( mu_pred_s))^2)
    }else{
      mspe <- mean((y_new - apply(apply(g_pred_s, c(1,3), sum), 2, mean) - mean( mu_pred_s))^2)
    }
  }
  
  if(!is.null(z_new)){
    return( list(mu_pred_s = mu_pred_s,  g_pred_s =  g_pred_s, g_pred_s_z = g_pred_s_z, mspe  =   mspe ))
  }else{
    return( list(mu_pred_s = mu_pred_s,  g_pred_s =  g_pred_s, mspe  =   mspe ))
  }
}


## for initialization
init_ppr <- function(x_train, y_train, K){
  
  p <- dim(x_train)[2]
  tmp <-  t(apply(x_train, c(1), function(u)(u[lower.tri(u, diag = TRUE)])))
  
  fit_ppr <- ppr(tmp, y_train, nterms = K)
  coef_init <- matrix(NA, p, K)
  for(i in 1:K){
    coef_tmp <- matrix(NA, dim(x_train)[2],  dim(x_train)[2])
    coef_tmp[lower.tri(coef_tmp, diag = TRUE)] <- fit_ppr$alpha[,i]
    for(j in 1: (p-1)){
      for(k in (j+1):p){
        coef_tmp[j,k] <- coef_tmp[k,j]
      }
    }
    for(ii in 1:p){
      for(jj in 1:p){
        if(ii!=jj){
          coef_tmp[ii,jj] <- coef_tmp[ii,jj]/2
        }
      }
    }
    coef_init[,i] <- eigen(coef_tmp)$vectors[,1]
  }
  
  return(coef_init)
}

