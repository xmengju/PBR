# sparse version of bppr, use spherical coordianates 
# K is the number of additive components 
# A sparse version of bppr

sbppr <- function(x_train, y_train, proposal_type = "von-Mises", prior_type = "uniform", knot_type = "quantile", lambda0_type = "no",  x_test = NULL, y_test = NULL, K, nknots, rho = 0, b = 0.1, b0 = 0.1, b1 = 1, prior_params, init_params, nwarmup, nsample, pred = FALSE, plot_opt = FALSE){
  
  set.seed(123)
  # parameters for the prior distribution 
  mu_0 <- prior_params$mu_0;
  sigma_0 <- prior_params$sigma_0;
  alpha_1 <- prior_params$alpha_1; # for sigma
  beta_1 <- prior_params$beta_1 # for sigma
  alpha_2 <- prior_params$alpha_2; # for "b" (LaplaceH or Laplace) or for "a" (Spike-and-slab)
  beta_2 <- prior_params$beta_2; # for "b" (LaplaceH or Laplace) or for "a" (Spike-and-slab)
  
  
  n_train <- dim(x_train)[1]
  n_test <- dim(x_test)[1]
  p <-  dim(x_train)[2]
  #lambda_0_begin <- 200
  #lambda_0_end <- 10000
  lambda_0 <- rep(10000,K)
  
  J <- nknots + 2  # number of basis functions 
  acc <- matrix(NA, nsample, K)
  acc[1,] <- rep(1,K)
  
  # save samples and initialize
  mu_s <- rep(NA, nsample); 
  mu_s[1] <- init_params$mu; 
  sigma2_s <- rep(NA, nsample);  
  sigma2_s[1] <- init_params$sigma2;
  
  beta_s <- array(NA, dim = c(nsample, J, K)) # intercept = true
  a_s <- array(NA, dim = c(nsample, p-1, K)) # save the angles 

  for(k in 1:K){
    beta_s[1,,k] <- rep(0, J)
  }
  
  gamma_s <-  array(NA, dim = c(nsample, p, K))
  b_s <-  array(NA, dim = c(nsample,K))
  knots_s <-  array(NA, dim = c(nsample, nknots + 2, K))
  

  if(prior_type == "spike_slab"){
    mix_a_s <- array(NA, dim = c(nsample,K)) # mixing proportions
    m_s <- array(NA, dim = c(nsample,p-1, K)) # save the membership of each angle
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
  
  if(prior_type == "spike_slab"){
    if(!is.null(init_params$mix_a)){ 
        mix_a_s[1,] <-  rep(init_params$mix_a,K)
    }else{
        mix_a_s[1,] <-  rep(0.5,K)
    }
  }
  
  for(k in 1:K){
    a_prop <- c_to_p(gamma_s[1,,k])
    a_s[1,,k] <-  a_prop
   }
  
  # save the estimated components at the current iteration
  comp_save <- matrix(0, nrow = n_train, ncol = K)
  g_center_save <- matrix(0, nrow = nsample, ncol = K) # the constant to center each g
  
  
  for(i in 2:nsample){
    
    #if(i < nwarmup){
       #  lambda_0 <- lambda_0_begin + ((i)/nwarmup)*(lambda_0_end - lambda_0_begin)
    #}
    if(i%%100 == 0){
      print(paste(i, "samples!"))
    }
    
    for(k in 1:K){  
      
      if(prior_type == "laplace_hyper"){
        if(i == 2){
          b_s[i,k] <- 1/rgamma(1,  alpha_2 + p - 1, (beta_2 + sum(abs(a_s[i-1,,k]))))
        }else{
          if(mh_tmp$acc){ # from the previous iterations
            b_s[i,k] <- 1/rgamma(1,  alpha_2 + p - 1, (beta_2 + sum(abs(a_s[i-1,,k]))))
          }else{
            b_s[i,k] <-  b_s[i-1,k] 
          }
        }
      }
      
      b_pre  <- b_s[i,k] # previous to sampling gamma
      gamma_pre <- gamma_s[i-1,,k]
      a_pre <- a_s[i-1,,k]
      
      # Step 1.1
      if(prior_type == "spike_slab"){
        pp1 <- (mix_a_s[i-1,k]*exp(-abs(a_pre)/b0))/b0
        pp2 <- ((1- mix_a_s[i-1,k])* exp(-abs(a_pre)/b1))/b1
        m_s[i,,k] <- rbinom(p-1,1, pp1/(pp1+pp2)) 
      }
      
      # Step 1.2 
      if(prior_type == "spike_slab"){
        mix_a_s[i,k] <- rbeta(1, sum(m_s[i,,k])+alpha_2,  sum(1-m_s[i,,k])+beta_2)         
        mix_a_pre <-  mix_a_s[i,k] 
        m_pre <- m_s[i,,k] 
      }else{
        mix_a_pre <- m_pre <- NULL
      }
        
      
      # Step 1.3 
      # compute the residuals 
      rr <- switch(as.character(K), "1" = y_train  -  mu_s[i-1],
                   "2" = y_train - comp_save[, -k] -  mu_s[i-1], 
                   y_train - apply(comp_save[, -k], 1, sum) -  mu_s[i-1])
      
      mh_tmp <- gamma.mh.s(x_train, proposal_type, prior_type, knot_type, gamma_pre,  b_pre, mix_a_pre, m_pre, nknots,  alpha_1, beta_1, alpha_2, beta_2, b0, b1, rho,  b, lambda_0[k], rr)
      acc[i,k] <-  mh_tmp$acc

      
      if(mh_tmp$acc){ # found an accepted sample
        gamma_s[i,,k] <-  mh_tmp$gamma
        a_s[i,,k] <-  mh_tmp$a
      }else{
        gamma_s[i,,k] <-  gamma_s[i-1,,k]
        a_s[i,,k] <-  a_s[i-1,,k] 
      }

      
      index_tmp <- apply(x_train, 1, function(u){t(gamma_s[i,,k])%*%u%*%gamma_s[i,,k]}) # compute the index 
      
      
      if(knot_type== "quantile_adj"){
        knots_tmp <- quantile(index_tmp, seq(0.05,0.95,length.out = nknots)) # include the boundary
        knots_delta <-0.05*(max(index_tmp)  -  min(index_tmp)) 
        knots <-  c(min(index_tmp)-  knots_delta, knots_tmp, max(index_tmp) +  knots_delta)
      }
      
      if(knot_type== "quantile"){
        knots_tmp <- quantile(index_tmp, seq(0,1,length.out = nknots+2)) # include the boundary
        knots <-    knots_tmp
      }
      
      if(knot_type == "even"){
        knots_tmp <- seq(min(index_tmp), max(index_tmp),length.out = nknots + 2) # include the boundary
        knots_delta <-0.05*(knots_tmp[nknots]  -  knots_tmp[1]) 
        knots <-  c(knots_tmp[1] -  knots_delta, knots_tmp[2:(nknots+1)], knots_tmp[nknots+2] +  knots_delta)
      }
      
      
      B_gamma <- ns(index_tmp,  Boundary.knots =c(knots[1], knots[nknots +2]), knots = knots[2:(2+nknots-1)],  intercept = TRUE)
      
      #B_gamma <- bSpline(index_tmp, Boundary.knots = c(knots[1],  knots[nknots + 2]), knots = knots[2:(2+nknots-1)], degree = 3, intercept = TRUE)
      B_gamma_T_B_gamma <- t(B_gamma)%*% B_gamma 
      
      beta_s[i,,k] <- solve(B_gamma_T_B_gamma + rho* diag(J))%*%t(B_gamma)%*%rr
      knots_s[i,,k] <- knots
      comp_save[,k] <-  B_gamma%*%beta_s[i,,k]
      g_center_save[i,k] <- mean(comp_save[,k])
      comp_save[,k] <-  comp_save[,k] - mean(comp_save[,k]) # center
      
      if(plot_opt){
        if(i > 100 & (i%%100 == 0)){
          plot( index_tmp,  comp_save[,k], main = c(i, k))
          points( index_tmp,  rr, col = "red")
          points( knots, rep(0, length(knots)), col = "green", pch = 3)
        }
      }
      
      
      if(lambda0_type == "adjusted"){
        if(i %% 100 == 0){
         if(mean(acc[(i-100+1):i,k]) < 0.2){
          if(lambda_0[k] <100000){
            lambda_0[k] <- lambda_0[k]*1.1
          }
        }
        if(mean(acc[(i-100+1):i,k]) > 0.4){
             lambda_0[k] <- lambda_0[k]/1.1
        }
        #print(lambda_0)
        #print(mean(acc[(i-100+1):i,k]))
        }
      }
     }
      # -- Step 2
      ri <- y_train - apply(comp_save, 1, sum)
      sigma2_tilde <- 1/(1/(sigma_0^2) + n_train/sigma2_s[i-1])
      mu_tilde <- sigma2_tilde*(mu_0/(sigma_0^2) + sum(ri)/sigma2_s[i-1])
      mu_s[i] <- rnorm(1,  mu_tilde, sd =   sqrt(sigma2_tilde))
      
      # --- Step 3 
      alpha_sigma2_post <- n_train/2 + alpha_1
      beta_sigma2_post <-  sum((ri - mu_s[i])^2)/2 + beta_1
      sigma2_s[i] <- 1/rgamma(1, shape = alpha_sigma2_post,  beta_sigma2_post) 
    }

  #obj_to_return <- list( acc = acc, g_center_save =  g_center_save, mu_s = mu_s, beta_s = beta_s, gamma_s = gamma_s,  knots_s = knots_s, sigma2_s =  sigma2_s, comp_save =  comp_save, 
  #                      nwarmup = nwarmup, a_s= a_s, nsample = nsample)
  
  obj_to_return <- list( acc = acc, g_center_save =  g_center_save, mu_s = mu_s, beta_s = beta_s, gamma_s = gamma_s,  knots_s = knots_s, sigma2_s =  sigma2_s,
                         nwarmup = nwarmup,  nsample = nsample)
  
  if(prior_type == "laplace_hyper"){
    obj_to_return <- c(obj_to_return, list(b_s = b_s))
  }
  
  if(prior_type == "spike_slab"){
    obj_to_return <- c(obj_to_return, list( mix_a_s =  mix_a_s,  m_s =  m_s))
  }
  if(prior_type == "laplace_hyper"){
    obj_to_return <- c(obj_to_return, list(b_s = b_s))
  }
  
  if(pred == TRUE){
    pred_train <- bppr.predict(obj_to_return, x_train, y_new = y_train)
    pred_test <- bppr.predict(obj_to_return, x_test, y_new = y_test)
    obj_to_return <- c(obj_to_return, list(pred_train = pred_train, pred_test = pred_test))
  }
  return( obj_to_return)
}

# sparse version of bppr, use spherical coordiantes 
# K is the number of additive components 
# A sparse version of bppr
# remove spike-and-slab, and do more inner iterations 

sbppr.new <- function(x_train, y_train, proposal_type = "von-Mises", prior_type = "uniform", knot_type = "quantile", x_test = NULL, y_test = NULL, K, nknots, rho = 0, b = 0.1, prior_params, init_params, nwarmup, nsample, nstep= 5,  pred = FALSE){
  
  # parameters for the prior distribution 
  mu_0 <- prior_params$mu_0;
  sigma_0 <- prior_params$sigma_0;
  alpha_1 <- prior_params$alpha_1; # for sigma
  beta_1 <- prior_params$beta_1 # for sigma
  alpha_2 <- prior_params$alpha_2; # for "b" (LaplaceH or Laplace) or for "a" (Spike-and-slab)
  beta_2 <- prior_params$beta_2; # for "b" (LaplaceH or Laplace) or for "a" (Spike-and-slab)
  
  
  n_train <- dim(x_train)[1]
  n_test <- dim(x_test)[1]
  p <-  dim(x_train)[2]
  lambda_0 <- 10000
  J <- nknots + 4  # number of basis functions 
  acc <- matrix(NA, nsample, K)
  acc[1,] <- rep(1,K)
  
  # save samples and initialize
  mu_s <- rep(NA, nsample); 
  mu_s[1] <- init_params$mu; 
  sigma2_s <- rep(NA, nsample);  
  sigma2_s[1] <- init_params$sigma2;
  
  beta_s <- array(NA, dim = c(nsample, J, K)) # intercept = true
  a_s <- array(NA, dim = c(nsample, p-1, K)) # save the angles 
  
  for(k in 1:K){
    beta_s[1,,k] <- rep(0, J)
  }
  
  gamma_s <-  array(NA, dim = c(nsample, p, K))
  b_s <-  array(NA, dim = c(nsample,K))
  knots_s <-  array(NA, dim = c(nsample, nknots + 2, K))
  
  
  if(!is.null(init_params$gammas)){
    print("init gamma provided!")
    gamma_s[1,,] <- init_params$gammas
  }else{
    for(k in 1:K){
      mu_gamma_0 <- rnorm(p)
      mu_gamma_0 <- mu_gamma_0/sqrt(sum(mu_gamma_0^2))  # some randomly initialized directions
      gamma_s[1,,k] <- rvmf(1, mu = mu_gamma_0, k = lambda_0)
    }
  }
  
  for(k in 1:K){
    a_prop <- c_to_p(gamma_s[1,,k])
    a_s[1,,k] <-  a_prop
  }
  
  # save the estimated components at the current iteration
  comp_save <- matrix(0, nrow = n_train, ncol = K)
  g_center_save <- matrix(0, nrow = nsample, ncol = K) # the constant to center each g
  
  for(i in 2:nsample){
    
    if(i%%100 == 0){
      print(paste(i, "samples!"))
    }
    
    for(k in 1:K){  
      if(prior_type == "laplace_hyper"){
        if(i == 2){
          b_s[i,k] <- 1/rgamma(1,  alpha_2 + p - 1, (beta_2 + sum(abs(a_s[i-1,,k]))))
        }else{
          if(mh_tmp$acc){ # from the previous iterations
            b_s[i,k] <- 1/rgamma(1,  alpha_2 + p - 1, (beta_2 + sum(abs(a_s[i-1,,k]))))
          }else{
            b_s[i,k] <-  b_s[i-1,k] 
          }
        }
      }
      
      b_pre  <- b_s[i,k] # previous to sampling gamma
      gamma_pre <- gamma_s[i-1,,k]
      a_pre <- a_s[i-1,,k]
      
      # compute the residuals 
      rr <- switch(as.character(K), "1" = y_train  -  mu_s[i-1],
                   "2" = y_train - comp_save[, -k] -  mu_s[i-1], 
                   y_train - apply(comp_save[, -k], 1, sum) -  mu_s[i-1])
      
      for(kk in 1:nstep){
        
        mh_tmp <- gamma.mh.s(x_train = x_train, proposal_type = proposal_type, prior_type = prior_type, knot_type = knot_type, gamma_pre = gamma_pre, 
                             b_pre = b_pre, mix_a_pre = NULL, m_pre = NULL, nknots = nknots, 
                             alpha_1 = alpha_1, beta_1 = beta_1, 
                             alpha_2 = alpha_2, beta_2 = beta_2, 
                             b_0 = NULL, b_1 = NULL, rho = rho,  
                             b = b, lambda_0 = lambda_0, rr = rr)
        if(mh_tmp$acc){ # found an accepted sample
          gamma_pre  <-  mh_tmp$gamma
          a_pre <-  mh_tmp$a
          b_pre <- 1/rgamma(1,  alpha_2 + p - 1, (beta_2 + sum(abs(a_pre))))
        }
      }
      
      acc[i,k] <-  mh_tmp$acc
      gamma_s[i,,k] <-gamma_pre
      a_s[i,,k] <-  a_pre 
      
      index_tmp <- apply(x_train, 1, function(u){t(gamma_s[i,,k])%*%u%*%gamma_s[i,,k]}) # compute the index 
      
      
      if(knot_type== "quantile_adj"){
        knots_tmp <- quantile(index_tmp, seq(0.05,0.95,length.out = nknots)) # include the boundary
        knots_delta <-0.05*(max(index_tmp)  -  min(index_tmp)) 
        knots <-  c(min(index_tmp)-  knots_delta, knots_tmp, max(index_tmp) +  knots_delta)
      }
      
      if(knot_type== "quantile"){
        knots_tmp <- quantile(index_tmp, seq(0,1,length.out = nknots+2)) # include the boundary
        knots <-    knots_tmp
      }
      
      if(knot_type == "even"){
        knots_tmp <- seq(min(index_tmp), max(index_tmp),length.out = nknots + 2) # include the boundary
        knots_delta <-0.05*(knots_tmp[nknots]  -  knots_tmp[1]) 
        knots <-  c(knots_tmp[1] -  knots_delta, knots_tmp[2:(nknots+1)], knots_tmp[nknots+2] +  knots_delta)
      }
      
      B_gamma <- ns(index_tmp,  Boundary.knots =c(knots[1], knots[nknots +2]), knots = knots[2:(2+nknots-1)],  intercept = TRUE)
      
      #B_gamma <- bSpline(index_tmp, Boundary.knots = c(knots[1],  knots[nknots + 2]), knots = knots[2:(2+nknots-1)], degree = 3, intercept = TRUE)
      B_gamma_T_B_gamma <- t(B_gamma)%*% B_gamma 
      
      beta_s[i,,k] <- solve(B_gamma_T_B_gamma + rho* diag(J))%*%t(B_gamma)%*%rr
      knots_s[i,,k] <- knots
      comp_save[,k] <-  B_gamma%*%beta_s[i,,k]
      g_center_save[i,k] <- mean(comp_save[,k])
      comp_save[,k] <-  comp_save[,k] - mean(comp_save[,k]) # center
    }
    
    # -- Step 2
    ri <- y_train - apply(comp_save, 1, sum)
    sigma2_tilde <- 1/(1/(sigma_0^2) + n_train/sigma2_s[i-1])
    mu_tilde <- sigma2_tilde*(mu_0/(sigma_0^2) + sum(ri)/sigma2_s[i-1])
    mu_s[i] <- rnorm(1,  mu_tilde, sd =   sqrt(sigma2_tilde))
    
    # --- Step 3 
    alpha_sigma2_post <- n_train/2 + alpha_1
    beta_sigma2_post <-  sum((ri - mu_s[i])^2)/2 + beta_1
    sigma2_s[i] <- 1/rgamma(1, shape = alpha_sigma2_post,  beta_sigma2_post) 
  }
  
  obj_to_return <- list( acc = acc, g_center_save =  g_center_save, mu_s = mu_s, beta_s = beta_s, gamma_s = gamma_s,  knots_s = knots_s, sigma2_s =  sigma2_s,
                         nwarmup = nwarmup,  nsample = nsample)
  
  if(prior_type == "laplace_hyper"){
    obj_to_return <- c(obj_to_return, list(b_s = b_s))
  }
  
  
  if(pred == TRUE){
    pred_train <- bppr.predict(obj_to_return, x_train, y_new = y_train)
    pred_test <- bppr.predict(obj_to_return, x_test, y_new = y_test)
    obj_to_return <- c(obj_to_return, list(pred_train = pred_train, pred_test = pred_test))
  }
  return( obj_to_return)
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



#gamma.mh.s(x_train, proposal_type, prior_type, gamma_pre,  b_pre, mix_a_pre, m_pre, nknots,  alpha_1, beta_1, alpha_2, beta_2, b0, b1, rho,  b, lambda_0, rr) 
gamma.mh.s <- function(x_train, proposal_type, prior_type,knot_type, gamma_pre, b_pre, mix_a_pre, m_pre,  nknots, alpha_1, beta_1, alpha_2, beta_2, b_0, b_1, rho, b, lambda_0, rr){
  

  n_train <- dim(x_train)[1]
  
  p <- dim(x_train)[2]
  #J <- nknots + 4
  J <- nknots + 2
  a_pre <- c_to_p(gamma_pre) # angles from previous iteration 
  

  cal_log_prop <- function(x_train, prior_type, a_tmp,  nknots, rr){
    
    gamma_tmp <- p_to_c(a_tmp)
    
    index_tmp <- apply(x_train, 1, function(u){t(gamma_tmp)%*%u%*%gamma_tmp}) # compute the index 
    
    if(knot_type== "quantile_adj"){
      knots_tmp <- quantile(index_tmp, seq(0.05,0.95,length.out = nknots)) # include the boundary
      knots_delta <-0.05*(max(index_tmp)  -  min(index_tmp)) 
      knots <-  c(min(index_tmp)-  knots_delta, knots_tmp, max(index_tmp) +  knots_delta)
    }
    
    if(knot_type== "quantile"){
      knots_tmp <- quantile(index_tmp, seq(0,1,length.out = nknots+2)) # include the boundary
      knots <-    knots_tmp
    }
    
    if(knot_type == "even"){
      knots_tmp <- seq(min(index_tmp), max(index_tmp),length.out = nknots + 2) # include the boundary
      knots_delta <-0.05*(knots_tmp[nknots]  -  knots_tmp[1]) 
      knots <-  c(knots_tmp[1] -  knots_delta, knots_tmp[2:(nknots+1)], knots_tmp[nknots+2] +  knots_delta)
    }
    
    B_gamma <- ns(index_tmp,  Boundary.knots =c(knots[1], knots[nknots +2]), knots = knots[2:(2+nknots-1)],  intercept = TRUE)
    
    #B_gamma <- bSpline(index_tmp, Boundary.knots = c(knots[1], knots[nknots +2]), knots = knots[2:(2+nknots-1)], degree = 3, intercept = TRUE)
    
    B_gamma_T_B_gamma <- t(B_gamma)%*% B_gamma 
    Sigma_rho <- solve(B_gamma_T_B_gamma + rho* diag(J))
    Sigma_0 <-  solve(B_gamma_T_B_gamma)
    inv_Sigma_0 <- B_gamma_T_B_gamma
    
    S_theta <- t(rr)%*%rr - t(rr)%*%B_gamma %*% (0.5 * Sigma_0 + Sigma_rho  - (Sigma_rho%*%inv_Sigma_0%*% Sigma_rho)/2) %*% t(B_gamma)%*%rr
    
    if(prior_type == "uniform"){
      log_D_theta <-  -(alpha_1 + n_train/2) * log(S_theta + 2*beta_1)  - sum((p-1-(1:(p-2)))*log(abs(cos(a_tmp[1:(p-2)]))))
    }

    if(prior_type == "laplace"){
      log_D_theta <-  -(alpha_1 + n_train/2) * log(S_theta + 2*beta_1) + sum(- abs(a_tmp)/b) -  sum((p-1-(1:(p-2)))*log(abs(cos(a_tmp[1:(p-2)]))))  
    }
    
    if(prior_type == "laplace_hyper"){
        log_D_theta <-  -(alpha_1 + n_train/2) * log(S_theta + 2*beta_1) + sum(-abs(a_tmp)/b) -  sum((p-1-(1:(p-2)))*log(abs(cos(a_tmp[1:(p-2)]))))  
    }
    if(prior_type == "spike_slab"){
      
      tmp_ss <- m_pre * (-abs(a_tmp)/b_0) + (1 - m_pre)*(-abs(a_tmp)/b_1)
      log_D_theta <-  -(alpha_1 + n_train/2) * log(S_theta + 2*beta_1) + sum(tmp_ss) -  sum((p-1-(1:(p-2)))*log(abs(cos(a_tmp[1:(p-2)]))))  
        
    }
    
    return(log_D_theta)
  }
  
  acc  <- 0  
  
  if(proposal_type == "von-Mises"){
    

    gamma_prop <-  generate.prop.von(gamma_pre, lambda_0)
    a_prop <- c_to_p(gamma_prop)
    acc_p <- min(1, exp(cal_log_prop(x_train, prior_type, a_prop,  nknots,  rr) -
                          cal_log_prop(x_train, prior_type, a_pre,  nknots, rr)))
    prob_tmp <- runif(1) 
    if( prob_tmp <= acc_p){
      acc <- 1
    }
  }
  
  
  if(acc == 1){
    return(list(acc = 1, gamma = gamma_prop, a = a_prop, a_prop= a_prop))
  }else{
    return(list(acc = 0))
  }
}

# make predictions for the test data 
bppr.predict <- function(bppr_obj, x_new, y_new = NULL){
  
  nsamples <- length(bppr_obj$mu_s)
  nwarmup <- bppr_obj$nwarmup
  n_new <- dim(x_new)[1]
  K <- dim(bppr_obj$gamma_s)[3]
  g_pred_s <- array(NA, dim = c(nsamples -  nwarmup, K, n_new))
  mu_pred_s <- rep(NA, nsamples -  nwarmup)
  nknots_all <- dim(bppr_obj$knots_s)[2]
  
  for(i in (nwarmup+1):nsamples){
    
    mu_pred_s[i- nwarmup] <-  bppr_obj$mu_s[i]
    
    for(k in 1:K){
      
      gamma_pre <- bppr_obj$gamma_s[i,,k]
      index_tmp <-  apply(x_new, 1, function(u){t(gamma_pre)%*%u%*%gamma_pre})
      B_gamma <- ns(index_tmp,  Boundary.knots =c(bppr_obj$knots_s[i,1,k], bppr_obj$knots_s[i,nknots_all,k]), knots =  bppr_obj$knots_s[i,,k][2:(nknots_all-1)],  intercept = TRUE)
      #B_gamma <- bSpline(index_tmp, Boundary.knots = c(bppr_obj$knots_s[i,1,k], bppr_obj$knots_s[i,nknots_all,k]), knots = bppr_obj$knots_s[i,,k][2:(nknots_all-1)], degree = 3, intercept = TRUE)
      g_pred_s[i - nwarmup,k, ]  <-   B_gamma%*% bppr_obj$beta_s[i,,k ] -  bppr_obj$g_center_save[i,k]
    }
  }
  
  if(!is.null(y_new)){
    mspe <- mean((y_new - apply(apply(g_pred_s, c(1,3), sum), 2, mean) - mean( mu_pred_s))^2)
  }
  
  return( list(mu_pred_s = mu_pred_s,  g_pred_s =  g_pred_s, mspe  =   mspe ))
}
