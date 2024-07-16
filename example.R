rm(list = ls())
source("pbr.R")
library(pracma)
library(Matrix)
library(Rfast)
library(splines)

## function to generate data 
dat.gen <- function(seed, n_train, n_test, p, d, C, gen_funs){
  
  set.seed(seed)
  
  # --- generate predictor variables 
  x_train <- array(NA, dim = c(n_train, p, p))
  x_test <- array(NA, dim = c(n_test, p, p))

  for(i in 1:n_train){
    tmp <- randortho(p) 
    x_train[i,,] <- tmp %*% diag(runif(p, -10, 10))%*%t(tmp)
  }
  
  for(i in 1:n_test){
    tmp <- randortho(p) 
    x_test[i,,] <- tmp %*% diag(runif(p, -10, 10))%*%t(tmp)
  }
  

  Ba <- randortho(p) # the basis 
  coefs <- matrix( runif(d*p), ncol = d)
  FF <- Ba %*% coefs   
  FF <- apply(FF, 2, function(u) u/sqrt(sum(u^2)))
  
  Gamma <- matrix(FF, ncol = d)
  sparse_num <- p - 4

  Gamma_tmp <- matrix(0,p, d)
  while(rankMatrix(Gamma_tmp)[1] != d){
    Gamma_tmp <- Gamma
    for(jj in 1:d){
      Gamma_tmp[sample(1:p,  sparse_num), jj] <- 0
    }
  }
  Gamma <- Gamma_tmp
  for(jj in 1:d){
    Gamma[,jj] <-  Gamma[,jj]/sqrt(sum(Gamma[,jj]^2))
  }
  
  # generate the response
  y_train <- rep(0, n_train)
  y_test <- rep(0, n_test)

  index_train_tmp <- matrix(NA, n_train, d)
  index_test_tmp <- matrix(NA, n_test, d)
  
  for(j in 1:d){
    index_train_tmp[, j] <- apply(x_train, 1, function(u){t(Gamma[,j])%*%u%*%Gamma[,j]})
    index_test_tmp[, j] <- apply(x_test, 1, function(u){t(Gamma[,j])%*%u%*%Gamma[,j]})
  }
  
  
  g_all_train <- matrix(NA, n_train, ncol = d)
  g_all_test <- matrix(NA, n_test, ncol = d)
  
  for(j in 1:d){
    g_all_train[,j] <- gen_funs[[j]](index_train_tmp[, j])
    tmp <- mean(g_all_train[,j])
    g_all_train[,j] <- g_all_train[,j]  - tmp
    g_all_test[,j] <- gen_funs[[j]](index_test_tmp[, j]) - tmp
  }
  
  g_train <-  apply(g_all_train, 1, sum)
  g_test <- apply(g_all_test, 1, sum)
  
  y_train <-g_train + rnorm(n_train, 0, sd = C)
  y_test <- g_test + rnorm(n_test, 0, sd = C)


  return(list(index_train_tmp = index_train_tmp,  g_all_train =  g_all_train, x_train = x_train,   g_train =  g_train,  g_test =   g_test, y_train = y_train, x_test = x_test, y_test = y_test, 
              Gamma = Gamma))

}




## an example running bpr 
gen_funs <- list(f1 = function(u){-u}, f2 =function(u){-u^2/4})

seed <- 1; p <- 15; d <- 2; C <- 1
n_train <- 400; n_test <- 1000
dat <- dat.gen(seed,  n_train, n_test, p, d, C, gen_funs)

x_train <- dat$x_train
x_test <- dat$x_test
y_train <- dat$y_train
y_test <- dat$y_test


K <- 2
mu_0 <- 0; sigma_0 = 3
alpha_1  <- 1; beta_1 <- 1;   
alpha_2  <- 1; beta_2 <- 1; 

init_gammas <- init_ppr(x_train, y_train, K)
init_params <- list(mu = mean(y_train), sigma2 = var(y_train), gammas = init_gammas)
prior_params <- list(mu_0 = mu_0, sigma_0 = sigma_0, alpha_1  = alpha_1, beta_1 = beta_1, alpha_2  =  alpha_2 , beta_2 = beta_2)
pred <- FALSE
init_gammas <- init_ppr(x_train, y_train, K)
nknots <- 3;  J <- nknots + 2
nsample <- 5000
nwarmup <- 3000
res <- pbr(x_train = x_train, z_train = NULL, y_train = y_train, prior_type = "SSL", 
           x_test = x_test, z_test = NULL,  y_test = y_test, 
           K = K, nknots = nknots, rho = 0, b0 = 0.1, b1 = 1,
           prior_params = prior_params, init_params = init_params, 
           nwarmup = nwarmup, nsample = nsample, pred = pred)
  
# posterior average of gamma
apply(res$gamma_s[(nwarmup+1):nsample,,], c(2,3), mean)
dat$Gamma

# posterior samples of the ridge functions 
indices_train_1 <- matrix(NA, nsample -nwarmup, n_train)
indices_train_2 <- matrix(NA, nsample -nwarmup, n_train)
g_train_1 <- matrix(NA, nsample -nwarmup, n_train)
g_train_2 <- matrix(NA, nsample -nwarmup, n_train)



for(ii in (nwarmup+1):nsample){
  indices_train_1[ii - nwarmup,] <- apply(x_train, 1, function(u){t(res$gamma_s[ii,,1])%*%u%*%res$gamma_s[ii,,1]})
  
  B_gamma <- ns(indices_train_1[ii - nwarmup,] ,  Boundary.knots =c(res$knots_s[ii,1,1], res$knots_s[ii,J,1]), knots =  res$knots_s[ii,,1][2:(J-1)],  intercept = TRUE)
  g_train_1[ii - nwarmup, ] <- B_gamma%*% res$beta_s[ii,,1] -  res$g_center_save[ii,1]
  
  indices_train_2[ii - nwarmup,] <- apply(x_train, 1, function(u){t(res$gamma_s[ii,,2])%*%u%*%res$gamma_s[ii,,2]})
  
  B_gamma <- ns(indices_train_2[ii - nwarmup,] ,  Boundary.knots =c(res$knots_s[ii,1,2], res$knots_s[ii,J,2]), knots =  res$knots_s[ii,,2][2:(J-1)],  intercept = TRUE)
  g_train_2[ii - nwarmup, ] <- B_gamma%*% res$beta_s[ii,,2] -  res$g_center_save[ii,2]
  
}

plot(apply(indices_train_1, 2, median), apply(g_train_1, 2, median))
plot(apply(indices_train_2, 2, median), apply(g_train_2, 2, median))


## with additional continuous predictors 
p_z <- 2
z_train <- matrix(runif(p_z*n_train, -1, 1), ncol = p_z) 
z_test <- matrix(runif(p_z*n_test, -1, 1), ncol = p_z) 

y_train_z <- y_train + z_train[,1] +  (z_train[,2])^2
y_test_z <- y_test + z_test[,1] +  (z_test[,2])^2


init_gammas <- init_ppr(x_train, y_train, K)
init_params <- list(mu = mean(y_train), sigma2 = var(y_train), gammas = init_gammas)

prior_params <- list(mu_0 = mu_0, sigma_0 = sigma_0, alpha_1  = alpha_1, beta_1 = beta_1, alpha_2  =  alpha_2 , beta_2 = beta_2)

res <- pbr(x_train = x_train, z_train = z_train, y_train = y_train_z, prior_type = "SSL", 
           x_test = x_test, z_test = z_test,  y_test = y_test_z, K = K, 
           nknots = nknots, rho = 0, b0 = 0.1, b1 = 1, 
           prior_params = prior_params, init_params = init_params,
           nwarmup = nwarmup, nsample = nsample, pred = pred)


# posterior average of gamma
apply(res$gamma_s[(nwarmup+1):nsample,,], c(2,3), mean)
dat$Gamma

# posterior samples of the ridge functions 
indices_train_1 <- matrix(NA, nsample -nwarmup, n_train)
indices_train_2 <- matrix(NA, nsample -nwarmup, n_train)
g_train_1 <- matrix(NA, nsample -nwarmup, n_train)
g_train_2 <- matrix(NA, nsample -nwarmup, n_train)


for(ii in (nwarmup+1):nsample){
  indices_train_1[ii - nwarmup,] <- apply(x_train, 1, function(u){t(res$gamma_s[ii,,1])%*%u%*%res$gamma_s[ii,,1]})
  
  B_gamma <- ns(indices_train_1[ii - nwarmup,] ,  Boundary.knots =c(res$knots_s[ii,1,1], res$knots_s[ii,J,1]), knots =  res$knots_s[ii,,1][2:(J-1)],  intercept = TRUE)
  g_train_1[ii - nwarmup, ] <- B_gamma%*% res$beta_s[ii,,1] -  res$g_center_save[ii,1]
  
  indices_train_2[ii - nwarmup,] <- apply(x_train, 1, function(u){t(res$gamma_s[ii,,2])%*%u%*%res$gamma_s[ii,,2]})
  
  B_gamma <- ns(indices_train_2[ii - nwarmup,] ,  Boundary.knots =c(res$knots_s[ii,1,2], res$knots_s[ii,J,2]), knots =  res$knots_s[ii,,2][2:(J-1)],  intercept = TRUE)
  g_train_2[ii - nwarmup, ] <- B_gamma%*% res$beta_s[ii,,2] -  res$g_center_save[ii,2]
}

plot(apply(indices_train_1, 2, median), apply(g_train_1, 2, median))
plot(apply(indices_train_2, 2, median), apply(g_train_2, 2, median))

# ridge function for z
g_z_train_1 <- matrix(NA, nsample -nwarmup, n_train)
g_z_train_2 <- matrix(NA, nsample -nwarmup, n_train)

for(ii in (nwarmup+1):nsample){
  B_z_gamma <- ns(z_train[,1],  Boundary.knots =c(res$knots_s_z[ii,1,1], res$knots_s_z[ii,J,1]), knots =  res$knots_s_z[ii,,1][2:(J-1)],  intercept = TRUE)
  g_z_train_1[ii - nwarmup, ] <- B_z_gamma%*% res$beta_s_z[ii,,1] -  res$g_center_save_z[ii,1]
  
  B_z_gamma <- ns(z_train[,2],  Boundary.knots =c(res$knots_s_z[ii,1,2], res$knots_s_z[ii,J,2]), knots =  res$knots_s_z[ii,,2][2:(J-1)],  intercept = TRUE)
  g_z_train_2[ii - nwarmup, ] <- B_z_gamma%*% res$beta_s_z[ii,,2] -  res$g_center_save_z[ii,2]
}

plot(z_train[,1],   apply(g_z_train_1, 2, median))
plot(z_train[,2],   apply(g_z_train_2, 2, median))

