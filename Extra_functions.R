######################################################
#This function simulates values from a GF (by Elias)
######################################################
simulate_GF <- function(coords, kappa, variance, lambda=1) {
  # Compute the number of locations
  n <- nrow(coords)
  # Compute the distance matrix  
  dist.m <- as.matrix(dist(coords))
  # Compute the Matern correlation matrix
  cor.m <- 2^(1-lambda)/gamma(lambda)*(kappa*dist.m)^lambda*
    besselK(x=dist.m*kappa, nu=lambda)
  diag(cor.m) <- 1
  # Compute the covariance matrix
  Sigma <- variance * cor.m
  # Simulate date using standard algorithm based on Cholesky fact.
  #c(chol(Sigma) %*% rnorm(n=n,mean=0,sd=1)) wrong! see Elias's comment
  c(crossprod(chol(Sigma), rnorm(n=n,mean=0,sd=1)))
} 


######################################################
#This function fit a Poisson model with INLA assuming a Normal prior for the exposure (by Vir)
######################################################
fit.inla <- function(d, b) {
  # d is the dataframe containing the expected values, the area mean and precision for exposure, the poisson data
  # b is the fixed value for the log RR
  
  #offset: log(E_i) + beta * muhat_i
  oset <- log(d$expected) + b * d$meanexp
  
  #weights for iid random effects of exposure
  ws <- b / sqrt(d$precexp)
  
  #Fit model with INLA
  form <- y_Poisson ~ 1 +  f(id2, model = "iid") +
    f(id, ws, model = "iid", 
      hyper = list(prec = list(fixed = TRUE, initial = log(1)))) 
  
  res <- inla(form, family="poisson", data = d, offset = oset,
              control.fixed = list(prec = 0.001, prec.intercept = 0.001)
              #control.inla = list(strategy = "laplace", int.strategy = "grid",
              #  fast = FALSE, h = 0.005, dz = 0.1, diff.logdens = 1, restart = 0,
              #  tolerance = 0.005/1)
  )
  
  return(res)
}

fit.inla.method1 <- function(d, b) {
  # d is the dataframe containing the expected values, the area mean and precision for exposure, the poisson data
  # b is the fixed value for the log RR
  
  #offset: log(E_i) + beta * muhat_i
  oset <- log(d$expected) + b * d$meanexp_method1
  
  #weights for iid random effects of exposure
  ws <- b / sqrt(d$precexp_method1)
  
  #Fit model with INLA
  form <- y_Poisson ~ 1 +  f(id2, model = "iid") +
    f(id, ws, model = "iid", 
      hyper = list(prec = list(fixed = TRUE, initial = log(1)))) 
  
  res <- inla(form, family="poisson", data = d, offset = oset,
              control.fixed = list(prec = 0.001, prec.intercept = 0.001)
              #control.inla = list(strategy = "laplace", int.strategy = "grid",
              #  fast = FALSE, h = 0.005, dz = 0.1, diff.logdens = 1, restart = 0,
              #  tolerance = 0.005/1)
  )
  
  return(res)
}


fit.inla.method2 <- function(d, b) {
  # d is the dataframe containing the expected values, the area mean and precision for exposure, the poisson data
  # b is the fixed value for the log RR
  
  #offset: log(E_i) + beta * muhat_i
  oset <- log(d$expected) + b * d$meanexp_method2
  
  #weights for iid random effects of exposure
  ws <- b / sqrt(d$precexp_method2)
  
  #Fit model with INLA
  form <- y_Poisson ~ 1 +  f(id2, model = "iid") +
    f(id, ws, model = "iid", 
      hyper = list(prec = list(fixed = TRUE, initial = log(1)))) 
  
  res <- inla(form, family="poisson", data = d, offset = oset,
              control.fixed = list(prec = 0.001, prec.intercept = 0.001)
              #control.inla = list(strategy = "laplace", int.strategy = "grid",
              #  fast = FALSE, h = 0.005, dz = 0.1, diff.logdens = 1, restart = 0,
              #  tolerance = 0.005/1)
  )
  
  return(res)
}
######################################################
#This function computes the relative bias, RMSE and coverage given some values generated
#from the posterior distribution of a PARAMETER and the true or pseudotrue values
#IT CONSIDERS all the iterations TOGETHER (FOR BIAS AND RMSE)
######################################################
my_bias_RMSE_coverage_OLD = function(vec_simulatedvalues, n_simulations, pseudotrue_vec){
  matrix_simulatedvalues = matrix(vec_simulatedvalues, ncol=n_simulations)
  cov = c()
  
  for(i in 1:ncol(matrix_simulatedvalues)){
    q0025 = quantile(matrix_simulatedvalues[,i],0.025)
    q0975 = quantile(matrix_simulatedvalues[,i],0.975)
    cov[i] = (q0025 < pseudotrue_vec[i]  & pseudotrue_vec[i] < q0975)
  }
  
  bias = mean (vec_simulatedvalues - pseudotrue_vec[1])
  rmse = sqrt (mean ( (vec_simulatedvalues - pseudotrue_vec[1]) ^2))
  return(list(bias=bias,
              rmse=rmse,
              cov=sum(cov)/n_simulations*100))
}
######################################################
#This function computes the mean bias, RMSE and coverage given some values generated
#from the posterior distribution of a PARAMETER and the true or pseudotrue values
#IT CONSIDERS EACH ITERATION SEPARATELY 
######################################################
my_bias_RMSE_coverage = function(vec_simulatedvalues, n_simulations, pseudotrue_vec){
  matrix_simulatedvalues = matrix(vec_simulatedvalues, ncol=n_simulations)
  bias = rmse = cov = c()
   
  for(i in 1:ncol(matrix_simulatedvalues)){
    #Relative indexes:
    #bias[i] = mean( (matrix_simulatedvalues[,i] - pseudotrue_vec[i]) /pseudotrue_vec[i] )
    #rmse[i] = sqrt(mean(((matrix_simulatedvalues[,i] - pseudotrue_vec[i])/pseudotrue_vec[i])^2))
    
    #Absolute indexes:
    bias[i] = mean( matrix_simulatedvalues[,i] - pseudotrue_vec[i] )
    rmse[i] = sqrt(mean( (matrix_simulatedvalues[,i] - pseudotrue_vec[i])^2 ))
    
    q0025 = quantile(matrix_simulatedvalues[,i],0.025)
    q0975 = quantile(matrix_simulatedvalues[,i],0.975)
    cov[i] = (q0025 < pseudotrue_vec[i]  & pseudotrue_vec[i] < q0975)
  }
  return(list(bias=mean(bias),
              rmse=mean(rmse),
              cov=sum(cov)/n_simulations*100))
}


# this function computes the bias and rmse in a vector form (length given by the number of simulations)
my_bias_RMSE_vec = function(vec_simulatedvalues, n_simulations, pseudotrue_vec){
  matrix_simulatedvalues = matrix(vec_simulatedvalues, ncol=n_simulations)
  bias = rmse = c()
  
  for(i in 1:ncol(matrix_simulatedvalues)){
    bias[i] = mean( matrix_simulatedvalues[,i] - pseudotrue_vec[i] )
    rmse[i] = sqrt(mean( (matrix_simulatedvalues[,i] - pseudotrue_vec[i])^2 ))
  }
  return(data.frame(bias=bias,
              rmse=rmse))
}

######################################################
#This function computes the relative bias
# and the coverage given some values generated
#from the exposure posterior distribution of each area
######################################################
area_bias_coverage_JOINT = function(list_simulatedvalues, n_simulations, true_exposure, n_areas){
  #list_simulatedvalues is a list with n_areas elements. Each element of the list is a vector with n_simulations*n_random_INLA values
  bias = c()
  cov = c()
  
  for(A in 1:n_areas){
    #Change the list in matrix
    matrix_simulatedvalues = matrix(list_simulatedvalues[[A]], ncol=n_simulations)
  
    temp_bias = c()
    temp_cov = c()
    #For each simulation compute the relative bias and the coverage (nb: the true values change)
    for(j in 1:n_simulations){
      #temp_bias[j] = mean((matrix_simulatedvalues[,j] - true_exposure[A]) / true_exposure[A])
      temp_bias[j] = mean(matrix_simulatedvalues[,j] - true_exposure[A])
      q0025 = quantile(matrix_simulatedvalues[,j],0.025)
      q0975 = quantile(matrix_simulatedvalues[,j],0.975)
      temp_cov[j] = q0025 < true_exposure[A]  & true_exposure[A] < q0975
    }     
    
    bias[A] = mean(temp_bias)  #mean over simulations
    cov[A] = sum(temp_cov)/n_simulations * 100
  }
  return(list(bias=bias,coverage=cov))
}

######################################################
#This function computes the relative bias
#and the coverage given some values generated
#from the exposure posterior distribution of each area
#This is used for values sampled from the joint posterior
######################################################
area_bias_coverage_PI_FF = function(list_simulatedvalues, n_simulations, true_exposure, n_areas){
  #list_simulatedvalues is a list with n_simulation elements. Each element of the list is a matrix with dimension n_areas*n_random_INLA values
  
  bias = matrix(NA,nrow=n_simulations,ncol=n_areas)
  cov = matrix(NA,nrow=n_simulations,ncol=n_areas)
  
  for(i in 1:n_simulations){
    temp_bias = c()
    temp_cov = c()
    
    matrix_simulatedvalues = list_simulatedvalues[[i]]
    
    #For each area compute the relative bias and the coverage (nb: the true values change)
    for(A in 1:n_areas){
      #temp_bias[A] = mean((matrix_simulatedvalues[A,] - true_exposure[A]) / true_exposure[A])
      temp_bias[A] = mean(matrix_simulatedvalues[A,] - true_exposure[A])
      q0025 = quantile(matrix_simulatedvalues[A,],0.025)
      q0975 = quantile(matrix_simulatedvalues[A,],0.975)
      temp_cov[A] = q0025 < true_exposure[A]  & true_exposure[A] < q0975
      }     
    
    bias[i,] = temp_bias
    cov[i,] = temp_cov
    }
  
  bias = apply(bias,2,mean)  #mean over simulations
  coverage = apply(cov,2,sum) / n_simulations * 100
  
  return(list(bias=bias,coverage=coverage))
}



######################################################
#This function computes the relative bias
#and the coverage given some values generated
#from the exposure posterior distribution of each area
#the true exposure changes for each iteration
#This is used for values sampled from the joint posterior
######################################################
area_bias_coverage_PI_FF_changingexposure = function(list_simulatedvalues, n_simulations, true_exposure_list, n_areas){
  #list_simulatedvalues is a list with n_simulation elements. Each element of the list is a matrix with dimension n_areas*n_random_INLA values
  
  bias = matrix(NA,nrow=n_simulations,ncol=n_areas)
  rmse = matrix(NA,nrow=n_simulations,ncol=n_areas)
  cov = matrix(NA,nrow=n_simulations,ncol=n_areas)
  
  
  for(i in 1:n_simulations){
    temp_bias = c()
    temp_rmse = c()
    temp_cov = c()
    true_exposure = true_exposure_list[[i]]$weighted_mean_true_exposure
    matrix_simulatedvalues = list_simulatedvalues[[i]]
    
    #For each area compute the relative bias and the coverage (nb: the true values change)
    for(A in 1:n_areas){
      #temp_bias[A] = mean((matrix_simulatedvalues[A,] - true_exposure[A]) / true_exposure[A])
      temp_bias[A] = mean(matrix_simulatedvalues[A,] - true_exposure[A])
      temp_rmse[A] = mean( (matrix_simulatedvalues[A,] - true_exposure[A])^2 )
      q0025 = quantile(matrix_simulatedvalues[A,],0.025)
      q0975 = quantile(matrix_simulatedvalues[A,],0.975)
      temp_cov[A] = q0025 < true_exposure[A]  & true_exposure[A] < q0975
    }     
    
    bias[i,] = temp_bias
    rmse[i,] = temp_rmse
    cov[i,] = temp_cov
  }
  
  bias = apply(bias,2,mean)  #mean over simulations
  rmse = sqrt(apply(rmse,2,mean))  #mean over simulations
  coverage = apply(cov,2,sum) / n_simulations * 100
  
  return(list(bias=bias,rmse=rmse,coverage=coverage))
}

