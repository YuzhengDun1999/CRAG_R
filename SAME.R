library("truncnorm")
library("MASS")
library("actuar")

sample_b = function(sigma_b, sum_anno, annot, z, number){
  #Sigma is a vector of the variance of coefficient of annotation
  #sum_anno is sum A_i^T * A_i
  
  Az = 0
  for (i in 1:nrow(annot)) {
    Az = Az + matrix(annot[i,]*z[i], ncol = 1)
  }
  Az = matrix(Az, ncol = 1)
  
  for (j in 1:(ncol(annot)-1)) {
    sum_anno[j,j] = sum_anno[j,j] + 1/sigma_b
  }
  # inver_sigma_b = diag(1/sigma_b, nrow = ncol(annot))
  # sum_anno = sum_anno + inver_sigma_b
  
  Var = solve(sum_anno)
  E = Var %*% Az
  if(number == 1){
    return(matrix(mvrnorm(n = number, mu = E, Sigma = Var),nrow = number))
  }
  else{
    return(mvrnorm(n = number, mu = E, Sigma = Var))
  }
  #return a matrix, nrow=number, ncol=the number of annotation
}

sample_beta_old = function(G, geno, gamma, sigma_r, sigma_c, sigma_e, number){
  N = length(G)#total number of samples
  R = length(gamma)
  C = length(geno[1,])-R
  zero = 10^(-10)
  Sigma_x = matrix(nrow = R+C, ncol = R+C)
  for (i in 1:R) {
    Sigma_x[i,i] = 1/(sigma_r^gamma[i] * zero^(1-gamma[i]))
  }
  for (j in 1:C) {
    Sigma_x[R+j,R+j] = 1/sigma_c
  }
  Sigma_x[is.na(Sigma_x)] = 0##inverse of sigma_1
  
  Var = solve(as.numeric(1/sigma_e) * t(geno) %*% geno +  Sigma_x)
  E = Var %*% t(geno) %*% G * as.numeric(1/sigma_e)
  
  if(number == 1){
    return(matrix(mvrnorm(n = number, mu = E, Sigma = Var),nrow = number))
  }
  else{
    return(mvrnorm(n = number, mu = E, Sigma = Var))
  }
  #return a matrix, nrow=number, ncol = R+C
}


sample_beta_chol = function(G, geno, gamma, sigma_r, sigma_c, sigma_e, number){
  N = length(G)#total number of samples
  R = length(gamma)
  C = length(geno[1,])-R
  zero = 10^(-10)
  # temp = 1/(sigma_r^gamma * zero^(1-gamma))
  # temp = c(temp,rep(1/sigma_c,C))
  # Sigma_x = diag(temp) ##inverse of sigma_1
  # 
  # Var = solve(as.numeric(1/sigma_e) * t(geno) %*% geno +  Sigma_x)
  
  ###Sherman¨CMorrison¨CWoodbury formula
  temp = sigma_r^gamma * zero^(1-gamma)
  temp = c(temp,rep(sigma_c,C))
  Sigma_x = diag(temp) ##sigma_1
  X = geno/sqrt(as.numeric(sigma_e))
  temp1 = Sigma_x%*%t(X)

  Var = Sigma_x-temp1%*%solve(X%*%Sigma_x%*%t(X)+diag(nrow=N))%*%t(temp1)
  
  E = Var %*% t(geno) %*% G * as.numeric(1/sigma_e)
  CH = chol(Var) ##cholesky decomposition of Var
  E_norm = t(E) %*% solve(CH) ##expectation of dependent multivariate normal
  
  new_beta = matrix(nrow = number, ncol = R+C)
  for (i in 1:number) {
    new_beta[i,] = (E_norm + matrix(rnorm(R+C,mean = 0,sd = 1), nrow = 1)) %*% CH
  }
  return(new_beta)
  #return a matrix, nrow=number, ncol = R+C
}


sample_gamma = function(index, G, rare, b, annotation, beta, gamma, sigma_r, number){
  N = length(G)#total number of samples
  R = length(gamma)
  
  zero = 10^(-10)
  post = c()
  for (m in 1:number) {
    Probit = pnorm(annotation[index,] %*% b[m,])
    ratio = dnorm(beta[m, index], 0, sqrt(sigma_r)) / (zero + dnorm(beta[m, index], 0, sqrt(zero)))
    odds = (Probit+zero) / (1-Probit+zero) * ratio
    if(odds == Inf){
      post[m] = 1
    }
    else{
      post[m] = odds/(1+odds)   
    }
  }
  
  gamma[index] = rbinom(1,1,sum(post)/number)
  return(gamma)
}


sample_z= function(annotation, b, gamma, number){
  R = length(gamma)
  z = c()
  sum_b = colSums(b)
  for (i in 1:R) {
    mean = 1/number * annotation[i,] %*% sum_b
    Var = 1/number
    if(gamma[i] == 1){
      z_i = rtruncnorm(n = 1, a = 0, b = Inf, mean = mean, sd = sqrt(Var))
    }
    if(gamma[i] == 0){
      z_i = rtruncnorm(n = 1, a = -Inf, b = 0, mean = mean, sd = sqrt(Var))
    }
  z = c(z, z_i)
  }
  return(z)#return a vector
}


sample_e= function(G, geno, beta, alpha_e, tau_e, number){
  N = length(G)
  shape = number * alpha_e + N*number/2 + number -1
  scale = number * tau_e
  for (m in 1:number) {
    scale = scale + sum((G - geno %*% beta[m,])^2)/2
  }
  return(rinvgamma(1, shape = shape, scale = scale))
}


sample_r= function(gamma, beta, alpha_r, tau_r, number){
  R = length(gamma)
  shape = number * alpha_r + sum(gamma) * number / 2 + number -1
  scale = number * tau_r
  for (i in 1:R) {
    for (m in 1:number) {
      scale = scale + gamma[i] * beta[m,i]^2 / 2
    }
  }
  return(rinvgamma(1, shape = shape, scale = scale))
}


sample_c= function(beta, R, alpha_c, tau_c, number){
  #input of beta should be the effect size of common variants
  C = ncol(beta) - R
  shape = number * alpha_c + C*number/2 + number -1
  scale = number * tau_c + sum(beta[,(R+1):(R+C)]^2)/2
  return(rinvgamma(1, shape = shape, scale = scale))
}


sample_sigma_b= function(b, alpha_b, tau_b, number){
  K = ncol(b)-1
  shape = number * alpha_b + K * number/2 + number -1
  scale = number * tau_b + sum(b[,1:K]^2)/2
  return(rinvgamma(1, shape = shape, scale = scale))
}


SAME = function(G, rare, common, annotation, gamma, b, alpha_e, tau_e, alpha_r, tau_r, 
                alpha_c, tau_c, alpha_b, tau_b, Ite, rate, generate){
  geno = cbind(rare, common)
  R = length(gamma)
  C = ncol(common)
  #calculate sum A_i^T * A_i
  sum_anno = matrix(nrow = ncol(annotation), ncol = ncol(annotation))
  sum_anno[is.na(sum_anno)] = 0
  for(i in 1:nrow(annotation)){
    sum_anno = sum_anno + annotation[i,]%o%annotation[i,]
  }
  ### initialization ###
  z = sample_z(annotation, b, gamma, 1)
  # sigma_e = rinvgamma(1, shape = alpha_e, scale = tau_e)
  # sigma_r = rinvgamma(1, shape = alpha_r, scale = tau_r)
  # sigma_c = rinvgamma(1, shape = alpha_c, scale = tau_c)
  # sigma_b = rinvgamma(1, shape = alpha_b, scale = tau_b)
  sigma_e = 1
  sigma_r = 0.01
  sigma_c = 0.001
  sigma_b = 1
  beta = matrix(nrow = 1, ncol = R+C)
  for (i in 1:R) {
    if(gamma[i] == 0){
      beta[1,i] = 0
    }else{
      beta[1,i] = rnorm(1,0,sqrt(sigma_r))
    }
  }
  for (j in 1:C) {
    beta[1,R+j] = rnorm(1,0,sqrt(sigma_c))
  }
  
  ###test initilization
  # gamma = generate$gamma
  # b = matrix(generate$b, nrow = 1)
  # z = sample_z(annotation, b, gamma, 1)
  # sigma_e = generate$sigma_e
  # sigma_r = generate$sigma_r
  # sigma_c = generate$sigma_c
  # sigma_b = generate$sigma_b
  # beta = matrix(generate$beta, nrow = 1)
  ###store the result, to break the for loop
  Sigma_e = c()
  Sigma_r = c()
  Sigma_c = c()
  Sigma_b = c()
  Gamma = matrix(nrow = Ite, ncol = length(gamma))
  
  ### run the algorithm ###
  
  #number=1,2,2,4,4,4,4,8,8,8,8
  for (ite in 1:Ite) {
    if(trunc(log(ite, rate)) - log(ite, rate) ==0){
      number = ite
    }else{
      number = rate^(floor(log(ite, rate)))
    }
    # number = 1

    b = sample_b(sigma_b, sum_anno, annotation, z, number)
    
    ###beta needs to be initialize
    # beta = sample_beta(G, geno, gamma, beta, sigma_r, sigma_c, sigma_e, number, rate)
    # beta = sample_beta_old(G, geno, gamma, sigma_r, sigma_c, sigma_e, number)
    beta = sample_beta_chol(G, geno, gamma, sigma_r, sigma_c, sigma_e, number)
    
    for (i in 1:R) {
      gamma = sample_gamma(i, G, rare, b, annotation, beta, gamma, sigma_r, number)
    }
    gamma_test = apply(c(1:R), 1, sample_gamma, G = G, rare = rare, b = b, annotation = annotation, beta = beta, gamma = gamma, sigma_r = sigma_r, number = number)
    
    z = sample_z(annotation, b, gamma, number)
    sigma_e = sample_e(G, geno, beta, alpha_e, tau_e, number)
    sigma_r = sample_r(gamma, beta, alpha_r, tau_r, number)
    sigma_c = sample_c(beta, R, alpha_c, tau_c, number)
    sigma_b = sample_sigma_b(b, alpha_b, tau_b, number)
    
    # sigma_e = 1
    
    Sigma_e = c(Sigma_e, sigma_e)
    Sigma_r = c(Sigma_r, sigma_r)
    Sigma_c = c(Sigma_c, sigma_c)
    Sigma_b = c(Sigma_b, sigma_b)
    Gamma[ite,] = gamma
    residual = c()
    if(ite >= 2){
      residual[1] = abs(Sigma_e[ite] - Sigma_e[ite-1])/Sigma_e[ite-1]
      residual[2] = abs(Sigma_r[ite] - Sigma_r[ite-1])/Sigma_r[ite-1]
      residual[3] = abs(Sigma_c[ite] - Sigma_c[ite-1])/Sigma_c[ite-1]
      residual[4] = abs(Sigma_b[ite] - Sigma_b[ite-1])/Sigma_b[ite-1]
      residual[5] = sum(abs(Gamma[ite] - Gamma[ite-1]))/ncol(Gamma)
    
      if(max(residual)<0.01)
        break
    }
    
  }
  
  return(list(beta = beta, b = b, gamma = gamma, sigma_e = sigma_e, 
              sigma_r = sigma_r, sigma_c = sigma_c, sigma_b = sigma_b, 
              ite = ite, residual = residual, Sigma_e = Sigma_e, 
              Sigma_r = Sigma_r, Sigma_c = Sigma_c, Sigma_b = Sigma_b))
}