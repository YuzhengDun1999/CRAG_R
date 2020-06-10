#!/usr/bin/R
# @ annotation: annotation matrix  with the intercept in the last column
# @ h2_anno: the variance of a SNP being effective or not explained by annotation
# @ h2_x: heritability of y/gene expression
# @ pi: the percentage of SNPs being effective/causal
# @ common: the genotype matrix of common SNPs
# @ rare: the genotype matrix of rare SNPs
generation_onlyCommon <- function(common, annotation,h2_anno,h2_x,pi){

  G = c()
  N = nrow(common)
  C = ncol(common)
  # first deal with NA in the genoytype matrix
  which(is.na(common),arr.ind=T) -> inds_na
  mean_col <- colMeans(common,na.rm = T)
  common[inds_na] <- mean_col[inds_na[,2]]
  
  geno = common
  # first scale the genotype matrix
  geno <- scale(geno)
  common <- scale(common)
#  geno = cbind(rare, common)
  
  
  ###real parameters
  # number of SNPs being causal 
  p_causal <- floor(pi*C)
  sigma_c = h2_x/p_causal
  sigma_b = h2_anno/(ncol(annotation)-1)
  b = rnorm(ncol(annotation)-1, mean = 0, sd = sqrt(sigma_b))
  # first determine random noise in the annotation layer
  anno_error <- ((1 - h2_anno)/h2_anno)*var(as.matrix(annotation[,-c(ncol(annotation))],ncol=length(b)) %*% as.matrix(b,ncol=1))
  b <- b/sqrt(anno_error[1,1])
  sigma_b <- h2_anno/(ncol(annotation)-1)/anno_error[1,1]
  # determine the intercept based on how many values larger than 0
  anno_noise <- rnorm(C,mean=0,sd=1)
  anno_val <- as.matrix(annotation[,-c(ncol(annotation))]) %*% as.matrix(b,ncol=1) + anno_noise
  anno_intercept <- -anno_val[order(anno_val,decreasing = T)[p_causal+1]]
  b <- c(b,anno_intercept)
  w <- as.matrix(annotation) %*% as.matrix(b,ncol=1) + anno_noise
  
  gamma = rep(0,C)
  beta = rep(0,C)
  gamma[which(w>0)] <- 1
  beta[which(w>0)] <- rnorm(sum(w>0),mean=0,sd=sqrt(sigma_c))
  sigma_e = (1-h2_x)/h2_x* var(geno %*% beta)[1,1]
  noise <- rnorm(N,mean = 0, sd=sqrt(sigma_e))
  noise <- noise/sd(noise)*sqrt(sigma_e)

  G <- geno %*% beta + noise

  return(list(G = G, common = common, rare = rare, gamma = gamma, beta = beta, sigma_c = sigma_c, sigma_e = sigma_e, sigma_b = sigma_b, b = b))
}


# generate for both common and rare
# @annotation: annotation for rare SNPs
# @h2_anno: the variance of rare SNPs getting selected explained by annnotation
generation_CommonRare <- function(common,rare, annotation,h2_anno,h2_r,h2_c,pi){
  
  G = c()
  #  N = nrow(rare)
  R = ncol(rare)
  N = nrow(common)
  C = ncol(common)
  # first deal with NA in the genoytype matrix
  which(is.na(common),arr.ind=T) -> inds_na
  mean_col <- colMeans(common,na.rm = T)
  common[inds_na] <- mean_col[inds_na[,2]]
  
  which(is.na(rare),arr.ind=T) -> inds_na
  mean_col <- colMeans(rare,na.rm = T)
  rare[inds_na] <- mean_col[inds_na[,2]]
  
  geno = cbind(rare,common)
  # first scale the genotype matrix
  geno <- scale(geno)
  common <- scale(common)
  rare <- scale(rare)
  #  geno = cbind(rare, common)
  
  
  ###real parameters
  # number of SNPs being causal 
  p_causal <- floor(pi*R)
  sigma_r = h2_r/p_causal

  sigma_b = h2_anno/(ncol(annotation)-1)
  b = rnorm(ncol(annotation)-1, mean = 0, sd = sqrt(sigma_b))
  #b = c(b,-1)
  # first determine random noise in the annotation layer
  anno_error <- ((1 - h2_anno)/h2_anno)*var(as.matrix(annotation[,-c(ncol(annotation))],ncol=length(b)) %*% as.matrix(b,ncol=1))
  b <- b/sqrt(anno_error[1,1])
  sigma_b <- h2_anno/(ncol(annotation)-1)/anno_error[1,1]
  
  # determine the intercept based on how many values larger than 0
  anno_noise <- rnorm(R,mean=0,sd=1)
  anno_val <- as.matrix(annotation[,-c(ncol(annotation))]) %*% as.matrix(b,ncol=1) + anno_noise
  anno_intercept <- -anno_val[order(anno_val,decreasing = T)[p_causal+1]]
  b <- c(b,anno_intercept)
  w <- as.matrix(annotation) %*% as.matrix(b,ncol=1) + anno_noise
  
  gamma = rep(0,R)
  beta_rare = rep(0,R)
  gamma[which(w>0)] <- 1
  beta_rare[which(w>0)] <- rnorm(sum(w>0),mean=0,sd=sqrt(sigma_r))
  
  sigma_c <- h2_c/C
  beta_common = rnorm(C,mean=0,sd=sqrt(sigma_c))
  
  sigma_e = (1-h2_c-h2_r)/(h2_c+h2_r)* var(rare %*% beta_rare + common %*% beta_common)[1,1]
  noise <- rnorm(N,mean = 0, sd=sqrt(sigma_e))
  noise <- noise/sd(noise)*sqrt(sigma_e)
  
  beta <- c(beta_rare,beta_common)
  G <- geno %*% beta + noise

  return(list(G = G, common = common, rare = rare, gamma = gamma, beta = beta, sigma_r = sigma_r, sigma_c = sigma_c, sigma_e = sigma_e, sigma_b = sigma_b, b = b))
}


generation_old = function(rare, common, annotation){
  G = c()
  N = nrow(rare)
  R = ncol(rare)
  C = ncol(common)
  geno = cbind(rare, common)
  
  ###real parameters
  sigma_r = 0.003
  sigma_c = 0.0005
  # sigma_r = 5
  # sigma_c = 3
  sigma_b = 1
  b = rnorm(ncol(annotation)-1, mean = 0, sd = sqrt(sigma_b))
  b = c(b,-1)
  gamma = c()
  beta = c()
  w = pnorm(annotation %*% b, 0, 1)
  k = 0
  for (i in 1:R) {
    if(w[i] < 0.5){
      gamma[i] = 0
      beta[i] = 0
    }
    if(w[i] >= 0.5){
      gamma[i] = 1
      k = k + 1
      beta[i] = rnorm(1, mean = 0, sd = sqrt(sigma_r))
    }
  }
  
  
  for (j in 1:C) {
    beta[R+j] = rnorm(1, mean = 0, sd = sqrt(sigma_c))
  }
  # sigma_e = 1 - sigma_r * sum(gamma) - sigma_c * ncol(common)
  sigma_e = 1 - var(geno %*% beta)
  # sigma_e = 0.1
  for (n in 1:N) {
    G_n = rnorm(1, geno[n,] %*% beta, sqrt(sigma_e))
    G = c(G, G_n)
  }
  return(list(G = G, gamma = gamma, beta = beta, sigma_r = sigma_r, sigma_c = sigma_c, sigma_e = sigma_e, sigma_b = sigma_b, b = b))
}
