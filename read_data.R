setwd('G:/thesis/code')
library('genio')
bim = read_bim("./data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_TRIM17_chr1_227407935-229416861.bim")
fam = read_fam("./data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_TRIM17_chr1_227407935-229416861.fam")
bed = read_bed("./data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_TRIM17_chr1_227407935-229416861.bed", bim$id, fam$id)
frq = read.table("GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_TRIM17_chr1_227407935-229416861.frq",header = TRUE)
ano = read.table("annotation_Whole_Blood.txt",header = TRUE)
ano = ano[,-c(1:3)]
rare_index = c()
common_index = c()
ultrarare = c()
for(i in 1:length(frq$MAF)){
  if(frq$MAF[i]>0.005 & frq$MAF[i]<0.05){
    rare_index = c(rare_index, i)
  }
  else if(frq$MAF[i]>=0.05){
    common_index = c(common_index, i)
  }
  if(frq$MAF[i]>0.005 & frq$MAF[i]<0.01){
    ultrarare = c(ultrarare, i)
  }
}


rare_geno = t(bed[rare_index,])#rare variants genotype matrix,row=samples,col=rare variant
common_geno = t(bed[common_index,])#common variants genotype matrix,row=samples,col=common variant
# annotation = matrix(unlist(ano[common_index,1:6]), nrow = length(common_index), ncol = 6)#annotation of common variants
# annotation = cbind(annotation, rep(1,nrow(annotation)))#add intercept


# new generation function
# generate <-generation_onlyCommon(common_geno,annotation,h2_anno = 0.5,h2_x = 0.2,pi = 0.1)

# for common and rare together
annotation = matrix(unlist(ano[rare_index,1:6]), nrow = length(rare_index), ncol = 6)#annotation of rare variants
annotation = cbind(annotation, rep(1,nrow(annotation)))#add intercept

###test
rare_geno = rare_geno[,1:500]
# rare_geno = common_geno[,1:1000]
common_geno = common_geno[,1:500]
annotation = annotation[1:500,]

# new generation function
generate <- generation_CommonRare(common = common_geno,rare = rare_geno,annotation = annotation, 
                                  h2_anno = 0.5,h2_r = 0.2,h2_c = 0.1,pi = 0.05)

G1 = generate$G
rare_geno1 = generate$rare
common_geno1 = generate$common

G = generate$G[1:650]
rare_geno = generate$rare[1:650,]
common_geno = generate$common[1:650,]
rare_training_duplicate = ncol(rare_geno) - ncol(rare_geno[,!duplicated(t(rare_geno))])
common_training_duplicate = ncol(common_geno) - ncol(common_geno[,!duplicated(t(common_geno))])
rare_all_duplicate = ncol(rare_geno1) - ncol(rare_geno1[,!duplicated(t(rare_geno1))])
common_all_duplicate = ncol(common_geno1) - ncol(common_geno1[,!duplicated(t(common_geno1))])

initial_gamma = function(x, y ,length){
  model = lm(y~x)
  model = summary(model)
  p = model$coefficients[2,4]
  # p=1-pf(model$fstatistic[1],model$fstatistic[2],model$fstatistic[3])
  if(p < (0.05)){
    return(1)
  }
  else{
    return(0)
  }
}

initial_p = function(x, y ,length){
  model = lm(y~x)
  model = summary(model)
  p = model$coefficients[2,4]
  return(p)
}

### initialization
p_value = apply(rare_geno, 2, initial_p, y = G, length = ncol(rare_geno))
p_order = rank(abs(p_value))
# top_gamma = generate$gamma[p_order<100]
gamma = rep(0,ncol(rare_geno))
gamma[which(p_order<130)] = 1
# gamma = apply(rare_geno, 2, initial_gamma, y = G, length = ncol(rare_geno))
b = c()
for (k in 1:6) {
  b[k] = 0
}
b[7] = -5
b = matrix(b, nrow = 1)
alpha_e = 0.01
tau_e = 0.01 
alpha_r = 0.01
tau_r = 0.01 
alpha_c = 0.01 
tau_c = 0.01 
alpha_b = 0.01 
tau_b = 0.01
Ite = 200
rate = 2

time1 = Sys.time()
results = SAME(G, rare_geno, common_geno, annotation, gamma, b, alpha_e, tau_e, alpha_r, tau_r,
                alpha_c, tau_c, alpha_b, tau_b, Ite, rate, generate)
# results = SAME_Test(G, rare_geno, common_geno, annotation, gamma, b, alpha_e, tau_e, alpha_r, tau_r,
#                alpha_c, tau_c, alpha_b, tau_b, Ite, rate, generate)
time2 = Sys.time()
time = time2 - time1
# plot(colMeans(results$beta)-generate$beta)

### test for different initialization


### evaluation ##
##R-square
geno = cbind(rare_geno,common_geno)
G_pre = geno %*% colMeans(results$beta)##predicted gene expression on training dataset
G_pre_new = geno %*% (colMeans(results$beta)*c(results$gamma,rep(1,ncol(common_geno))))
SSR = sum((G-G_pre)^2)
SST = sum((G-mean(G))^2)
R_square = 1-SSR/SST
SSR1=sum((G-geno %*% generate$beta)^2)
true_R = 1-SSR1/SST

rare_geno_test = rare_geno1[nrow(rare_geno):nrow(rare_geno1),]##rare variant for test
common_geno_test = common_geno1[nrow(common_geno):nrow(common_geno1),]##common variant for test
geno_test = cbind(rare_geno_test,common_geno_test)##genotype matrix for test
G_pre_test = geno_test %*% colMeans(results$beta)##predicted gene expression for test
G_test = generate$G[nrow(rare_geno):nrow(rare_geno1)]##real gene expression for test
cor_generate = cor(G_test,geno_test%*%generate$beta)##correlation between real gene expression and real beta
cor_results_mean = cor(G_test,geno_test%*%colMeans(results$beta))##correlation between real gene expression and results beta
cor_true_training = cor(G,cbind(rare_geno,common_geno)%*%generate$beta)
cor_pre_training = cor(G,cbind(rare_geno,common_geno)%*%colMeans(results$beta))
##h2 on predicted dataset
h2_true = var(geno_test%*%generate$beta)/var(G_test)
h2_pre = var(geno_test%*%colMeans(results$beta))/var(G_test)
##h2 on traning dataset
h2_true_training = var(cbind(rare_geno,common_geno)%*%generate$beta)/var(G)
h2_pre_training = var(cbind(rare_geno,common_geno)%*%colMeans(results$beta))/var(G)

log_likelihood = sum(log((dnorm((G-geno %*% colMeans(results$beta))/sqrt(results$sigma_e)))))
zero = 10^(-10)
log_likelihood = log_likelihood + sum(log(dnorm(colMeans(results$beta)[1:ncol(rare_geno)]/sqrt(results$gamma*results$sigma_r+zero))))
log_likelihood = log_likelihood + sum(log(dnorm(colMeans(results$beta)[(1+ncol(rare_geno)):ncol(results$beta)]/sqrt(results$sigma_c))))
Probit = pnorm(annotation%*%colMeans(results$b))
log_likelihood = log_likelihood + sum(results$gamma*log(Probit)+(1-results$gamma)*log(1-Probit))

write.csv(list(predicted=G_pre,observed=G,R_square=R_square,true_R=true_R),file = "expression.csv")
write.csv(list(cor_true_predicted=cor_generate,cor_pre_predicted=cor_results_mean,h2_true_training=h2_true_training,h2_pre_training=h2_pre_training),file = "cor.csv")
write.csv(list(h2_true_predicted=h2_true,h2_pre_predicted=h2_pre,h2_true_training=h2_true_training,h2_pre_training=h2_pre_training),file = "h2.csv")


# p_value = apply(rare_geno, 2, initial_p, y = G, length = ncol(rare_geno))
# p_order = rank(abs(p_value))
# top_gamma = generate$gamma[p_order<100]

# cor_generate = cor(G,cbind(rare_geno,common_geno)%*%generate$beta)
# cor_results_mean = cor(G,cbind(rare_geno,common_geno)%*%colMeans(results$beta))
# cor_results = cor(G,cbind(rare_geno,common_geno)%*%results$beta[1,])
# write.csv(list(generate=cor_generate,results_mean=cor_results_mean,results=cor_results),file = "cor.csv")
# write.csv(list(sigma_r=generate$sigma_r,sigma_c=generate$sigma_c,sigma_e=generate$sigma_e,sigma_b=generate$sigma_b),file = "true.csv")
# write.csv(list(sigma_r=results$sigma_r,sigma_c=results$sigma_c,sigma_e=results$sigma_e,sigma_b=results$sigma_b),file = "result.csv")
# MSE_mean = sum((G-cbind(rare_geno,common_geno)%*%colMeans(results$beta))^2)/length(G)
# MSE_all = c()
# for (i in 1:nrow(results$beta)) {
#   MSE_all = c(MSE_all, sum((G-cbind(rare_geno,common_geno)%*%results$beta[i,])^2)/length(G))
# }
# write.csv(list(mean = MSE_mean, min = min(MSE_all)), file = "mse.csv")