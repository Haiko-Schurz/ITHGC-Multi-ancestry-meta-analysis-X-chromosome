############################################
fLLR<- function (outcome,sex,genotype,covariates,L0,c) {
## Log-likelihood for random XCI and skewed XCI
inter<-1 # 3 models for the random and skewed XCI
#group <-c(1:nrow(outcome))
geno_snp<-rowSums(genotype) #sum the two genotype columns together accross each row
geno_snp <- cbind(geno_snp)
gammaval<-seq(0,2,by = inter)
n<-length(gammaval)
gammaval[4] <-3
LL<-matrix(0,(n+1),1)
L<-matrix(0,(n+1),1)
t<-cbind(geno_snp, sex, covariates,deparse.level = 0)  #combine columbs into tables
mm<-ncol(t)
B<-matrix(0,(n+1),(mm+1))
pw<-matrix(0,(n+1),(mm+1))
std<-matrix(0,(n+1),(mm+1))
ci<-matrix(0,n+1,2)
age.covar <-covariates[,1]
san <-covariates[,4]
afr<-covariates[,3]
eur<-covariates[,5]
asi<-covariates[,2]
#calculate log-likelihood ratio for XCI-S and random XCI
# k=1 -> gammval = 0 XCI-S toward normal allele
#k=2 -> gammaval = 1 random XCI
#k=3 -> gammaval = 2 XCI-S toward deleterious allele
for (k in 1:n){
  x<-geno_snp-2
  x[geno_snp==3]<-gammaval[k]
  LO.RD <- glm(outcome ~ sex + age.covar + san + afr + eur + asi + x ,family = binomial(),method = "glm.fit")
  L[k] <- LO.RD$deviance
  beta <- LO.RD$coefficients
  LL[k]<-L[k]/-2
  B[k,]<-beta

  if(length(coef(summary(LO.RD))[,4])==8) {
    pw[k,]<-coef(summary(LO.RD))[,4]
    std[k,]<-coef(summary(LO.RD))[,2]
  }
  if(c==1){
    ci.temp<-confint(LO.RD)
    ci[k,]<-ci.temp[mm+1,]
    }
  
}

## Log-likelihood for escape of XCI
x<-geno_snp-2
x[sex==1]<-(geno_snp[sex==1]-2)/2
LO.RD <- glm(outcome ~ sex + age.covar + san + afr + eur + asi + x ,family = binomial(),method = "glm.fit")
L[n+1] <- LO.RD$deviance
beta <- LO.RD$coefficients
LL[n+1]<-L[n+1]/-2
B[n+1,]<-beta

if(length(coef(summary(LO.RD))[,4])==8) {
  pw[n+1,]<-coef(summary(LO.RD))[,4]
  std[n+1,]<-coef(summary(LO.RD))[,2]
}

if(c==1){
  ci.temp<-confint(LO.RD)
  ci[n+1,]<-ci.temp[mm+1,]
  }


## check the max Log-likelihood ratio to determine which XCI model works best
max_LL<-max(LL)
a<-gammaval[LL==max_LL]
max_gamma<-a[1]
if (max_gamma<1){
  m<- 1 #'skewed XCI to normal allele'
  conf95<-1
} else if (max_gamma==1){
  m<- 2 #'random XCI'
  conf95<-2
} else if (max_gamma==3){
  m<- 3 #'escape of XCI'
  conf95 <-4
} else {
  m<- 4 #'skewed XCI to risk allele'
  conf95<-3
}
i<-match(max_LL,LL)
b<-B[i[1],]
se<-std[i[1],]
pwald<-pw[i[1],]
llr<-L0-(-2*max_LL)
llrp<-pchisq(llr, 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)

## Output the results for model with max likelihood ratio
#which values for b and se are returned here
# b value must correspond to covariate X and same for se
xb <- length(b)
xse <-length(se)
pval<-length(pwald)
results<-c(b[xb], se[xse], llr,m, ci[conf95,],pwald[pval],llrp)
## results for test
results.test<-c(2*max_LL, a, max_LL, b, pwald, llr, llrp, se)
return(results)
}