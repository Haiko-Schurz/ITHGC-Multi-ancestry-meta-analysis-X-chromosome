fperm = function(outcome,sex,genotype,covariates,N,c) {
## Permutation of the case-control status within females and males respectively
data_hold=cbind(outcome,sex,genotype,covariates)
m=ncol(data_hold)
data_hold1<-data_hold[data_hold$sex==1,]
nn1<-nrow(data_hold1)
outcome_hold1<-data_hold1$outcome
data_hold0<-data_hold[data_hold$sex==2,]
nn0<-nrow(data_hold0)
outcome_hold0<-data_hold0$outcome
results2<-matrix(0,N,7)
for ( i in 1:N){
#rand('seed',sum(clock)*1001);
outcome1<-sample(outcome_hold1,nn1,replace = FALSE)
data1<-data_hold1
data1$outcome<-outcome1
outcome0=sample(outcome_hold0,nn0,replace=FALSE)
data0=data_hold0
data0$outcome=outcome0
data=rbind(data1,data0)
outcome2<-data$outcome
sex2<-data$sex
genotype<-data[,3:4]
covariates<-data[,5:m]
age.covar2 <-covariates[,1]
san2 <-covariates[,4]
afr2<-covariates[,3]
eur2<-covariates[,5]
asi2<-covariates[,2]
LO.RD <- glm(outcome2 ~ sex2 + age.covar2 + san2 + afr2 + eur2 + asi2 ,family = binomial(),method = "glm.fit")
L0 <- LO.RD$deviance
#send permutated data to the fllr function to calculate likelihood ratios
source( 'fLLR.R' )
results2[i,]<-fLLR(outcome,sex,genotype,covariates,L0,c)
}

return (results2[,c(1,2,3)])
}
