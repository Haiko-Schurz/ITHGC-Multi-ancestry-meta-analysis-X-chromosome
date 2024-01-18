#setwd("")
#function Output<-fxchrom(N,datafile,covfile)

# Fxchrom: Function for performing association test for X-chromosome
# genetic variants using the approach proposed in the paper: 
# Wang J and Shete S, 2014, Genetic Epidemiology
# 
# usage: OUTPUT<-FXCHROM(N,datafile,covfile)
# 
# This function will call function fLLR.m and fperm.m
#
# arguments:
#  N - number of permutations for evaluation of the p value
#
#  datafile - input data file. 
#       The data file has four columns: case-control, sex and genotypes (2 columns)
#
#  covfile - input covariate files. This file is optional. 
#
# Author: Jian Wang
# E-mail: jianwang@mdanderson.org
# Release: 1.0
# Release date: 03/17/2014

#Read in .map and .ped files
map<-read.table('.map',colClasses = "character")
data<-read.table('.ped')
g1<-7
g2<-8

#Create variables to store results
mod<-as.character(c(1:nrow(map)))
beta1<-c(1:nrow(map))
se<-c(1:nrow(map))
LLR<-c(1:nrow(map))
OR<-c(1:nrow(map))
CI95_L<-c(1:nrow(map))
CI95_H<-c(1:nrow(map))
p_value<-c(1:nrow(map))
llrp<-c(1:nrow(map))
case<-c(1:nrow(map))
control<-c(1:nrow(map))

for(i in 1:(nrow(map)) ){
print(i)
  outcome<-data[,6]   #disease status
  for(j in 1:length(outcome)){
    if (outcome[j]==1){
      outcome[j]<-0
    }
    else if(outcome[j]==2){
      outcome[j]<-1
    }
  }

#Create variables for covariates     
sex<-data[,5]
genotype<-data[,c(g1,g2)]  
covariates<-read.table('Wang_xci_covar20190829.txt',header=T)
age.covar <-covariates[,5]
san <-covariates[,3]
afr<-covariates[,2]
eur<-covariates[,4]
asi<-covariates[,1]

#remove individuals with missing genotype from all files
ind.rem1<-which((genotype[,1]==0) | (genotype[,1]==-9) | is.na(genotype[,1]) )
ind.rem2<-which((genotype[,2]==0) | (genotype[,2]==-9) | is.na(genotype[,2]))
ind.rem<-c(ind.rem1,ind.rem2)
ind.rem<-unique(ind.rem)
if (length(ind.rem)!=0){
  genotype<-genotype[-ind.rem,]
  outcome<-outcome[-ind.rem]
  sex<-sex[-ind.rem]
  age.covar<-age.covar[-ind.rem]
  san<-san[-ind.rem]
  afr<-afr[-ind.rem]
  eur<-eur[-ind.rem]
  asi<-asi[-ind.rem]
  covariates<-covariates[-ind.rem,]
}
control[i]<-sum(outcome==0, na.rm =T)
case[i]<-sum(outcome==1, na.rm =T)

## Calculate the observed Likelihood Ratio
LO.RD <- glm(outcome ~ sex + age.covar + san + afr + eur + asi ,family = binomial(),method = "glm.fit")
L0 <- LO.RD$deviance

#Call function fLLR
#Determine the XCI model that gives the highest likelihood ratio with corresponding beta0, log-likelihood ratio and standard error
c=1
source( 'fLLR.R' )
results<-fLLR(outcome,sex,genotype,covariates,L0,c)
#Record max log-likelihood ratio corresponding to a specific model
LLR_star<-results[3]
if (results[4]<1){
  model<- 'Skewed_XCI_to_normal_allele'
} else if (results[4]==1){
  model<- 'Random_XCI'
} else if (results[4]==3){
  model<- 'Escape_of_XCI'
} else {
  model<- 'Skewed_XCI_to_risk_allele'
}

#The permutations section was updated to calculate FDR instead of permutation. Uncomment the code below to use permutation function
#Call function fperm
## Permutations for p values. calculate p-value by sex stratified permutation 
#N=1000000
#c=0
#source( 'fperm.R' )
#perm_results<-fperm(outcome,sex,genotype,covariates,N,c)
#record all log likelihood ratios
#LLR_perm<-perm_results[,3]
## Calculate the empirical p values by summing LLR from permutation that are bigger than LLR for the optimal model calclated by fllr function
#p<-sum(LLR_perm>LLR_star)/N
## Output the results
#Output<-struct('model',model,'beta',results(1),'se',results(2),...
 #              'LLR',results(3),'p',p);

g1<-g1+2
g2<-g2+2

#Save results to variables
#Model
model
mod[i]<-model

#Beta
results[1]
beta1[i]<-results[1]
#se
results[2]
se[i]<-results[2]
#LLR
results[3]
LLR[i]<-results[3]
#llrp
llrp[i]<-results[8]
#p-value

p_value[i]<-results[7]
#OR #make sure about this
exp(results[1])
OR[i]<-exp(results[1])
#95%CI
exp(results[5:6])
CI95_L[i]<-exp(results[5])
CI95_H[i]<-exp(results[6])

#Create and save output file
XWAS_results<-cbind(map,mod,case,control,beta1,se,LLR,OR,CI95_L,CI95_H,p_value,llrp)
colnames(XWAS_results)<-c("Chr","SNP","Gen_dist","Pos","Model","Case","Control","Beta1","SE","LLR","OR","95CI-L","95CI-H","p","llrp")
write.table(XWAS_results, file="XWAS_wang_results_combined_MEGA_impute20190829.txt", col.names = T, row.names = F, quote = F )
}

# Add FDR adjusted p-value to results 
print(paste("Nr tests:", length(XWAS_results$p)))
#source("http://bioconductor.org/biocLite.R") #TODO: Uncomment to install qvalue package
#biocLite("qvalue") #TODO: Uncomment to install qvalue package
library(qvalue)
q.vals <- sort(qvalue(XWAS_results$p, pi0.method="bootstrap")$qvalues)
XWAS_results <- XWAS_results[order(XWAS_results$p),]
XWAS_results$Q <- q.vals[1:dim(XWAS_results)[1]]
write.table(XWAS_results, file="XWAS_wang_results_combined_MEGA_impute20190829.txt", row.names=F, quote=F, sep='\t', col.names = T)