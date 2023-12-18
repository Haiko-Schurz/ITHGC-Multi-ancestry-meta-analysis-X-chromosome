
####################################################################################################
#   Male only X chromosome analysis
####################################################################################################

#  Europe as reference  ####################################################################

#Read in Xwas logistic regression association testing output file.
#XWAS software is available here https://academic.oup.com/jhered/article/106/5/666/2960036
eur<-read.table("russia_fish.xstrat.logistic",stringsAsFactors=FALSE,header=TRUE)
eur<-eur[,c(1,2,3,4,7,8)]
#Required columbs
#SNP: SNP_ID
#BP: Base pair number
#CHR: Chromosome number
#P: P-value

#Run the script 4 times (for males), each time with a different p-value range (different data.p dataset)
#Iteration 1 
data.p<-eur[eur$P<=0.001,]
#Iteration 2
data.p<-eur[eur$P>0.001 & eur$P<=0.01,]
#Iteration 3
data.p<-eur[eur$P>0.01 & eur$P<=0.5,]
#Iteration 4
data.p<-eur[eur$P>0.5,]


#Find list of SNPs within p-value range for each chromosome (in this case only X) within european data
chr<-sort(unique(data.p$CHR))
snp<-c()
j<-1
for(i in 1:length(chr)){
  cdata<-data.p[data.p$CHR==chr[i],]
  cdata<-cdata[order(cdata$BP),]
  
  while(nrow(cdata)>0){
    snp[j]<-cdata$SNP[1]
    r<-which(cdata$BP<=(cdata$BP[1]+500000))
    if(length(r)!=0){
      cdata<-cdata[-r,]
    }
    j<-j+1
  }
}

##Extract selected European SNPs from African and Asian XWAS meta-analysis output file
afr<-read.table("fisher_male_SAC.male.meta",stringsAsFactors=FALSE,header=TRUE)
afr<-afr[,c(1,2,3,4,7,9)]
asi<-read.table("fisher_male_china.male.meta",stringsAsFactors=FALSE,header=TRUE)
asi<-asi[,c(1,2,3,4,7,9)]

for (i in 1:length(snp)){
  e<-which(eur$SNP==snp[i])
  af<-which(afr$SNP==snp[i])
  as<-which(asi$SNP==snp[i])
  if (i ==1){
    eur.snp<-eur[e,]
    afr.snp<-afr[af,]
    asi.snp<-asi[as,]
  } else {
    eur.snp<-rbind(eur.snp,eur[e,])
    afr.snp<-rbind(afr.snp,afr[af,])
    asi.snp<-rbind(asi.snp,asi[as,])
  }
}

# Concordance analysis European vs. Africa
eur.snp<-eur.snp[,c(2,4,5,6)]
colnames(eur.snp)<-c("rs_number","reference_allele_eur","OR_eur","p.value_eur")
afr.snp<-afr.snp[,c(3,4,6,5)]
colnames(afr.snp)<-c("rs_number","reference_allele_afr","OR_afr","p.value_afr")
eur.afr<-merge(eur.snp,afr.snp,by.x="rs_number",by.y="rs_number")
R<-length(which((eur.afr$OR_eur < 1 & eur.afr$OR_afr < 1) | (eur.afr$OR_eur > 1 & eur.afr$OR_afr > 1)))
N<-nrow(eur.afr)
binom.test(R,N,0.5,alternative="greater") 

# Concordance analysis European vs. Asia
asi.snp<-asi.snp[,c(3,4,6,5)]
colnames(asi.snp)<-c("rs_number","reference_allele_asi","OR_asi","p.value_asi")
eur.asi<-merge(eur.snp,asi.snp,by.x="rs_number",by.y="rs_number")
R<-length(which((eur.asi$OR_eur < 1 & eur.asi$OR_asi < 1) | (eur.asi$OR_eur > 1 & eur.asi$OR_asi > 1)))
N<-nrow(eur.asi)
binom.test(R,N,0.5,alternative="greater") 

####  Africa as reference  ##########################################################
#Read in XWAS meta-analysis association testing output file. 
afr<-read.table("fisher_male_africa.male.meta",stringsAsFactors=FALSE,header=TRUE)
afr<-afr[,c(1,2,3,4,7,9)]
#Required columbs
#SNP: SNP_ID
#BP: Base pair number
#CHR: Chromosome number
#P: P-value

#Run the script 4 times, each time with a different p-value range (different data.p dataset)
#Iteration 1 
data.p<-afr[afr$P<=0.001,]
#Iteration 2
data.p<-afr[afr$P>0.001 & afr$P<=0.01,]
#Iteration 3
data.p<-afr[afr$P>0.01 & afr$P<=0.5,]
#Iteration 4
data.p<-afr[afr$P>0.5,]

#Find list of SNPs within p-value range for each autosomal chromosome within African data
chr<-sort(unique(data.p$CHR))
snp<-c()
j<-1
for(i in 1:length(chr)){
  cdata<-data.p[data.p$CHR==chr[i],]
  cdata<-cdata[order(cdata$BP),]
  
  while(nrow(cdata)>0){
    snp[j]<-cdata$SNP[1]
    r<-which(cdata$BP<=(cdata$BP[1]+500000))
    if(length(r)!=0){
      cdata<-cdata[-r,]
    }
    j<-j+1
  }
}

#Extract selected African SNPs from European and Asian XWAS meta-analysis output file
eur<-read.table("russia_fish.xstrat.logistic",stringsAsFactors=FALSE,header=TRUE)
eur<-eur[,c(1,2,3,4,7,8)]
asi<-read.table("fisher_male_china.male.meta",stringsAsFactors=FALSE,header=TRUE)
asi<-asi[,c(1,2,3,4,9,7)]
for (i in 1:length(snp)){
  e<-which(eur$SNP==snp[i])
  af<-which(afr$SNP==snp[i])
  as<-which(asi$SNP==snp[i])
  if (i ==1){
    eur.snp<-eur[e,]
    afr.snp<-afr[af,]
    asi.snp<-asi[as,]
  } else {
    eur.snp<-rbind(eur.snp,eur[e,])
    afr.snp<-rbind(afr.snp,afr[af,])
    asi.snp<-rbind(asi.snp,asi[as,])
  }
}


#Concordance analysis Africa vs. Europe
eur.snp<-eur.snp[,c(2,4,5,6)]
colnames(eur.snp)<-c("rs_number","reference_allele_eur","OR_eur","p.value_eur")
afr.snp<-afr.snp[,c(3,4,6,5)]
colnames(afr.snp)<-c("rs_number","reference_allele_afr","OR_afr","p.value_afr")
eur.afr<-merge(eur.snp,afr.snp,by.x="rs_number",by.y="rs_number")
R<-length(which((eur.afr$OR_eur < 1 & eur.afr$OR_afr < 1) | (eur.afr$OR_eur > 1 & eur.afr$OR_afr > 1)))
N<-nrow(eur.afr)
binom.test(R,N,0.5,alternative="greater") 

#Concordance analysis Africa vs. Asia
asi.snp<-asi.snp[,c(3,4,5,6)]
colnames(asi.snp)<-c("rs_number","reference_allele_asi","OR_asi","p.value_asi")
afr.asi<-merge(afr.snp,asi.snp,by.x="rs_number",by.y="rs_number")
R<-length(which((afr.asi$OR_afr < 1 & afr.asi$OR_asi < 1) | (afr.asi$OR_afr > 1 & afr.asi$OR_asi > 1)))
N<-nrow(afr.asi)
binom.test(R,N,0.5,alternative="greater") 

###  Asia as reference  #############################################################
#Read in XWAS meta-analysis association testing output file. 
asi<-read.table("fisher_male_africaSAC.male.meta",stringsAsFactors=FALSE,header=TRUE)
asi<-asi[,c(1,2,3,4,7,9)]
#Required columbs
#SNP: SNP_ID
#BP: Base pair number
#CHR: Chromosome number
#P: P-value

#Run the script 4 times, each time with a different p-value range (different data.p dataset)
#Iteration 1 
data.p<-asi[asi$P<=0.001,]
#Iteration 2
data.p<-asi[asi$P>0.001 & asi$P<=0.01,]
#Iteration 3
data.p<-asi[asi$P>0.01 & asi$P<=0.5,]
#Iteration 4
data.p<-asi[asi$P>0.5,]

#Find list of SNPs within p-value range for each autosomal chromosome within Asian data
chr<-sort(unique(data.p$CHR))
snp<-c()
j<-1
for(i in 1:length(chr)){
  cdata<-data.p[data.p$CHR==chr[i],]
  cdata<-cdata[order(cdata$BP),]
  
  while(nrow(cdata)>0){
    snp[j]<-cdata$SNP[1]
    r<-which(cdata$BP<=(cdata$BP[1]+500000))
    if(length(r)!=0){
      cdata<-cdata[-r,]
    }
    j<-j+1
  }
}

#Extract selected Asian SNPs from European and African plink v2 meta-analysis output file
eur<-read.table("russia_fish.xstrat.logistic",stringsAsFactors=FALSE,header=TRUE)
eur<-eur[,c(1,2,3,4,7,8)]
afr<-read.table("fisher_male_africaSAC.male.meta",stringsAsFactors=FALSE,header=TRUE)
afr<-afr[,c(1,2,3,4,7,9)]

for (i in 1:length(snp)){
  e<-which(eur$SNP==snp[i])
  af<-which(afr$SNP==snp[i])
  as<-which(asi$SNP==snp[i])
  if (i ==1){
    eur.snp<-eur[e,]
    afr.snp<-afr[af,]
    asi.snp<-asi[as,]
  } else {
    eur.snp<-rbind(eur.snp,eur[e,])
    afr.snp<-rbind(afr.snp,afr[af,])
    asi.snp<-rbind(asi.snp,asi[as,])
  }
}

#Concordance analysis Asia vs. Europe
eur.snp<-eur.snp[,c(2,4,5,6)]
colnames(eur.snp)<-c("rs_number","reference_allele_eur","OR_eur","p.value_eur")
asi.snp<-asi.snp[,c(3,4,6,5)]
colnames(asi.snp)<-c("rs_number","reference_allele_asi","OR_asi","p.value_asi")
eur.asi<-merge(eur.snp,asi.snp,by.x="rs_number",by.y="rs_number")
R<-length(which((eur.asi$OR_eur < 1 & eur.asi$OR_asi < 1) | (eur.asi$OR_asi > 1 & eur.asi$OR_eur > 1)))
N<-nrow(eur.asi)

binom.test(R,N,0.5,alternative="greater") 

#Concordance analysis Asia vs. Africa
afr.snp<-afr.snp[,c(3,4,6,5)]
colnames(afr.snp)<-c("rs_number","reference_allele_afr","OR_afr","p.value_afr")
afr.asi<-merge(afr.snp,asi.snp,by.x="rs_number",by.y="rs_number")
R<-length(which((afr.asi$OR_afr < 1 & afr.asi$OR_asi < 1) | (afr.asi$OR_asi > 1 & afr.asi$OR_afr > 1)))
N<-nrow(afr.asi)
binom.test(R,N,0.5,alternative="greater") 


####################################################################################################
#   Female only X chromosome analysis
####################################################################################################

#  Europ as reference  #############################################################################

#Read in Xwas logistic regression association testing output file.
#XWAS software is available here https://academic.oup.com/jhered/article/106/5/666/2960036
eur<-read.table("russia_fish.xstrat.logistic",stringsAsFactors=FALSE,header=TRUE)
eur<-eur[,c(1,2,3,4,12,13)]
#Required columbs
#SNP: SNP_ID
#BP: Base pair number
#CHR: Chromosome number
#P: P-value

#Run the script 4 times, each time with a different p-value range (different data.p dataset)
#Iteration 1 
data.p<-eur[eur$P<=0.001,]
#Iteration 2
data.p<-eur[eur$P>0.001 & eur$P<=0.01,]
#Iteration 3
data.p<-eur[eur$P>0.01 & eur$P<=0.5,]
#Iteration 4
data.p<-eur[eur$P>0.5,]


#Find list of SNPs within p-value range for each autosomal chromosome within european data
chr<-sort(unique(data.p$CHR))
snp<-c()
j<-1
for(i in 1:length(chr)){
  cdata<-data.p[data.p$CHR==chr[i],]
  cdata<-cdata[order(cdata$BP),]
  
  while(nrow(cdata)>0){
    snp[j]<-cdata$SNP[1]
    r<-which(cdata$BP<=(cdata$BP[1]+500000))
    if(length(r)!=0){
      cdata<-cdata[-r,]
    }
    j<-j+1
  }
}

#Extract selected European SNPs from African and Asian XWAS meta-analysis output file
afr<-read.table("fisher_female_africaSAC.female.meta",stringsAsFactors=FALSE,header=TRUE)
afr<-afr[,c(1,2,3,4,7,9)]
asi<-read.table("fisher_female_china.female.meta",stringsAsFactors=FALSE,header=TRUE)
asi<-asi[,c(1,2,3,4,7,9)]

for (i in 1:length(snp)){
  e<-which(eur$SNP==snp[i])
  af<-which(afr$SNP==snp[i])
  as<-which(asi$SNP==snp[i])
  if (i ==1){
    eur.snp<-eur[e,]
    afr.snp<-afr[af,]
    asi.snp<-asi[as,]
  } else {
    eur.snp<-rbind(eur.snp,eur[e,])
    afr.snp<-rbind(afr.snp,afr[af,])
    asi.snp<-rbind(asi.snp,asi[as,])
  }
}

# Concordance analysis European vs. Africa
eur.snp<-eur.snp[,c(2,4,5,6)]
colnames(eur.snp)<-c("rs_number","reference_allele_eur","OR_eur","p.value_eur")
afr.snp<-afr.snp[,c(3,4,6,5)]
colnames(afr.snp)<-c("rs_number","reference_allele_afr","OR_afr","p.value_afr")
eur.afr<-merge(eur.snp,afr.snp,by.x="rs_number",by.y="rs_number")
R<-length(which((eur.afr$OR_eur < 1 & eur.afr$OR_afr < 1) | (eur.afr$OR_eur > 1 & eur.afr$OR_afr > 1)))
N<-nrow(eur.afr)
binom.test(R,N,0.5,alternative="greater") 

# Concordance analysis European vs. Asia
asi.snp<-asi.snp[,c(3,4,6,5)]
colnames(asi.snp)<-c("rs_number","reference_allele_asi","OR_asi","p.value_asi")
eur.asi<-merge(eur.snp,asi.snp,by.x="rs_number",by.y="rs_number")
R<-length(which((eur.asi$OR_eur < 1 & eur.asi$OR_asi < 1) | (eur.asi$OR_eur > 1 & eur.asi$OR_asi > 1)))
N<-nrow(eur.asi)
binom.test(R,N,0.5,alternative="greater") 

###  Africa as reference  #############################################################################
#Read in XWAS meta-analysis association testing output file. 
afr<-read.table("fisher_female_africaSAC.female.meta",stringsAsFactors=FALSE,header=TRUE)
afr<-afr[,c(1,2,3,4,7,9)]
#Required columbs
#SNP: SNP_ID
#BP: Base pair number
#CHR: Chromosome number
#P: P-value

#Run the script 4 times, each time with a different p-value range (different data.p dataset)
#Iteration 1 
data.p<-afr[afr$P<=0.001,]
#Iteration 2
data.p<-afr[afr$P>0.001 & afr$P<=0.01,]
#Iteration 3
data.p<-afr[afr$P>0.01 & afr$P<=0.5,]
#Iteration 4
data.p<-afr[afr$P>0.5,]

#Find list of SNPs within p-value range for each autosomal chromosome within African data
chr<-sort(unique(data.p$CHR))
snp<-c()
j<-1
for(i in 1:length(chr)){
  cdata<-data.p[data.p$CHR==chr[i],]
  cdata<-cdata[order(cdata$BP),]
  
  while(nrow(cdata)>0){
    snp[j]<-cdata$SNP[1]
    r<-which(cdata$BP<=(cdata$BP[1]+500000))
    if(length(r)!=0){
      cdata<-cdata[-r,]
    }
    j<-j+1
  }
}

#Extract selected African SNPs from European and Asian plink v2 meta-analysis output file
eur<-read.table("russia_fish.xstrat.logistic",stringsAsFactors=FALSE,header=TRUE)
asi<-read.table("fisher_female_china.female.meta",stringsAsFactors=FALSE,header=TRUE)

for (i in 1:length(snp)){
  e<-which(eur$SNP==snp[i])
  af<-which(afr$SNP==snp[i])
  as<-which(asi$SNP==snp[i])
  if (i ==1){
    eur.snp<-eur[e,]
    afr.snp<-afr[af,]
    asi.snp<-asi[as,]
  } else {
    eur.snp<-rbind(eur.snp,eur[e,])
    afr.snp<-rbind(afr.snp,afr[af,])
    asi.snp<-rbind(asi.snp,asi[as,])
  }
}


#Concordance analysis Africa vs. Europe
eur.snp<-eur.snp[,c(2,4,11,12)]
colnames(eur.snp)<-c("rs_number","reference_allele_eur","OR_eur","p.value_eur")
afr.snp<-afr.snp[,c(3,4,6,5)]
colnames(afr.snp)<-c("rs_number","reference_allele_afr","OR_afr","p.value_afr")
eur.afr<-merge(eur.snp,afr.snp,by.x="rs_number",by.y="rs_number")
R<-length(which((eur.afr$OR_eur < 1 & eur.afr$OR_afr < 1) | (eur.afr$OR_eur > 1 & eur.afr$OR_afr > 1)))
N<-nrow(eur.afr)
binom.test(R,N,0.5,alternative="greater") 

#Concordance analysis Africa vs. Asia
asi.snp<-asi.snp[,c(3,4,9,7)]
colnames(asi.snp)<-c("rs_number","reference_allele_asi","OR_asi","p.value_asi")
afr.asi<-merge(afr.snp,asi.snp,by.x="rs_number",by.y="rs_number")
R<-length(which((afr.asi$OR_afr < 1 & afr.asi$OR_asi < 1) | (afr.asi$OR_afr > 1 & afr.asi$OR_asi > 1)))
N<-nrow(afr.asi)
binom.test(R,N,0.5,alternative="greater") 
 
####  Asia as reference  ##############################################################
#Read in XWAS meta-analysis association testing output file. 
asi<-read.table("fisher_female_china.female.meta",stringsAsFactors=FALSE,header=TRUE)
asi<-asi[,c(1,2,3,4,7,9)]
#Required columbs
#SNP: SNP_ID
#BP: Base pair number
#CHR: Chromosome number
#P: P-value

#Run the script 4 times, each time with a different p-value range (different data.p dataset)
#Iteration 1 
data.p<-asi[asi$P<=0.001,]
#Iteration 2
data.p<-asi[asi$P>0.001 & asi$P<=0.01,]
#Iteration 3
data.p<-asi[asi$P>0.01 & asi$P<=0.5,]
#Iteration 4
data.p<-asi[asi$P>0.5,]

#Find list of SNPs within p-value range for each autosomal chromosome within Asian data
chr<-sort(unique(data.p$CHR))
snp<-c()
j<-1
for(i in 1:length(chr)){
  cdata<-data.p[data.p$CHR==chr[i],]
  cdata<-cdata[order(cdata$BP),]
  
  while(nrow(cdata)>0){
    snp[j]<-cdata$SNP[1]
    r<-which(cdata$BP<=(cdata$BP[1]+500000))
    if(length(r)!=0){
      cdata<-cdata[-r,]
    }
    j<-j+1
  }
}

#Extract selected Asian SNPs from European and African XWAS meta-analysis output file
eur<-read.table("russia_fish.xstrat.logistic",stringsAsFactors=FALSE,header=TRUE)
eur<-eur[,c(1,2,3,4,12,13)]
afr<-read.table("fisher_male_africaSAC.male.meta",stringsAsFactors=FALSE,header=TRUE)
afr<-afr[,c(1,2,3,4,7,9)]

for (i in 1:length(snp)){
  e<-which(eur$SNP==snp[i])
  af<-which(afr$SNP==snp[i])
  as<-which(asi$SNP==snp[i])
  if (i ==1){
    eur.snp<-eur[e,]
    afr.snp<-afr[af,]
    asi.snp<-asi[as,]
  } else {
    eur.snp<-rbind(eur.snp,eur[e,])
    afr.snp<-rbind(afr.snp,afr[af,])
    asi.snp<-rbind(asi.snp,asi[as,])
  }
}

#Concordance analysis Asia vs. Europe
eur.snp<-eur.snp[,c(2,4,5,6)]
colnames(eur.snp)<-c("rs_number","reference_allele_eur","OR_eur","p.value_eur")
asi.snp<-asi.snp[,c(3,4,6,5)]
colnames(asi.snp)<-c("rs_number","reference_allele_asi","OR_asi","p.value_asi")
eur.asi<-merge(eur.snp,asi.snp,by.x="rs_number",by.y="rs_number")
R<-length(which((eur.asi$OR_eur < 1 & eur.asi$OR_asi < 1) | (eur.asi$OR_asi > 1 & eur.asi$OR_eur > 1)))
N<-nrow(eur.asi)
binom.test(R,N,0.5,alternative="greater") 

#Concordance analysis Asia vs. Africa
afr.snp<-afr.snp[,c(3,4,6,5)]
colnames(afr.snp)<-c("rs_number","reference_allele_afr","OR_afr","p.value_afr")
afr.asi<-merge(afr.snp,asi.snp,by.x="rs_number",by.y="rs_number")
R<-length(which((afr.asi$OR_afr < 1 & afr.asi$OR_asi < 1) | (afr.asi$OR_asi > 1 & afr.asi$OR_afr > 1)))
N<-nrow(afr.asi)
binom.test(R,N,0.5,alternative="greater") 
########################################################################################################################
