#################################################################################################################################
##################################    Sex-biased admixture    ###################################################################
#################################################################################################################################
# Analyse diff between males and females useing glm and anova
auto<-read.table('AutoCovar20170316.txt',header=T)
xchr<-read.table('XchrCovar20170316.txt',header=T)
auto.m<-auto[auto$SEX==1,c(5:9,3)]
auto.f<-auto[auto$SEX==2,c(5:9,3)]
xchr.m<-xchr[xchr$SEX==1,c(5:9,3)]
xchr.f<-xchr[xchr$SEX==2,c(5:9,3)]

#update sex 0 for male 1 for female
auto.m$SEX<-0
xchr.m$SEX<-0
auto.f$SEX<-1
xchr.f$SEX<-1

#Male X chromosome vs. female X chromosome

xm.v.xf<-rbind(xchr.m,xchr.f)
fit<-glm(SEX~ AFR.LWK + SAS.GIH + EUR.GBR + SAN.AFR ,family = "binomial", data=xm.v.xf)
fit.null<-glm(SEX~1 ,family = "binomial", data=xm.v.xf)
fit.anova<-anova(fit.null,fit,test="LRT")
fit.anova$`Pr(>Chi)`[2]
summary(fit)
par(mfrow=c(3,2))
plot(fit, which=1:6)

#Male X chromosome vs. female autosome

xm.v.am<-rbind(xchr.m,auto.f)
fit<-glm(SEX~ AFR.LWK + SAS.GIH + EUR.GBR+ SAN.AFR ,family = "binomial", data=xm.v.am)
fit.null<-glm(SEX~1 ,family = "binomial", data=xm.v.am)
fit.anova<-anova(fit.null,fit,test="LRT")
fit.anova$`Pr(>Chi)`[2]
summary(fit)
plot(fit, which=1:6)

#Male autosome vs. female autosome

am.v.af<-rbind(auto.m,auto.f)
fit<-glm(SEX~ AFR.LWK + SAS.GIH + EUR.GBR+ SAN.AFR ,family = "binomial", data=am.v.af)
fit.null<-glm(SEX~1 ,family = "binomial", data=am.v.af)
fit.anova<-anova(fit.null,fit,test="LRT")
fit.anova$`Pr(>Chi)`[2]
summary(fit)
plot(fit, which=1:6)

#Male autosome vs. female X chromosome

xf.v.am<-rbind(auto.m,xchr.f)
fit<-glm(SEX~ AFR.LWK + SAS.GIH + EUR.GBR+ SAN.AFR ,family = "binomial", data=xf.v.af)
fit.null<-glm(SEX~1 ,family = "binomial", data=xf.v.af)
fit.anova<-anova(fit.null,fit,test="LRT")
fit.anova$`Pr(>Chi)`[2]
summary(fit)
plot(fit, which=1:6)

#Male autosome vs. male X chromosome
auto.m$SEX<-1
xm.v.am<-rbind(auto.m,xchr.m)
fit<-glm(SEX~ AFR.LWK + SAS.GIH + EUR.GBR+ SAN.AFR ,family = "binomial", data=xm.v.am)
fit.null<-glm(SEX~1 ,family = "binomial", data=xm.v.am)
fit.anova<-anova(fit.null,fit,test="LRT")
fit.anova$`Pr(>Chi)`[2]
summary(fit)
plot(fit, which=1:6)

#Female autosome vs. Female X chromosome
auto.f$SEX<-0
xf.v.af<-rbind(auto.f,xchr.f)
fit<-glm(SEX~ AFR.LWK + SAS.GIH + EUR.GBR+ SAN.AFR ,family = "binomial", data=xf.v.af)
fit.null<-glm(SEX~1 ,family = "binomial", data=xf.v.af)
fit.anova<-anova(fit.null,fit,test="LRT")
fit.anova$`Pr(>Chi)`[2]
summary(fit)
plot(fit, which=1:6)
