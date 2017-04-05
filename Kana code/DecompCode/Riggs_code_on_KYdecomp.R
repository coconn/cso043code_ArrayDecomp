library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)
library(tidyr)
library(bbmle)

## bringing in data
pathdata <- "~/_Berkeley/_thesis/R stats/"
# path for christine
pathdata <- "/Users/christine/Desktop/RESEARCH PROJECTS/cso043_ArrayDecomp/Kana code/DecompCode/"

## bring in decomposition data
decomp <- as.data.frame(read.csv(paste(pathdata,"decomp_weights_Rtest02",".csv",sep=""), stringsAsFactors=FALSE))

## change variables to factors
decomp$TopographicLocation <- as.factor(decomp$TopographicLocation)
decomp$CollectionGroup <- as.factor(decomp$CollectionGroup)
decomp$DateDeployed <- as.POSIXct(decomp$DateDeployed, format="%Y.%m.%d")
decomp$DateCollected <- as.POSIXct(decomp$DateCollected, format="%Y.%m.%d")

## add time columns
decomp <- mutate(decomp, Sec = DateCollected - DateDeployed)
decomp <- mutate(decomp, Days = Sec/(3600*24))
decomp <- mutate(decomp, Year = Days/365)


# check out data
str(decomp)

#####################################################################

#Code to fit single, double, asymptotic models on pooled replicates (for model selection)

#vectors for k values (col 1 or 1&2) and A (col 2 or 3)
single.k=array(0,dim=c(max(decomp$ID),1)) 
double.k=array(0,dim=c(max(decomp$ID),3)) 
asymptotic.k=array(0,dim=c(max(decomp$ID),2))


#arrays to store AICc values
a.aicc=array(0,dim=c(max(decomp$ID),6)) #3aicc; 3weights; single, asym, double


### alternatively
# single.k=array(0,dim=c(max(decomp$TopographicLocation),1)) 
# double.k=array(0,dim=c(max(decomp$TopographicLocation),3)) 
# asymptotic.k=array(0,dim=c(max(decomp$TopographicLocation),2))
# 
# 
# #arrays to store AICc values
# a.aicc=array(0,dim=c(max(decomp$TopographicLocation),6)) #3aicc; 3weights; single, asym, double


#fit models in a loop
for (i in 1:4){
  
  s1=subset(decomp,decomp$ID==i)
  
  
  t=s1$Year
  t_orig=s1$Year
  t=as.numeric(t_orig)
  Mt=s1$percent_remaining
  x <- data.frame(t, Mt)
  #omit lines with missing data (NAs)
  xNA <-na.exclude(x)
  
  #Single pool: M(t)=1*exp(-k*t)
  LLS = function(y,k){
    Mhat=1*exp(-k*xNA$t) #creates a vector of=length to obs of preds
    ssq = sum((Mhat - xNA$Mt)^2)
    sigma = sqrt(ssq/length(xNA$Mt))
    return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
  }
  
  singleLL = mle2(minuslogl = LLS, start = list(k = .5), data = list(y=xNA$Mt),method="L-BFGS-B", lower=c(0),upper=c(100))
  #add 1 to K (in AIC caculation) for estimating sigma
  attr(singleLL ,"df") = attr(singleLL,"df")+1
  summary(singleLL)
  
  #Extract k value
  single.k[i,1]=coef(singleLL)[1]
  
  #single pool predictions
  if(i==1){predmle2.s=1*exp(-single.k[i,1]*x[,1])}
  if(i>1){
    predmle2temp=1*exp(-single.k[i,1]*x[,1])
    predmle2.s=append(predmle2.s,predmle2temp)
  }
  
  
  #Asymptotic: M(t)=A+(1-A)*exp(-k*t)
  #Liklihood function  
  LLA = function(y,k,A){
    Mhat=A+(1-A)*exp(-k*xNA$t) #creates a vector of=length to obs of preds
    ssq = sum((Mhat - xNA$Mt)^2)
    sigma = sqrt(ssq/length(xNA$Mt))
    return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
  }
  
  asymLL = mle2(minuslogl = LLA, start = list(k = .9,A=0.2), data = list(y=xNA$Mt),method="L-BFGS-B", lower=c(0,0),upper=c(100,1))#,control=list(maxint=10e6))
  #I get errors here as I hit bounds; modify starting parameters
  attr(asymLL ,"df") = attr(asymLL,"df")+1
  summary(asymLL)
  
  #Extract k & A value
  asymptotic.k[i,1]=coef(asymLL)[1]
  asymptotic.k[i,2]=coef(asymLL)[2]
  
  #Asymptotic preds
  if(i==1){predmle2.a=asymptotic.k[i,2]+(1-asymptotic.k[i,2])*exp(-asymptotic.k[i,1]*x[,1])}
  if(i>1){
    predmle2temp=asymptotic.k[i,2]+(1-asymptotic.k[i,2])*exp(-asymptotic.k[i,1]*x[,1])
    predmle2.a=append(predmle2.a,predmle2temp)
  }
  
  #Double pool: propinitmass = A*exp(-ks*t)+(1-A)*exp(-k*t);
  #Liklihood function	
  LLD = function(y,ks,k,A){
    Mhat=A*exp(-ks*xNA$t)+(1-A)*exp(-k*xNA$t) #creates a vector of=length to obs of preds
    ssq = sum((Mhat - xNA$Mt)^2)
    #browser()
    sigma = sqrt(ssq/length(xNA$Mt))
    return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
  }
  
  start=list(k=1,ks=.3,A=.4)
  doubleLL = mle2(minuslogl = LLD, start = start, data = list(y=xNA$Mt),method="L-BFGS-B", lower=c(0,0,0),upper=c(100,100,1),control=list(maxit=100, trace=TRUE, parscale=abs(unlist(start))))
  #I get errors here as I hit bounds; modify starting parameters
  attr(doubleLL ,"df") = attr(singleLL,"df")+1
  summary(doubleLL)
  
  #Extract k value; ks&A, k&(1-A)
  double.k[i,1]=coef(doubleLL)[1]#ks
  double.k[i,2]=coef(doubleLL)[2]#k
  double.k[i,3]=coef(doubleLL)[3]#A
  
  #Double preds
  if(i==1){predmle2.d=double.k[i,3]*exp(-double.k[i,1]*x[,1])+(1-double.k[i,3])*exp(-double.k[i,2]*x[,1])}
  if(i>1){
    predmle2temp=double.k[i,3]*exp(-double.k[i,1]*x[,1])+(1-double.k[i,3])*exp(-double.k[i,2]*x[,1])
    predmle2.d=append(predmle2.d,predmle2temp)
  }
  
  #Extract AICC values
  xx<-AICctab(singleLL,asymLL,doubleLL, nobs=nrow(x), sort=FALSE, delta=TRUE, weights=TRUE)
  xx
  a.aicc[i,1]=xx$dAICc[1]
  a.aicc[i,2]=xx$dAICc[2]
  a.aicc[i,3]=xx$dAICc[3]
  a.aicc[i,4]=xx$weight[1]
  a.aicc[i,5]=xx$weight[2]
  a.aicc[i,6]=xx$weight[3]
  a.aicc[i,]
}

#get values after all reps are run
write.csv(single.k,file="KY_singlek_modelselection.csv")
write.csv(asymptotic.k,file="KY_asympk_modelselection.csv")
write.csv(double.k,file="KY_doublek_modelselection.csv")
write.csv(a.aicc,file="KY_aic_modelselection.csv")


#evaluate AIC output; select model that minimizes AIC; use that model to fit parameters
#for whole dataset (below)






#########################################

#Code to fit asymptotic model for each species x habitat x site combination


#vectors for k values (col 1) and A (col 2)
asymptotic.k=array(0,dim=c(max(willow$ID_sihasp),2)) 

#arrays for obs,pred,res
#predmle2.a=array(0,dim=c(max(willow$ID_sihasp),4))	#col: 1=Mt,2=t,3=pred,4=res	


#fit models in a loop
for (i in 1:198){
  
  s1=subset(willow,willow$ID_sihasp==i)
  t=s1$Years
  Mt=s1$PropInitMass
  x <- data.frame(t, Mt)
  #omit lines with missing data (NAs)
  xNA <-na.exclude(x)
  
  #Asymptotic: M(t)=A+(1-A)*exp(-k*t)
  #Liklihood function  
  LLA = function(y,k,A){
    Mhat=A+(1-A)*exp(-k*xNA$t) #creates a vector of=length to obs of preds
    ssq = sum((Mhat - xNA$Mt)^2)
    sigma = sqrt(ssq/length(xNA$Mt))
    return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
  }
  
  asymLL = mle2(minuslogl = LLA, start = list(k = 1.4,A=0.5), data = list(y=xNA$Mt),method="L-BFGS-B", lower=c(0,0),upper=c(100,1))#,control=list(maxint=10e6))
  #I get errors here as I hit bounds; modify starting parameters
  #10/21/13: changed start from list(k = .9,A=0.2) to list(k = 1.4,A=0.5)
  attr(asymLL ,"df") = attr(asymLL,"df")+1
  summary(asymLL)
  
  #Extract k & A value
  asymptotic.k[i,1]=coef(asymLL)[1]
  asymptotic.k[i,2]=coef(asymLL)[2]
  
  #Asymptotic preds
  if(i==1){predmle2.a=asymptotic.k[i,2]+(1-asymptotic.k[i,2])*exp(-asymptotic.k[i,1]*x[,1])}
  if(i>1){
    predmle2temp=asymptotic.k[i,2]+(1-asymptotic.k[i,2])*exp(-asymptotic.k[i,1]*x[,1])
    predmle2.a=append(predmle2.a,predmle2temp)
  }
  
}

#get values after all reps are run
write.csv(predmle2.a,"predmle2a.csv")  

#get values after all reps are run
write.csv(asymptotic.k,"asymp_kvalues_alldata.csv")