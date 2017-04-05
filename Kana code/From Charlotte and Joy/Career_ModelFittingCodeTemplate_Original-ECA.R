#Template to fit single, double, asymptotic models 

#NOTES:
#EC Adair 9/8/2011

#fitting one rep at a time 
#there are two ways that I've used bbmle. One is by writing my own liklihood function
#the other uses bbmle to do this I've shown both versions for the single pool model
#they should give you the same estimates. You can choose which version you prefer. 
#One may be more "stable" or behave better than the other.

#To clear memory
rm(list=ls())

#To find csv file:
file.choose()

#Load aic-liklihood library
library(bbmle)

#read in data
data=read.csv("C:\\Documents and Settings\\adair\\Desktop\\carol\\BioCon\\Data\\Sarah CAREER\\CareerData_Jan09\\CareerMassLossAllHarvests2.csv", header=TRUE)
#look at row headings
head(data)

#vectors for k values (col 1 or 1&2) and sse (col 2)
single.k=array(0,dim=c(max(data$ID),1))
double.k=array(0,dim=c(max(data$ID),3))
asymptotic.k=array(0,dim=c(max(data$ID),3))

#arrays for obs,pred,res
single.pred=array(0,dim=c(nrow(data),4))		#col: 1=Mt,2=t,3=pred,4=res
double.pred=array(0,dim=c(nrow(data),4))		#col: 1=Mt,2=t,3=pred,4=res
asymptotic.pred=array(0,dim=c(nrow(data),6))	#col: 1=Mt,2=t,3=pred,4=res

#arrays to store AICc values
a.aicc=array(0,dim=c(max(data$ID),6)) #3aicc; 3weights; single, asym, double

#Here, models are fit one at a time. You could make this a loop if you 
#are convinced all your models will converge using the same starting values.
#I never have such luck.

#choose rep or set of data points to fit(here identified by the data column ID)

i=1

	s1=subset(data,data$ID==i)
	#browser()
	#create new data frame (with appropriate column names) to put in 
	# only necessary columns, so you don't get messed up
	# by excluding NA's from unneeded columns (e.g., if Final_N data is missing R will 
	# delete the entire row whether you need Final_N for your analysis or not
	t=s1$Years
	Mt=s1$propinit
	x <- data.frame(t, Mt)
	#omit lines with missing data (NAs)
	xNA <-na.exclude(x)

	#Single pool 
		#Method 1: Liklihood function	
		LLS = function(y,k){
	         Mhat=1*exp(-k*xNA$t) #creates a vector of=length to obs of preds
	         ssq = sum((Mhat - xNA$Mt)^2)
		   sigma = sqrt(ssq/length(xNA$Mt))
	         return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
		}

		singleLL = mle2(minuslogl = LLS, start = list(k = 0.5), data = list(y=xNA$Mt),method="L-BFGS-B", lower=c(0),upper=c(100))
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
		
		#Method 2:
		singlepool <- mle2(Mt ~ dnorm(mean = exp(-k * t), sd = sd), 
			data = x, start = list(k = 0.4, sd = .05), control = list(parscale = c(k = 0.01, sd = 0.01) ))
			summary(singlepool)


	#Asymptotic (method 1 only): M(t)=A+(1-A)*exp(-k*t)
		#Liklihood function	
		LLA = function(y,k,A){
         		Mhat=A+(1-A)*exp(-k*xNA$t) #creates a vector of=length to obs of preds
         		ssq = sum((Mhat - xNA$Mt)^2)
         		sigma = sqrt(ssq/length(xNA$Mt))
         		return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
		}

		asymLL = mle2(minuslogl = LLA, start = list(k = .9,A=0.2), data = list(y=xNA$Mt),method="L-BFGS-B", lower=c(0,0),upper=c(100,1))#,control=list(maxint=10e6))
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


	#Double pool (method 1 only): propinitmass = A*exp(-ks*t)+(1-A)*exp(-k*t);
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


#To copy arrays/vectors to excel spreadsheet (just paste into excel after running each line)
#get values after every rep
write.table(predmle2.s,"clipboard-384", sep = "\t",row.names=FALSE,col.names=FALSE) 
write.table(predmle2.a,"clipboard-384", sep = "\t",row.names=FALSE,col.names=FALSE) 
write.table(predmle2.d,"clipboard-384", sep = "\t",row.names=FALSE,col.names=FALSE) 

#get values after all reps are run
write.table(single.k,"clipboard",sep="\t",row.names=FALSE,col.names=FALSE)
write.table(asymptotic.k,"clipboard",sep="\t",row.names=FALSE,col.names=FALSE)
write.table(double.k,"clipboard",sep="\t",row.names=FALSE,col.names=FALSE)
write.table(a.aicc,"clipboard",sep="\t",row.names=FALSE,col.names=FALSE)




