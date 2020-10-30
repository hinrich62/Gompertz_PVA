#Population viability analysis of populations using Gompertz production function
#uses the SPS spring/summer Chinook Salmon data 
#Richard A. Hinrichsen, Ph.D.
#rich@hinrichsenenvironmental.com
#(206) 715-2859
#10/29/2020

#Functions to get extinction probability estimates
#estimates2.table.gompertz()
#Get extinction risk estimates for QET=1, QET=10, QET=30 and QET=50


#estimates.table.gompertz()
#Get the extinction risk confidence intervals for a particular QET and RTF
#returns parameter estimates of the Gompertz model, extinction risks and confidence intervals
#and AIC of fitted Gompertz Model

#function to get carrying capacity estiates and confidence intervals
#get.capacity.gompertz()


library(MASS)
library(boot)

#change the below variable PATH to equal the location of the input file
PATH="./SPS_Download_DEC042017.csv"
#use NA for missing values instead of "-99"
sps.tab=read.table(PATH,header=T,as.is=T,sep=",",na.strings = "-99")
iii=sps.tab$FracWild==-99
sps.tab$FracWild[iii]=NA

#limit analysis to Snake River spring/summer Chinook
iii=(sps.tab$ESU=="Salmon, Chinook (Snake River spring/summer-run ESU)")
sps.tab=sps.tab[iii,]

#Limit analysis to 1978 forward
sps.tab=sps.tab[sps.tab$Year>=1978,]

#define some useful variables
sps.tab$p_HOS=1-sps.tab$FracWild
sps.tab$AGE3=sps.tab$Age.3.Returns
sps.tab$AGE4=sps.tab$Age.4.Returns
sps.tab$AGE5=sps.tab$Age.5.Returns
sumage=sps.tab$AGE3+sps.tab$AGE4+sps.tab$AGE5
iii=sumage<=0
sps.tab$AGE3[iii]=sps.tab$AGE4[iii]=sps.tab$AGE5[iii]=NA
iii=is.na(sumage)
sps.tab$AGE3[iii]=sps.tab$AGE4[iii]=sps.tab$AGE5[iii]=NA
sps.tab$S_tot_obs=sps.tab$Spawners
sps.tab$pop=sps.tab$Common.Population.Name
sps.tab$adultspawners=sps.tab$Spawners*(1-sps.tab$AGE3)
sps.tab=sps.tab[order(sps.tab$NMFS_POPID,sps.tab$Year),]
POPIDS=unique(sps.tab$NMFS_POPID)
POPIDS=sort(POPIDS)


#Do run reconstruction for a single population
get.recruits=function(sps2.tab){
 n=dim(sps2.tab)[1]
 recruits=rep(NA,n)
 for(ii in 1:(n-5)){
 recruits[ii]=sps2.tab$AGE3[ii+3]*sps2.tab$S_tot_obs[ii+3]*(1-sps2.tab$p_HOS[ii+3])+
                       sps2.tab$AGE4[ii+4]*sps2.tab$S_tot_obs[ii+4]*(1-sps2.tab$p_HOS[ii+4])+
                       sps2.tab$AGE5[ii+5]*sps2.tab$S_tot_obs[ii+5]*(1-sps2.tab$p_HOS[ii+5])
 }
 return(recruits)
}

n=dim(sps.tab)[1]
sps.tab$recruits=rep(NA,n)
for(ii in 1:length(POPIDS)){
 iii=sps.tab$NMFS_POPID==POPIDS[ii]
 sps2.tab=sps.tab[iii,]
 sps.tab$recruits[iii]=get.recruits(sps2.tab)
}


#get rid of NAs
sps.tab=sps.tab[!is.na(sps.tab$recruits),]
sps.tab=sps.tab[!is.na(sps.tab$S_tot_obs),]


#Gompertz model fit
mydat2=data.frame(Y=log(sps.tab$recruits/sps.tab$adultspawners),LSPAWNERS=log(sps.tab$adultspawners),
   YEARS=as.factor(sps.tab$Year),STOCK=as.factor(sps.tab$pop))
iii=is.infinite(mydat2$Y)|is.na(mydat2$Y)
mydat2=mydat2[!iii,]
res2=glm(Y~-1+STOCK+LSPAWNERS:STOCK+C(YEARS,sum),data=mydat2,na.action=na.omit)



#Estimate extinction probabilities using Gompertz production function
#input variables
#qet = quasi-extinction threshhold
#nyears = extinction time window
#dela = change in Gompertz-a parameter
#rft = reproductive failure threshhold
#NBOOT = number of bootstrap replications (if NBOOT=0, there is no bootstrapping)
#output variable
#mymat= matrix of coefficients (MLEs) from Gompertz model fit
#and extinction probabilities with confidence intervals if NBOOT>0
estimates.table.gompertz=function(COEF=NULL,S2=NULL,qet=50,nyears=100,dela=0.0,rft=10,NBOOT=0)
{
 stocknames=sort(unique(sps.tab$pop))
 nstocks=length(stocknames)
 mymat=summary(res2)$coefficients
 if(!is.null(COEF)){mymat[,1]=COEF} #handle bootstrap
 if(is.null(COEF)){COEF=res2$coefficients}
#mymat contains the parameter estimates
#append extinction probablity and s2 (dispersion)
#to this matrix
 np=dim(mymat)[1]
 mymat2=matrix(NA,ncol=4,nrow=nstocks+1)
 dimnames(mymat2)=list(c("Sigma^2",paste(stocknames,"Prob{Extinct}",sep="")),NULL)
 mymat=rbind(mymat,mymat2)
 iii=grep("YEARS",names(COEF))
 alphastart=COEF[iii]
 alphastart=c(alphastart,-sum(alphastart))
 s2=summary(res2)$dispersion
 if(!is.null(S2)){s2=S2} #boostrap version of s2
 astart=COEF[1:nstocks]
 bstart=COEF[(np-nstocks+1):np]
 mymat[np+1,1]<-s2

#loop over populations
 for(ii in 1:nstocks){
  print(stocknames[ii])
  iii=sps.tab$pop==stocknames[ii]
  sps2.tab=sps.tab[iii,]
  nobs=dim(sps2.tab)[1]

#calculate initial spawners
  SINIT=rep(NA,5)
  lastyear=max(sps2.tab$Year)
  yrs=sps2.tab$Year-(lastyear-5);iii=(yrs>0);yrs=yrs[iii]
  SINIT[yrs]=sps2.tab$adultspawners[iii]
  means=mean(SINIT,na.rm=T)
  iii=is.na(SINIT)
  SINIT[iii]=means


#calculate average age distribution
  sps2.tab$AGE1=rep(0,nobs)
  sps2.tab$AGE2=rep(0,nobs)
  age.mat=cbind(sps2.tab$AGE1,
		sps2.tab$AGE2,
		sps2.tab$AGE3,
		sps2.tab$AGE4,
		sps2.tab$AGE5)
  age=apply(age.mat,MARGIN=c(2),FUN=mean,na.rm=T)

#hold carrying capacity constant
   bstart[ii]=bstart[ii]*(astart[ii]+dela)/astart[ii]
   
   myext=extinct.gompertz(qet=qet,rft=rft,nyears=nyears,
		NTRAJ=5000,SINIT=SINIT,age=age,
		a=astart[ii]+dela,b=bstart[ii],alpha=alphastart,
		s2=s2)
   mymat[np+1+ii,1]=myext

}#populations loop

#next take on bootstrap confidence intervals
 if(NBOOT>0){
#first get the variance-covariance matrix
  nobs=length(residuals(res2))
  k=length(res2$coefficients)
  VCOV=vcov(res2)
  bootmat=matrix(NA,ncol=NBOOT,nrow=nstocks)
  for(ii in 1:NBOOT){
   print(paste("BOOT",ii,"OF",NBOOT))
#get a Monte Carlo Sample from the parameter space
   COEF.NEW=mvrnorm(n=1,mu=res2$coefficients,Sigma=VCOV)
   S2.NEW=s2*rchisq(n=1,df=nobs-k)/(nobs-k)
   res=estimates.table.gompertz(COEF=COEF.NEW,S2=S2.NEW,qet=qet,nyears=nyears,rft=rft,NBOOT=0)
   bootmat[,ii]=res$extmat[,1]
  }#ii
  CI1=apply(bootmat,c(1),quantile,probs=c(.025))
  CI2=apply(bootmat,c(1),quantile,probs=c(.975))
  np=dim(mymat)[1]
  indx=(np-nstocks+1):np
  mymat[indx,2]=CI1
  mymat[indx,3]=CI2
}#if

#get extinction risk estimates
  np=length(mymat[,1])
  iii=(np-nstocks+1):np
  ext.mat=mymat[iii,1:3]
  dimnames(ext.mat)=list(stocknames,c("Pr{ext}","LOWER95","UPPER95"))

  res=list(parmat=mymat[-iii,],extmat=ext.mat,aic=AIC(res2), qet=qet, 
     nyears=nyears,dela=dela,rft=rft,NBOOT=NBOOT)
  return(res)
}


#write table of extinction probabilities
#using the gompertz production function
#output is matrix of extinction probabilites
#with columns representing different QETs
#and rows representing different stock
estimates2.table.gompertz=function()
{
 stocknames=sort(unique(sps.tab$pop))
 nstocks=length(stocknames)
 nyrs=24
 qets=c(1,10,30,50)
 rfts=c(2,10,10,10)
 mymat1=matrix(NA,nrow=nstocks,ncol=length(qets))
 dimnames(mymat1)=list(stocknames,paste("QET=",qets,sep=""))
 i=0
  for(jj in 1:length(qets)){
    print("Progress...")
    print("QET")
    print(qets[jj])
    i=i+1
#uncomment this line for sensitivity analysis of intrinsic productivity
#    res=estimates.table.gompertz(qet=qets[jj],nyears=nyrs,rft=rfts[jj],dela=log(1+.25))
    res=estimates.table.gompertz(qet=qets[jj],nyears=nyrs,rft=rfts[jj],dela=0)
    mymat1[,i]=res$extmat[,1]
}

 return(mymat1)
}

#identity for use in block bootstrap
IDENTITY.fun=function(x){
return(x)
}

#estimate extinction probability using Gompertz model
#input variables
#qet = quasi-extinction threshhold
#rft = reproductive failure thresshold
#nyears = length of extinction window
#NTRAJ = number of random trajectories used to estimate extinction probability
#SINIT = initial spawners used in spawner population trajectories
#age = age distribution of recruits 
#a = Gompertz-a 
#b = Gompertz-b 
#alpha = vector of common year effects
#s2 = variance of residuals errors
extinct.gompertz=function(qet,rft,nyears,NTRAJ,SINIT,age,a,b,alpha,s2)
{
 ext=0
 for(jj in 1:NTRAJ){
  s=rep(0,nyears+10)
  s[1:5]=SINIT
  alphaboot=c(tsboot(tseries=alpha,statistic=IDENTITY.fun,R=1,sim="fixed",l=4,n.sim=nyears+6)$t)
  phi=rnorm(1)*sqrt(s2)+alphaboot[1]
  for(ii in 1:(nyears+5)){
   if(ii>5){
     seq=s[ii]
     iii=(s[(ii-3):ii]<qet)
     if(is.na(s[ii]))return(NA)
     if(sum(iii)==4){ext=ext+1;break}
   }#if
   if(ii<=5){
     seq=s[ii]
   }
   y=a+b*log(s[ii])
   r=seq*exp(y+phi)
   if(s[ii]<rft){r=0}
   if((ii+4)>5)s[ii+4]=s[ii+4]+r*age[4]
   s[ii+5]=s[ii+5]+r*age[5]
   phi=alphaboot[ii+1]+rnorm(1)*sqrt(s2)
  }#for ii
 }#for jj
 ext=ext/NTRAJ
 return(ext)
}

#Get carrying capacity estimate for Gompertz model
#with 95% confidence intervals
get.capacity.gompertz=function(){
 stocknames=sort(unique(sps.tab$pop))
 nstocks=length(stocknames)
 COEF=coef(res2)
 p=length(COEF)
 a=COEF[1:nstocks]
 b=(COEF[(p-nstocks+1):p])
 LCAPACITY=-a/b
 A1=diag(-1/b)
 A2=diag(a/(b*b))
 A3=matrix(0,ncol=nstocks,nrow=p-2*nstocks)
 AMAT=rbind(A1,A3,A2)
 VCOV=vcov(res2)
 VAR=t(AMAT)%*%VCOV%*%AMAT
 PARMAT=matrix(NA,ncol=4,nrow=nstocks)
 PARMAT[,1]=exp(LCAPACITY)
 PARMAT[,2]=sqrt(diag(VAR))
 PARMAT[,3]=LCAPACITY+qnorm(.025)*PARMAT[,2]
 PARMAT[,4]=LCAPACITY+qnorm(.975)*PARMAT[,2]
 PARMAT[,3]=exp(PARMAT[,3])
 PARMAT[,4]=exp(PARMAT[,4])
 dimnames(PARMAT)=list(stocknames,c("Capacity","SE of LOGE(CAPACITY)","LOWER95","UPPER95"))
 return(PARMAT)
}



