library(glasso)
library(MASS)
library(space)
library(Matrix)
library(PDSCE)
library(flare)

######Preload all the source below.######
source("FPandFN.R")
source("Bagus-algo.R")
source("tune-Bagus.R")
source("BIC_Bagus.R")
source("10fold-glasso-cv.R")
source("space-cv.R")

#####Generate setting for sample inverse covariance matrix######
###case test
p_n=p=30
n=50
###case 1
p_n=p=50 #number of variables
n=100 #number of observations

###case 2
p_n=p=100
n=100

###case 3
n=100
p_n=p=200

##AR(2)####
C=toeplitz(c(1,0.5,0.25,rep(0,p_n-3)))
Sigma=solve(C)


###circle case#####
C=toeplitz(c(2,1,rep(0,p_n-3),0.9))
Sigma = solve(toeplitz(c(2,1,rep(0,p_n-3),0.9)))


###star case####
C=diag(1,p)
C[1,-1]=C[-1,1]=1/sqrt(p_n)
diag(C)=1
C=C*1
Sigma=solve(C)


###Random Select 5% as edges
sigma=3
C=diag(1,p_n)*sigma
#set.seed(92)
id2=sample(1:(p_n*(p_n-1)),1.5*p_n)
#set.seed(415)
C[-c(1:p_n,1:p_n)][id2]=sample(c(-1,1),length(id2),replace=TRUE)*runif(length(id2),0.4,1)*sigma
for(i in 1:p_n){
  if(sum(abs(C[i,-i]))>0){
    C[i,-i]=C[i,-i]/(1.1*sum(abs(C[i,-i]))/C[i,i])
  }
}
C=(t(C)+C)/2
diag(C)=sigma
Sigma=solve(C)


###selection accuracy Init#########

MCC1=Sens1=Spec1=Fnorm1=matrix(NA,nrow=1,ncol=50)
MCC2=Sens2=Spec2=Fnorm2=matrix(NA,nrow=1,ncol=50)
MCC3=Sens3=Spec3=Fnorm3=matrix(NA,nrow=1,ncol=50)
MCC4=Sens4=Spec4=Fnorm4=matrix(NA,nrow=1,ncol=50)
MCC5=Sens5=Spec5=Fnorm5=matrix(NA,nrow=1,ncol=50)
MCC6=Sens6=Spec6=Fnorm6=matrix(NA,nrow=1,ncol=50)

####tune parameters for Bagus
eta=0.5
v0=tau=sqrt(1/(n*log(p)))*c(0.4,2,4,20)


for(i in 1:50){
  
  #generate samples from multivariate Gaussian distribution
  Y<-mvrnorm(n,rep(0,p_n),Sigma)
  S<-cov(Y)  #sample covariance
  
  ####the approach of Friedman glasso#####
  rho1=glasso.cv(Y,rho=c(0.07,0.1,0.3,0.2)) #the approach of Friedman
  output1=glasso(S,rho1,penalize.diagonal=FALSE)
  
  #####the approach of Meinhausen-Buhlmann######
  alpha=0.1
  #rho2=neighbor.cv(Y,rho=n^(-1/2)*qnorm(1-alpha/(2*p^2))*c(0.01,0.1,0.5,1,2,5,n))
  rho2=n^(-1/2)*qnorm(1-alpha/(2*p^2))*norm(Y,"F")/n
  #l1=1/sqrt(n)*qnorm(1-alpha/(2*p^2))
  output2=space.neighbor(Y,lam1=rho2)  #the approach of Meinhausen-Buhlmann
  
  #####SPACE##########
  iter=5
  rho3=space.cv(Y,rho=n^(1/2)*qnorm(1-alpha/(2*p^2))*c(0.01,0.03,0.05,0.07,0.1,0.5,0.7,1,1.5,2,5,10))
  #l2=n^(3/2)*qnorm(1-alpha/(2*p^2))
  output3=space.joint(Y,lam1=rho3,iter=iter)  #SPACE by Peng

  Fnorm1[i]<-norm(output1$wi-C,type="F")
  Fnorm2[i]<-norm(output2$ParCor*output2$sig.fit-C,type="F")
  
  Fnorm3[i]<-norm(output3$ParCor*output3$sig.fit-C,type="F")

  
  Spec1[i]<-TN(output1$wi,C)/(TN(output1$wi,C)+FP(output1$wi,C))
  Spec2[i]<-TN(output2$ParCor,C)/(TN(output2$ParCor,C)+FP(output2$ParCor,C))
  Spec3[i]<-TN(output3$ParCor,C)/(TN(output3$ParCor,C)+FP(output3$ParCor,C))
 
  Sens1[i]<-TP(output1$wi,C)/(TP(output1$wi,C)+FN(output1$wi,C))
  Sens2[i]<-TP(output2$ParCor,C)/(TP(output2$ParCor,C)+FN(output2$ParCor,C))
  Sens3[i]<-TP(output3$ParCor,C)/(TP(output3$ParCor,C)+FN(output3$ParCor,C))
 
  
  MCC1[i]<-(TP(output1$wi,C)*TN(output1$wi,C)-FP(output1$wi,C)*FN(output1$wi,C))/
    (sqrt((TP(output1$wi,C)+FP(output1$wi,C)))*sqrt(TP(output1$wi,C)+FN(output1$wi,C))*
       sqrt(TN(output1$wi,C)+FP(output1$wi,C))*sqrt(TN(output1$wi,C)+FN(output1$wi,C)))
  MCC2[i]<-(TP(output2$ParCor,C)*TN(output2$ParCor,C)-FP(output2$ParCor,C)*FN(output2$ParCor,C))/
    (sqrt(TP(output2$ParCor,C)+FP(output2$ParCor,C))*sqrt(TP(output2$ParCor,C)+FN(output2$ParCor,C))*
       sqrt(TN(output2$ParCor,C)+FP(output2$ParCor,C))*sqrt(TN(output2$ParCor,C)+FN(output2$ParCor,C)))
  MCC3[i]<-(TP(output3$ParCor,C)*TN(output3$ParCor,C)-FP(output3$ParCor,C)*FN(output3$ParCor,C))/
    (sqrt(TP(output3$ParCor,C)+FP(output3$ParCor,C))*sqrt(TP(output3$ParCor,C)+FN(output3$ParCor,C))*
       sqrt(TN(output3$ParCor,C)+FP(output3$ParCor,C))*sqrt(TN(output3$ParCor,C)+FN(output3$ParCor,C)))
 
  
  ######Bagus#######
  eta=0.5
 Tune=Tune_Bagus(v0,tau,S,n,p_n,eta)
  maxiter=20
  v0_t=Tune$v0
  v1_t=Tune$v1
  tau_t=Tune$tau
  result1<-Bagus(S,n,v0_t,v1_t,maxiter,eta,tau_t)
  
  P1=diag(1,p_n)
  P1[result1$P>0.5]=1

  Fnorm4[i]<-norm(result1$Theta-C,type="F") 
  Spec4[i]<-TN1(result1$P,C)/(TN1(result1$P,C)+FP1(result1$P,C))
  
  Sens4[i]<-TP1(result1$P,C)/(TP1(result1$P,C)+FN1(result1$P,C))
  
  MCC4[i]<-(TP1(result1$P,C)*TN1(result1$P,C)-FP1(result1$P,C)*FN1(result1$P,C))/
    (sqrt(TP1(result1$P,C)+FP1(result1$P,C))*sqrt(TP1(result1$P,C)+FN1(result1$P,C))*
       sqrt(TN1(result1$P,C)+FP1(result1$P,C))*sqrt(TN1(result1$P,C)+FN1(result1$P,C)))
  
  
  
  
  ##clime estimator()
    out=sugm(Y,method="clime",lambda.min.ratio=0.005,nlambda=30,verbose=FALSE)
  out1.select = sugm.select(out, criterion = "cv",verbose = FALSE)
  clime=out1.select$opt.icov
  
  Fnorm5[i]<-norm(clime-C,type="F")
  
  Spec5[i]<-TN(clime,C)/(TN(clime,C)+FP(clime,C))
  
  Sens5[i]<-TP(clime,C)/(TP(clime,C)+FN(clime,C))
  
  
  MCC5[i]<-(TP(clime,C)*TN(clime,C)-FP(clime,C)*FN(clime,C))/
    (sqrt((TP(clime,C)+FP(clime,C)))*sqrt(TP(clime,C)+FN(clime,C))*
       sqrt(TN(clime,C)+FP(clime,C))*sqrt(TN(clime,C)+FN(clime,C)))
   
  
}


###Results
round(c(mean(Spec1),sd(Spec1),mean(Sens1),sd(Sens1),mean(MCC1),sd(MCC1)),digits=3)
round(c(mean(Spec2),sd(Spec2),mean(Sens2),sd(Sens2),mean(MCC2),sd(MCC2)),digits=3)
round(c(mean(Spec3),sd(Spec3),mean(Sens3),sd(Sens3),mean(MCC3),sd(MCC3)),digits=3)
round(c(mean(Spec4),sd(Spec4),mean(Sens4),sd(Sens4),mean(MCC4),sd(MCC4)),digits=3)
round(c(mean(Spec5),sd(Spec5),mean(Sens5),sd(Sens5),mean(MCC5),sd(MCC5)),digits=3)

round(c(mean(Fnorm1),sd(Fnorm1)),digits=3)
round(c(mean(Fnorm2),sd(Fnorm2)),digits=3)
round(c(mean(Fnorm3),sd(Fnorm3)),digits=3)
round(c(mean(Fnorm4),sd(Fnorm4)),digits=3)
round(c(mean(Fnorm5),sd(Fnorm5)),digits=3)

round(c(mean(Fnorm6),sd(Fnorm6)),digits=3)

P1_estimate=P1_trace/50






