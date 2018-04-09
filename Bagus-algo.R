####Main algorithm for Bagus Algorithm
Soft_thres<-function(a,b,lambda){
  return(sign(-b/a)*max(abs(-b/a)-lambda/a,0))
}

Cor_desc<-function(Theta_11_inverse,W_22,Theta_12,S_12,v0,v1,P_12,n){
  j=nrow(Theta_11_inverse)
  o=1
  Theta_12_1=10000
  Theta_12_2=Theta_12
  while(max(abs(Theta_12_1-Theta_12))>10^-3&o<=1000){
    o=o+1
    Theta_12_1=Theta_12
    for(i in 1:j){
      b=n/2*((Theta_11_inverse[i,]%*%Theta_12)*W_22-
        (Theta_11_inverse[i,i]%*%Theta_12[i])*W_22+S_12[i,])
      a=n/2*Theta_11_inverse[i,i]*W_22
      lambda=P_12[i]/v1+(1-P_12[i])/v0
      Theta_12[i]=Soft_thres(a=a,b=b,lambda=lambda)
    }
    if(is.nan(sum(Theta_12))|max(abs(Theta_12))>10){
      return(Theta_12_2)
    }
  }
  if(o==1001){
    return(Theta_12_2)
  }else{
    return(Theta_12) 
  }
}
###EM function######
Bagus<-function(S,n,v0,v1,maxiter,eta,tau){
  #####Initialization######
  loss=NULL
  eigen=eigen1=NULL
  p_n=nrow(S)
  #W=S
  #diag(W)=diag(W)+2/n*tau
  #W=diag(diag(W))
  W=diag(1,p_n)
  #Theta=diag(1,p_n)
  Theta=diag(1,p_n)
  o1=0
  Theta1=0
  ##iteration till convergence 
  while(max(abs(Theta1-Theta))>0.01&o1<=maxiter){
    o1=o1+1
    Theta1=Theta
    
    ####E-step######
    ###Calculate P#####
    P=1/(1+v1/v0*exp(-abs(Theta)/v0+abs(Theta)/v1)*(1-eta)/eta) 
    Theta2=0
    o=1
    ###Determine adaptive shrinkage parameter######
    while(max(abs(Theta2-Theta))>0.01&o<=4){
      o=o+1
      Theta2=Theta
      for(i in 1:p_n){
        if(i==1){
          id=(i+1):p_n
        }else if(i!=p_n){
          id=c(1:(i-1),(i+1):p_n)
        }else{
          id=1:(i-1)
        }
        W_11=W[id,id]
        W_12=matrix(W[i,id])
        W_22=W[i,i]
        
        Theta_12=matrix(Theta[i,id])
        Theta_11=Theta[id,id]
        Theta_22=Theta[i,i]
        
        P_12=matrix(P[i,id])
        ####M-step#####
        S_12=matrix(S[i,id])
        
        #W_22=S[i,i]+2/n*tau
        
        Theta_11_inverse=W_11-W_12%*%t(W_12)/W_22
          #W_11-W_12%*%t(W_12)/W_22
        #solve(Theta_11)
        
        W_22=S[i,i]+2/n*tau
        ###Coordinate descent to update theta_12#####
        Theta_12=Cor_desc(Theta_11_inverse,W_22,Theta_12,S_12,v0,v1,P_12,n)
        ###Update theta_22, W#####
        
        Theta_22=1/(W_22)+t(Theta_12)%*%Theta_11_inverse%*%Theta_12
        
        Theta[i,id]=Theta[id,i]=Theta_12
        Theta[i,i]=Theta_22
        
        tmp=as.vector(Theta_22-t(Theta_12)%*%Theta_11_inverse%*%Theta_12)
        tmp1=Theta_11_inverse%*%Theta_12
        W_11=Theta_11_inverse+
          tcrossprod(tmp1)/tmp
        W_12=-tmp1/tmp
        
        W[id,id]=W_11
        W[i,id]=W[id,i]=W_12
        W[i,i]=W_22
      }
      
      
    }
    
    
  }

    
  P=1/(1+v1/v0*exp(-abs(Theta)/(v0)+abs(Theta)/(v1))*(1-eta)/eta)
  
  return(list(Theta=Theta,P=P,tau=tau,p=p,W=W))
}










