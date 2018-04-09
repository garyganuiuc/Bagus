####CV for SPACE method
space.cv<-function(Y,rho){
  n=nrow(Y)
  p=ncol(Y)
  id1=1:n
  id=matrix(0,nrow=10,ncol=n/10)
  for(i in 1:10){
    id[i,]=sample(id1,n/10)
    id1=id1[!id1 %in% id[i,]]
  }
  L=0
  L1=NULL
  for(j in 1:length(rho)){
    for(i in 1:10){
      output=space.joint(Y[-id[i,],],lam1=rho[j],iter=5)
      I=output$ParCor*(
        sqrt(matrix(output$sig.fit,ncol=1)%*%matrix(output$sig.fit,nrow=1))/output$sig.fit)
      diag(I)<-0
      RSS<-apply((Y[id[i,],]-Y[id[i,],]%*%I)^2,2,sum)
      
      L=L+n/10*sum(log(RSS))+2*log(n/10)*
        sum(output$ParCor[upper.tri(output$ParCor)]!=0)
    }
    L1=c(L1,L)  
    L=0
  } 
  rho[which.min(L1)]
  
}