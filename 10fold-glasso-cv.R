#ten-fold CV to choose penalty for glasso "likelihood approach"
glasso.cv=function(Y,rho){
  n=nrow(Y)
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
                  S=cov(Y[-id[i,],])
                  C=glasso(S,rho[j],penalize.diagonal=TRUE)$wi
                  S1=cov(Y[id[i,],])
                  L=L+log(det(C))-sum(diag(S1%*%C))-rho[j]*sum(abs(C)) 
            }
      L1=c(L1,L)  
      L=0
} 
rho[which.max(L1)]
  #loss
}