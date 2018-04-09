####Tune Bagus method
Tune_Bagus=function(v0,tau,S,n,p_n,eta){
  pb <- txtProgressBar(min = 0, max = length(v0)*4, style = 3)
  bic=NULL
  for(m in 1:length(eta)){
    for(i in 1:length(v0)){
      v1=v0[i]*c(1.5,3,5,10)
      for(j in 1:length(v1)){
       # for(k in 1:length(tau)){
          ##Bayes EM####
          #v1=sqrt(n/(log(p_n)))/5
          #v0=sqrst(n/(p_n*log(p_n)))/5
          
          #tau=sqrt((log(p_n)/n))
          w=1
          l=1
          maxiter=30
          result1<-Bagus(S,n,v0[i],v1[j],maxiter,eta[m],tau=v0[i])
          
          bic=rbind(bic,list(v0=v0[i],v1=v1[j],tau=v0[i],eta=eta[m],BIC=BIC_Bagus(result1$Theta,S,result1$P,n)))
          
          setTxtProgressBar(pb, (m-1)*length(v0)*length(v1)+(i-1)*length(v1)+j)
          
        }
      #} 
    }
    
  }
  
  close(pb)
  return(bic[which.min(bic[,5]),])
}
