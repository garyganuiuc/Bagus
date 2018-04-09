#number of false positive(falsely identified as edge) for Bagus
FP1<-function(C1,C){
  id=which(abs(C[upper.tri(C)])<0.00001)
  id1=which(C1[upper.tri(C1)]>0.5)
  length(which(id1%in%id))
}

#number of FN(missed edge) for Bagus
FN1<-function(C1,C){
  id=which(abs(C[upper.tri(C)])>=0.00001)
  id1=which(C1[upper.tri(C1)]<0.5)
  length(which(id1%in%id))
}

###number of TP for Bagus
TP1<-function(C1,C){
  id=which(abs(C[upper.tri(C)])>=0.00001)
  id1=which(C1[upper.tri(C1)]>=0.5)
  length(which(id1%in%id)) 
}

##number of TN for Bagus
TN1<-function(C1,C){
  id=which(abs(C[upper.tri(C)])<0.00001)
  id1=which(C1[upper.tri(C1)]<0.5)
  length(which(id1%in%id)) 
}
#number of false positive(falsely identified as edge)
FP<-function(C1,C){
  id=which(abs(C[upper.tri(C)])<0.00001)
  id1=which(C1[upper.tri(C1)]!=0)
  length(which(id1%in%id))
}

#number of FN(missed edge)
FN<-function(C1,C){
  id=which(abs(C[upper.tri(C)])>0.00001)
  id1=which(C1[upper.tri(C1)]==0)
  length(which(id1%in%id))
}

###number of TP 
TP<-function(C1,C){
  id=which(abs(C[upper.tri(C)])>0.00001)
  id1=which(C1[upper.tri(C1)]!=0)
  length(which(id1%in%id)) 
}

##number of TN
TN<-function(C1,C){
  id=which(abs(C[upper.tri(C)])<0.00001)
  id1=which(C1[upper.tri(C1)]==0)
  length(which(id1%in%id)) 
}