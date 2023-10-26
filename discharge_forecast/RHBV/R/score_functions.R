MSEq<-function(qobs,qsim){
  mean((qobs-qsim)^2)
}

NSq<-function(qobs,qsim){
  1-(sum((qobs-qsim)^2))/(sum((qobs-mean(qobs))^2))
}

LNNSq<-function(qobs,qsim,delta=0.001){
  qobs[!is.na(qobs) & qobs==0]<-delta
  1-(sum((log(qobs)-log(qsim))^2))/(sum((log(qobs)-mean(log(qobs)))^2))
}


VEq<-function(qobs,qsim){
  (sum(qsim)-sum(qobs))/sum(qobs)

}


RBq<-function(qobs,qsim){
  (sum(abs(qobs-qsim)))/sum(qobs)
}

ANEq<-function(qobs,qsim){
  mean(abs(qobs-qsim)/sum(qobs))
}




