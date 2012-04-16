best.reg <-
function(crocmiss){
estimated.matrix<-matrix(ncol=ncol(crocmiss),nrow=nrow(crocmiss))
for (i in 1:ncol(crocmiss)){
estimated.matrix[,i]<-est.reg1(crocmiss,i)
}
return(estimated.matrix)
}
