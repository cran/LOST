unbend.spine <-
function(coords,land.pairs,deg=3,restricted=NULL){
  aligned<-bilat.align(coords,land.pairs,average=FALSE,restricted)
  aligned<-aligned[,c(2,3,1)]
  colnames(aligned)<-c("X","Y","Z")
  aligned<-as.data.frame(aligned)
  
  fit<-lm(Y~poly(X,deg),data=aligned)
  new.x<-seq(min(aligned$X),max(aligned$X),length.out=1000)
  
  
  preds<-predict(fit,list(X=new.x))
  
  
  unbent<-coords
  unbent[]<-NA
  for(i in 1:nrow(aligned)){
    new.dat<-cbind(new.x,preds)
    dists<-apply(new.dat,1,function(x){sqrt((x[1]-aligned[i,1])^2+(x[2]-aligned[i,2])^2)})
    x.pos<-which(dists==min(dists))
    
    resids<-dists[x.pos]
    
    pred.current<-predict(fit)[i]
    resids<-resids*ifelse((aligned[i,2]-pred.current)>0,1,-1)
    
    
    x.segs<-vector(length=(x.pos-1))
    for(j in 1:(x.pos-1)){
      x.segs[j]<-sqrt(sum((new.dat[j,]-new.dat[j+1,])^2))   
    }
    
    unbent[i,1]<-sum(x.segs)
    unbent[i,2]<-resids
    unbent[i,3]<-aligned[i,3]
  }
  return(list(bilat.aligned=aligned[,1:3],unbent=unbent))
}
