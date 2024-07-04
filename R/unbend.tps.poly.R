unbend.tps.poly <-
function(coords,reference,axes=NULL,deg=3){
  
  if(!is.null(axes)){
    spec<-coords[,axes]
  } else {spec<-coords
  }
  
  spec<-cbind(spec[,1],spec[,2])
  colnames(spec)<-c("X","Y")
  spec<-as.data.frame(spec)
  
  
  fit<-lm(Y~poly(X, deg),data=spec[reference,])
  new.x<-seq(min(spec$X),max(spec$X),length.out=1000)
  
  preds<-predict(fit,list(X=new.x))
  
  
  unbent<-coords
  unbent[]<-NA
  for(i in 1:nrow(spec)){
    new.dat<-cbind(new.x,preds)
    dists<-apply(new.dat,1,function(x){sqrt((x[1]-spec[i,1])^2+(x[2]-spec[i,2])^2)})
    x.pos<-which(dists==min(dists))
    resids<-dists[x.pos]
    pred.current<-predict(fit,list(X=spec$X))[i]
    resids<-resids*ifelse((spec[i,2]-pred.current)>0,1,-1)
    x.segs<-vector(length=(x.pos-1))
    for(j in 1:(x.pos-1)){
      x.segs[j]<-sqrt(sum((new.dat[j,]-new.dat[j+1,])^2))   
    }
    
    unbent[i,1]<-sum(x.segs)
    unbent[i,2]<-resids
    unbent[i,3]<-coords[i,setdiff(1:3,axes)]
  }
  return(unbent)
}
