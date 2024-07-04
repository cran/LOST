bilat.align<-function(coords,land.pairs, average=TRUE, restricted=NULL){
  
  
  
  
  if(length(dim(coords))==2){
    middle<-coords
    if(!is.null(restricted)){
      middle[]<-rep(apply(coords[restricted,],2,mean),each=nrow(coords))
    } else {
      middle[]<-rep(apply(coords,2,mean),each=nrow(coords))
    }
    coords<-coords-middle
  } 
  
  
  midline<-setdiff(1:dim(coords)[1],unique(c(land.pairs[,1],land.pairs[,2])))
  
  if(!is.null(restricted)){
    
    midline<-intersect(midline,restricted)
    land.pairs<-land.pairs[which(land.pairs[,1] %in% restricted | land.pairs[,2] %in% restricted),]}
  
  
  if(length(dim(coords))==2){
    flat.coords<-array(dim=c(length(midline)+nrow(land.pairs),3,1))
    for(i in 1:nrow(land.pairs)){
      j<-1
      flat.coords[i,,j]<-apply(coords[as.matrix(land.pairs)[i,1:2],],2,mean)
    }
    flat.coords[(i+1):nrow(flat.coords),,]<-coords[midline,]
  } else {
    flat.coords<-array(dim=c(length(midline)+nrow(land.pairs),3,dim(coords)[3]))
    for(j in 1:dim(coords)[3]){
      for(i in 1:nrow(land.pairs)){
        flat.coords[i,,j]<-apply(coords[as.matrix(land.pairs)[i,1:2],,j],2,mean)
      }
    }
    flat.coords[(i+1):nrow(flat.coords),,]<-coords[midline,,]
  }
  
  
  
  optim.bilat<-function(rots,flat.coords,axis.choice){
    rot.x<-matrix(c(1,0,0,0,cos(rots[1]),sin(rots[1]),0,-sin(rots[1]),cos(rots[1])),ncol=3,nrow=3)
    rot.y<-matrix(c(cos(rots[2]),0,-sin(rots[2]),0,1,0,sin(rots[2]),0,cos(rots[2])),ncol=3,nrow=3)
    rot.z<-matrix(c(cos(rots[3]),sin(rots[3]),0,-sin(rots[3]),cos(rots[3]),0,0,0,1),ncol=3,nrow=3)
    new.coords<-flat.coords
    for(i in 1:dim(flat.coords)[3]){
      new.coords[,,i]<-new.coords[,,i]%*%rot.x%*%rot.y%*%rot.z
    }
    bilat.dev<-sum(new.coords[,3,]^2)
    return(bilat.dev)
  }
  
  rots<-optim(par=c(0,0,0),fn=optim.bilat,
              gr="L-BFGS-B",flat.coords=flat.coords,axis.choice=3)$par
  
  
  rot.x<-matrix(c(1,0,0,0,cos(rots[1]),sin(rots[1]),0,-sin(rots[1]),cos(rots[1])),ncol=3,nrow=3)
  rot.y<-matrix(c(cos(rots[2]),0,-sin(rots[2]),0,1,0,sin(rots[2]),0,cos(rots[2])),ncol=3,nrow=3)
  rot.z<-matrix(c(cos(rots[3]),sin(rots[3]),0,-sin(rots[3]),cos(rots[3]),0,0,0,1),ncol=3,nrow=3)
  
  
  aligned.coords<-array(dim=c(nrow(coords),ncol(coords),ifelse(length(dim(coords))==2,1,dim(coords)[3])))
  if(length(dim(coords))==3){
    for(i in 1:dim(coords)[3]){
      aligned.coords[,,i]<-coords[,,i]%*%rot.x%*%rot.y%*%rot.z
    }
  } else {
    aligned.coords<-coords%*%rot.x%*%rot.y%*%rot.z
    
  }
  
  if(length(dim(coords)==3)){
    dimnames(aligned.coords)[3]<-dimnames(coords)[3]
  }
  new.align<-array(dim=dim(flat.coords))
  if(length(dim(coords)==3)){
    dimnames(new.align)[3]<-dimnames(coords)[3]
  }
  if(average){
    for(j in 1:dim(aligned.coords)[3]){
      for(i in 1:nrow(land.pairs)){
        pair.mat<-aligned.coords[c(land.pairs[i,1],land.pairs[i,2]),,j]
        pair.mat[2,3]<-pair.mat[2,3]*(-1)
        new.align[i,,j]<-apply(pair.mat,2,mean)
        
      }
    }
    new.align[(1+nrow(land.pairs)):nrow(new.align),,]<-aligned.coords[midline,,]
    return(new.align)
  } else {
    return(aligned.coords)
  }
}
