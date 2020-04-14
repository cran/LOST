missing.data <-function (x, remperc=NULL,remsp=NULL,land.vec=NULL,land.identity=NULL) {

 
if(class(x)=="array"){

######### remove a percentage of data ####################  
  if(is.null(remsp)){
      if(is.null(land.identity)){
      totaldata <- dim(x)[[1]] * dim(x)[[3]]
      n <- round(totaldata * remperc)
      ndat <- 1:totaldata
      slindex<-cbind(rep(1:dim(x)[[1]],dim(x)[[3]]),rep(1:dim(x)[[3]],each=dim(x)[[1]]))
      remove <- sample(ndat, n, replace = FALSE)
      
      for (k in 1:n) {
        i <- slindex[remove[k],][1]
        j <- slindex[remove[k],][2]
        
        x[i,,j] <- rep(NA,dim(x)[[2]])
      }
      return(x)
      } else {
        
################# remove a percentage of data, but only from certain landmarks ################
            lands<-length(unique(land.identity))
            tot.possible<-lands*dim(x)[[3]]
              
              if((tot.possible/length(x))<remperc){stop(paste("Warning, the chosen landmarks represent less than ",100*remperc,"%",
                                                   "of the data available, choose a lower value for remperc.",sep=" "))}
              
            n <- round(dim(x)[[1]] * dim(x)[[3]] * remperc)
            ndat <- 1:tot.possible
            slindex<-cbind(rep(land.identity,dim(x)[[3]]),rep(1:dim(x)[[3]],each=length(land.identity)))
            remove <- sample(ndat, n, replace = FALSE)
            
            for (k in 1:n) {
              i <- slindex[remove[k],][1]
              j <- slindex[remove[k],][2]
              x[i,,j] <- rep(NA,dim(x)[[2]])
            }
        return(x)
        }
      
} else { 
  
  
####### remove landmarks from a percentage of specimens##############  
        if(is.null(land.identity)){
      removesps <- sample(1:dim(x)[[3]], round(dim(x)[[3]] * remsp), replace = FALSE)
      removels<- sample(land.vec,length(removesps),replace=TRUE)
      
      for (k in 1:length(removesps)) {
       spec<-removesps[k]
        n.ls<-removels[k]
        remove<-sample(1:nrow(x),n.ls,replace=TRUE)
        for(m in 1:length(remove)){
          i<-remove[m]
          x[i,,spec] <- rep(NA,dim(x)[[2]])
            }
        
      }
      return(x)
    
    } else {
      
################ remove certain landmarks from a percentage of specimens##############      
      
    lands<-length(unique(land.identity))
    if(max(land.vec)>lands){stop("Warning, the chosen landmarks represent less than the maximum value of land.vec")}
  
    removesps <- sort(sample(1:dim(x)[[3]], round(dim(x)[[3]] * remsp), replace = FALSE))
    removels<- sample(land.vec,length(removesps),replace=TRUE)
    
    
     for (k in 1:length(removesps)) {
      spec<-removesps[k]
      n.ls<-removels[k]
      remove<-sample(land.identity,n.ls,replace=FALSE)
      for(m in 1:length(remove)){
        i<-remove[m]
        x[i,,spec] <- rep(NA,dim(x)[[2]])
      }
      
    }
      
    
    return(x)
  }
}
      
} else {
  
  x<-as.matrix(x)
    totaldata <- nrow(x) * ncol(x)
    n <- round(totaldata * remperc)
    ndat <- 1:totaldata
    remove <- sample(ndat, n, replace = FALSE)
    for (k in 1:n) {
      i <- remove[k]
      x[i] <- NA
      
    }
    return(x)
  }
}
     

