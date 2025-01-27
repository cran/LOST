unbend.tps.poly <-function (coords, reference, axes = NULL, deg = 2){
  
  if(length(dim(coords))==2){
    old.coords<-coords
    if(ncol(coords)==2){coords<-cbind(coords,rep(0,nrow(coords)))}
    if (!is.null(axes)) {
      spec <- coords[, axes]
    }    else {
      spec <- coords
    }
    
    
    
    spec <- cbind(spec[, 1], spec[, 2])
    colnames(spec) <- c("X", "Y")
    spec <- as.data.frame(spec)
    fit <- lm(Y ~ poly(X, deg), data = spec[reference, ])
    new.x <- seq(min(spec$X), max(spec$X), length.out = 1000)
    preds <- predict(fit, list(X = new.x))
    unbent <- coords
    unbent[] <- NA
    for (i in 1:nrow(spec)) {
      new.dat <- cbind(new.x, preds)
      dists <- apply(new.dat, 1, function(x) {
        sqrt((x[1] - spec[i, 1])^2 + (x[2] - spec[i, 2])^2)
      })
      x.pos <- which(dists == min(dists))
      resids <- dists[x.pos]
      pred.current <- predict(fit, list(X = spec$X))[i]
      resids <- resids * ifelse((spec[i, 2] - pred.current) > 
                                  0, 1, -1)
      x.segs <- vector(length = (x.pos - 1))
      for (j in 1:(x.pos - 1)) {
        x.segs[j] <- sqrt(sum((new.dat[j, ] - new.dat[j + 
                                                        1, ])^2))
      }
      unbent[i, 1] <- sum(x.segs)
      unbent[i, 2] <- resids
      unbent[i, 3] <- coords[i, setdiff(1:3, axes)]
    }
    if(ncol(old.coords)==2){unbent<-unbent[,-3]}
    return(unbent)
    
    
  }
  
  if(length(dim(coords))==3){
    coords.old<-coords
    if(ncol(coords)==2){
      new.coords<-array(dim=c(nrow(coords),3,dim(coords)[3]))
      new.coords[,1,]<-coords[,1,]
      new.coords[,2,]<-coords[,2,]
      new.coords[,3,]<-0
      coords<-new.coords
    }
    
    unbent<-coords
    unbent[]<-NA
    for(z in 1:dim(coords)[3]){
      
      
      if (!is.null(axes)) {
        spec <- coords[, axes,z]
      }    else {
        spec <- coords[,,z]
      }
      
      
      spec <- cbind(spec[, 1], spec[, 2])
      colnames(spec) <- c("X", "Y")
      spec <- as.data.frame(spec)
      fit <- lm(Y ~ poly(X, deg), data = spec[reference, ])
      new.x <- seq(min(spec$X), max(spec$X), length.out = 1000)
      preds <- predict(fit, list(X = new.x))
      for (i in 1:nrow(spec)) {
        new.dat <- cbind(new.x, preds)
        dists <- apply(new.dat, 1, function(x) {
          sqrt((x[1] - spec[i, 1])^2 + (x[2] - spec[i, 2])^2)
        })
        x.pos <- which(dists == min(dists))
        resids <- dists[x.pos]
        pred.current <- predict(fit, list(X = spec$X))[i]
        resids <- resids * ifelse((spec[i, 2] - pred.current) > 
                                    0, 1, -1)
        x.segs <- vector(length = (x.pos - 1))
        for (j in 1:(x.pos - 1)) {
          x.segs[j] <- sqrt(sum((new.dat[j, ] - new.dat[j + 
                                                          1, ])^2))
        }
        unbent[i, 1, z] <- sum(x.segs)
        unbent[i, 2, z] <- resids
        unbent[i, 3, z] <- coords[i, setdiff(1:3, axes),z]
      }
    }
    if(ncol(coords.old)==2){unbent<-unbent[,-3,]}
    return(unbent)
    
  }
  
  
  
  
} 

