byclade <-
function (x, remperc, groups) {
  
  ngroups<-length(unique(groups))
  
  
  if(class(x)=="matrix"){
  remove.dat <- function(specimen, removes) {
        ndat <- length(specimen)
        rems <- sample(ndat, removes, replace = FALSE)
        for (k in 1:removes) {
            m <- rems[k]
            specimen[m] <- NA
        }
        return(specimen)
    }
    newx1 <- as.matrix(x)
    grouping <- as.factor(groups)
    newx2 <- as.matrix(x)
    totaldata <- nrow(x) * ncol(x)
    n <- round(totaldata * remperc)
    ndat <- 1:totaldata
    remove <- sample(ndat, n, replace = FALSE)
    for (k in 1:n) {
        i <- remove[k]
        newx1[i] <- NA
    }
    binary <- ifelse(is.na(newx1), 1, 0)
    numberper <- apply(binary, 1, sum)
    rows <- 1:nrow(x)
    numbersp <- ifelse(numberper == 0, 0, 1) * rows
    nsp <- length(numbersp)
    sorted <- sort(numberper, decreasing = TRUE)
    splitgroups <- split(as.data.frame(x), grouping)
    npergroup <- sapply(splitgroups, nrow, simplify = TRUE)
    counts <- 1:nsp
    for (i in 1:nsp) {
        m <- groups[i]
        a <- npergroup[m]
        counts[i] <- a
    }
    counts <- counts
    sums <- sum(npergroup)
    ratio <- sums/counts
    probs <- ratio/sum(ratio)
    orders <- sample(1:nsp, nsp, replace = FALSE, prob = probs)
    for (k in 1:nsp) {
        removes <- sorted[k]
        spnumber <- orders[k]
        specimen <- newx2[spnumber, ]
        newsp <- remove.dat(specimen, removes)
        newx2[spnumber, ] <- newsp
    }
    return(newx2)}
  
  
  
  
  if(class(x)=="array"){
    remove.dat <- function(specimen, removes) {
      ndat <- nrow(specimen)
      rems <- sample(ndat, removes, replace = FALSE)
      for (k in 1:removes) {
        specimen[rems[k],] <- rep(NA,dim(specimen)[[2]])
      }
      return(specimen)
    }
    
    newx1 <- x
    grouping <- as.factor(groups)
    totaldata <- nrow(x) * dim(x)[[3]]
    n <- round(totaldata * remperc)
    all.spl<-cbind(rep(1:dim(x)[3],each=nrow(x)),rep(1:nrow(x),dim(x)[[3]]))
    remove <- all.spl[sample(1:totaldata, n, replace = FALSE),]
    outs <- table(remove[,1])
    remove<-rep(0,dim(x)[[3]])
    remove[as.numeric(names(outs))]<-outs
   
    sorted <- sort(remove, decreasing = TRUE)
    npergroup <- table(groups)
    counts <- rep(0,dim(x)[[3]])
    for (i in 1:length(remove)) {
      m <- groups[i]
      a <- npergroup[m]
      counts[i] <- a
    }
    counts <- counts
    sums <- sum(npergroup)
    ratio <- sums/counts
    probs <- ratio/sum(ratio)
    orders <- sample(1:dim(x)[[3]], dim(x)[[3]], replace = FALSE, prob = probs)
    for (k in 1:length(sorted)){
      removes <- sorted[k]
      spnumber <- orders[k]
      specimen <- newx1[,,spnumber]
      if(removes==0){newsp<-specimen
      } else {
      newsp <- remove.dat(specimen, removes)}
      newx1[,,spnumber] <- newsp
      
    }
    dimnames(newx1)[[3]]<-dimnames(x)[[3]]
    return(newx1)}
  
  
  
  
}
