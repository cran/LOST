obliteratorGM <-
function (x, remperc, expo = 1) 
{
    probability.generator <- function(newsp, distances, removes, expo = 1) {
        if (sum(removes) == 0) {
            nrcur <- 1
            current <- 0
        }        else {
            current <- setdiff(removes, 0)
            nrcur <- length(current)
        }
        probmatrix <- matrix(1,ncol=nrow(newsp),nrow=length(removes))
        for (m in 1:(nrcur)) {
       
            if (sum(removes) == 0) {
                anchor <- 0
            } else {
                anchor <- current[m]
            }
            if(!anchor==0){
            inv <- 1/(distances[anchor,]^expo)
            inv[anchor] <- 0
            } else {inv<-rep(1,nrow(newsp))}
            
            sums <- sum(inv)
            ones <- rep(1, length(inv))
            if (anchor == 0) {
                probs <- ones
            }            else {
                probs <- inv/sums
            }
            probs <- ifelse(is.na(newsp[,1]), 0, probs)
            probmatrix[m,] <- probs
        }
        probs <- apply(probmatrix, 2, prod)
        sums <- sum(probs)
        probs <- probs/sums
        return(probs)
    }
    remove.points <- function(specimen, r, distances, expo) {
        removes <- rep(0, r)
        newsp <- specimen
        
        for (k in 1:r) {
            probs <- probability.generator(newsp, distances, 
                removes, expo)
            a <- sample(1:nrow(newsp), 1, prob = probs)
            newsp[a,] <- rep(NA,ncol(newsp))
            removes[k] <- a
        }
        return(newsp)
    }
  
    newx1<-x
    totaldata <- nrow(x) * dim(x)[[3]]
    n <- round(totaldata * remperc)
    all.spl<-cbind(rep(1:dim(x)[3],each=nrow(x)),rep(1:nrow(x),dim(x)[[3]]))
    remove <- all.spl[sample(1:totaldata, n, replace = FALSE),]
    outs <- table(remove[,1])
    remove<-rep(0,dim(x)[[3]])
    remove[as.numeric(names(outs))]<-outs
    
    for (i in 1:dim(x)[[3]]) {
       
        specimen <- x[,,i]
        r <- remove[i]
        if(r==0){next()}
        distances<-as.matrix(dist(specimen,upper=TRUE,diag=TRUE))
        newx1[,,i] <- remove.points(specimen, r, distances, expo)
    }
    dimnames(newx1)[[3]]<-dimnames(x)[[3]]
    return(newx1)
}
