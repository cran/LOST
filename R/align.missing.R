align.missing<-function (X) {

completes<-which(apply(ifelse(is.na(X),1,0),3,sum)==0)

if(length(completes)==0){stop("you have zero complete specimens, function cannot be run")}
incompletes<-which(!apply(ifelse(is.na(X),1,0),3,sum)==0)

complete.specimens<-X[,,completes]
incomplete.specimens<-X[,,incompletes]
GPA.com <- procGPA(complete.specimens, pcaoutput = FALSE, distances = FALSE)
X.com<-X

aligned.com <- GPA.com$rotated
mean.aligned <- GPA.com$mshape


for (i in 1:length(incompletes)) {
  cur.spec <- incomplete.specimens[,,i]
  cur.miss <- which(is.na(cur.spec[,1]))
  cur.spec.miss <- cur.spec[-cur.miss, ]
  mean.miss <- mean.aligned[-cur.miss, ]
  OPA <- procOPA(mean.miss, as.matrix(cur.spec.miss))
  transpose <- mean.miss[1, ] - OPA$Ahat[1, ]
  new.specX <- OPA$Bhat[, 1] + transpose[1]
  new.specY <- OPA$Bhat[, 2] + transpose[2]
  if(ncol(X)==3){
    new.specZ <- OPA$Bhat[, 3] + transpose[3]
    new.spec <- cbind(new.specX, new.specY, new.specZ)
    
  } else { new.spec <- cbind(new.specX, new.specY)}
  new.spec.in <- new.spec
  for (m in 1:length(cur.miss)) {
    insert <- cur.miss[m]
    new.spec.in <- insertRow(new.spec.in, insert)
  }
  X.com[,,incompletes[i]] <- new.spec.in
}

X.com[,,completes]<-aligned.com
return(X.com)
}