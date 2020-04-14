MissingGeoMorph<-function (x, method = "BPCA", original.scale=FALSE) {
  
  if(!requireNamespace("pcaMethods")){
    print("some packages must be downloaded from bioconductor, use setRepositories() to select 'BioC' options")
  }
  
  if(ncol(x)==3 & method=="TPS"){
    warning("method TPS can only support 2d data, see package geomorph for a 3d implementation")
  }
  
  
  tps2d <- function(M, matr, matt) {
    p <- dim(matr)[1]
    q <- dim(M)[1]
    n1 <- p + 3
    P <- matrix(NA, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        r2 <- sum((matr[i, ] - matr[j, ])^2)
        P[i, j] <- r2 * log(r2)
      }
    }
    P[which(is.na(P))] <- 0
    Q <- cbind(1, matr)
    L <- rbind(cbind(P, Q), cbind(t(Q), matrix(0, 3, 3)))
    m2 <- rbind(matt, matrix(0, 3, 2))
    coefx <- ginv(L) %*% m2[, 1]
    coefy <- ginv(L) %*% m2[, 2]
    fx <- function(matr, M, coef1) {
      Xn <- numeric(q)
      for (i in 1:q) {
        Z <- apply((matr - matrix(M[i, ], p, 2, byrow = T))^2, 
                   1, sum)
        U <- sum(coef1[1:p] * (Z * log(Z + 1)))
        U[which(is.na(U))] <- 0
        Xn[i] <- coef1[p + 1] + coef1[p + 2] * M[i, 1] + 
          coef1[p + 3] * M[i, 2] + U
      }
      Xn
    }
    matg <- matrix(NA, q, 2)
    matg[, 1] <- fx(matr, M, coefx)
    matg[, 2] <- fx(matr, M, coefy)
    return(matg)
  }
  
  
  
  missing.tps <- function(X) {
    completes<-which(apply(ifelse(is.na(X),1,0),3,sum)==0)
    incompletes<-which(!apply(ifelse(is.na(X),1,0),3,sum)==0)
    
    complete.specimens<-X[,,completes]
    incomplete.specimens<-X[,,incompletes]
    GPA.com <- procGPA(complete.specimens, pcaoutput = FALSE, distances = FALSE)
    X.com<-X
    
    mean.aligned <- GPA.com$mshape
    
    for (i in 1:length(incompletes)) {
      cur.spec <- incomplete.specimens[,,i]
      cur.miss <- which(is.na(cur.spec[,1]))
      cur.spec.miss <- cur.spec[-cur.miss, ]
      mean.miss <- mean.aligned[-cur.miss, ]
      tps.aligned <- tps2d(mean.aligned, mean.miss, cur.spec.miss)
      X.com[,,incompletes[i]] <- tps.aligned
      
    }
    X.com[,,completes]<-complete.specimens
    return(X.com)
  }
  
  
  est.reg1 <- function(X, col_indep) {
    cols <- ncol(X)
    rows <- nrow(X)
    estimated_matrix <- matrix(ncol = cols, nrow = rows)
    if (col_indep == 1) {
      deps <- 2:cols
    }
    else if (col_indep == cols) {
      deps <- 1:(cols - 1)
    }
    else {
      deps <- c(1:(col_indep - 1), (col_indep + 1):cols)
    }
    ndeps <- length(deps)
    indep <- X[, col_indep]
    rcoefs <- numeric()
    for (i in 1:ndeps) {
      vari <- deps[i]
      lm_fit <- lm(X[, vari] ~ indep)
      lm_sum <- summary.lm(lm_fit)
      rcoefs[i] <- sqrt(lm_sum$r.squared)
    }
    ranks <- rank(rcoefs)
    newindep <- indep
    for (m in 1:length(deps)) {
      a <- m - 1
      whichvari <- ifelse(ranks == (length(rcoefs) - a), 
                          1, 0)
      strongest <- sum(deps * whichvari)
      indeps_lm <- lm(newindep ~ X[, strongest])
      indeps_coef <- indeps_lm$coefficients
      logestimate_indep <- indeps_coef[1] + indeps_coef[2] * 
        (X[, strongest])
      estimate_indep <- logestimate_indep
      missings <- ifelse(is.na(newindep), 1, 0)
      nonmissings <- ifelse(is.na(newindep), 0, 1)
      fillnonmissing <- ifelse(is.na(newindep), 0, newindep)
      fillmissing <- ifelse(is.na(newindep), estimate_indep, 
                            0)
      newindep <- fillmissing + fillnonmissing
    }
    return(newindep)
  }
  best.reg <- function(x) {
    estimated.matrix <- matrix(ncol = ncol(x), nrow = nrow(x))
    for (i in 1:ncol(x)) {
      estimated.matrix[, i] <- est.reg1(x, i)
    }
    return(estimated.matrix)
  }
  
  
  
  if (method == "BPCA") {
    aligned <- align.missing(x)
    new.matrix <- two.d.array(aligned)
    estimator <- pca(new.matrix, method = "bpca")
    estimated.values <- completeObs(estimator)
    results <- arrayspecs(estimated.values, p=nrow(x),k=ncol(x))
  }
  else if (method == "mean") {
    aligned <- align.missing(x)
    new.matrix <- two.d.array(aligned)
    estimated.values <- impute(new.matrix, what = "mean")
    results <- arrayspecs(estimated.values, p=nrow(x),k=ncol(x))
  }
  else if (method == "reg") {
    aligned <- align.missing(x)
    new.matrix <- two.d.array(aligned)
    estimated.values <- best.reg(new.matrix)
    results <- arrayspecs(estimated.values, p=nrow(x),k=ncol(x))
  }
  else if (method == "TPS") {
    results <- missing.tps(aligned)
  }
  
  
if(original.scale==TRUE){
results2<-x  
incoms<-which(apply(ifelse(is.na(x),1,0),3,sum)>0)

    for(j in 1:length(incoms)){
    olds<-x[,,incoms[j]]
    news<-results[,,incoms[j]]
    which.com<-which(!is.na(olds[,1])) 
    procs<-procOPA(olds[which.com,],news[which.com,])
    transposer<-(olds[which.com,]-procs$Ahat)[1,]  
    trans.mat<-matrix(transposer,ncol=ncol(news),nrow=nrow(news),byrow=TRUE)
    rot<-fcnt(news)%*%procs$R*procs$s
    fixed<-rot+trans.mat-matrix((rot[which.com,]-procs$Bhat)[1,],ncol=ncol(x),nrow=nrow(rot),byrow=TRUE)
    results2[,,incoms[j]]<-fixed
  }
  return(results2)
} else {
  return(results)}
  
  
}