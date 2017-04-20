flipped<-function(specimen,land.pairs,show.plot=FALSE){
  
  reversed<-specimen
  reversed[,1]<- -reversed[,1]
  new.reversed<-reversed

  for(i in 1:nrow(land.pairs)){
  new.reversed[land.pairs[i,1],]<-reversed[land.pairs[i,2],]
  new.reversed[land.pairs[i,2],]<-reversed[land.pairs[i,1],]
}  
  
  sp.rev<-na.omit(cbind(specimen,new.reversed))
  
  pr<-procOPA(sp.rev[,1:3],sp.rev[,4:6])  
  
  transposer<-(sp.rev[,1:3]-pr$Ahat)[1,]  
  trans.mat<-matrix(transposer,ncol=ncol(new.reversed),nrow=nrow(na.omit(new.reversed)),byrow=TRUE)
  
  rot<-fcnt(na.omit(new.reversed))%*%pr$R
  drops<-which(is.na(specimen[which(!is.na(new.reversed[,1])),][,1]))
  if(length(drops)==0){fixed<-rot+trans.mat-matrix(apply((rot-pr$Bhat),2,mean),ncol=3,nrow=nrow(rot),byrow=TRUE)
  } else{
  fixed<-rot+trans.mat-matrix(apply((rot[-drops,]-pr$Bhat),2,mean),ncol=3,nrow=nrow(rot),byrow=TRUE)
}
  
  missing<-which(is.na(new.reversed[,1]))
  for(i in 1:length(missing)){
    fixed<-insertRow(fixed,missing[i],c(NA,NA,NA))
  }

  fills<-which(is.na(specimen[,1]) & !is.na(fixed[,1]))
  specimen.filled<-specimen
  specimen.filled[fills,]<-fixed[fills,]

  if(show.plot){ plot3d(specimen,aspect=FALSE)
    points3d(specimen.filled,col="red")}
  
  return(specimen.filled)
  
  
  }
  
    
  

  
  