obliterator <-
function(x,remperc,distances,expo=1) {
newx1<-as.matrix(x)
distances<-as.matrix(distances)
totaldata<-nrow(x)*ncol(x)
n<-round(totaldata*remperc)
ndat<-1:totaldata
remove<-sample(ndat,n,replace=FALSE)
for (k in 1:n) {
  i<-remove[k]
newx1[i]<-NA
}
binary<-ifelse(is.na(newx1),1,0)
numberper<-apply(binary,1,sum)
rows<-1:nrow(x)
numbersp<-setdiff((ifelse(numberper==0,0,1)*rows),0)
nsp<-length(numbersp)
sa<-rep(0,ncol(x))
newx<-x
for (k in 1:nsp){
i<-numbersp[k]
specimen<-x[i,]
r<-numberper[k]
newsp<-remove.points(specimen,r,distances,expo)
newx[i,]<-newsp
}
return(newx)
}
