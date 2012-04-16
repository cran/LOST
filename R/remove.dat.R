remove.dat <-
function(specimen,removes){
ndat<-length(specimen)
rems<-sample(ndat,removes,replace=FALSE)
for (k in 1:removes){
m<-rems[k]
specimen[m]<-NA
}
return(specimen)
}
