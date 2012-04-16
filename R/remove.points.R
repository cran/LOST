remove.points <-
function(specimen,r,distances,expo){
removes<-rep(0,r)
newsp<-specimen
l<-length(specimen)
nl<-1:(2*l); anchor<-0
site<-c(1:l,1:l)
startorstop<-c(rep(1,l),rep(2,l))
sa<-removes
for (k in 1:r){
probs<-probability.generator(newsp,distances,removes,expo,sa)
a<-sample(nl,1,prob=probs)
b<-site[a]
cc<-startorstop[a]
newsp[b]<-NA
removes[k]<-b
sa[k]<-cc
}
return(newsp)
}
