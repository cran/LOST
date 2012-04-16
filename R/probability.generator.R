probability.generator <-
function(newsp,distances,removes,expo=1,sa){
zeros<-0
if (sum(removes)==0){nrcur<-1;current=0} else
{current<-setdiff(removes,zeros); nrcur<-length(current)}
probmatrix<-rep(1,2*length(newsp))
for (m in 1:(nrcur)){
distancesx1<-distances[1,]
distancesx2<-distances[2,]
distancesy1<-distances[3,]
distancesy2<-distances[4,]
distancesz1<-distances[5,]
distancesz2<-distances[6,]
distancesx<-rbind(distancesx1,distancesx2)
distancesy<-rbind(distancesy1,distancesy2)
distancesz<-rbind(distancesz1,distancesz2)
if (sum(removes)==0){anchor<-0} else
{anchor<-current[m]}
if (anchor==0){sss<-1} else
{sss<-sa[m]}
if (anchor==0){basex<-0} else
{basex<-distancesx[sss,anchor]}
if (anchor==0){basey<-0} else
{basey<-distancesy[sss,anchor]}
if (anchor==0){basez<-0} else
{basez<-distancesz[sss,anchor]}
distsx1<-sqrt((distancesx1-basex)^2)
distsx2<-sqrt((distancesx2-basex)^2)
#
distsy1<-sqrt((distancesy1-basey)^2)
distsy2<-sqrt((distancesy2-basey)^2)
#
distsz1<-sqrt((distancesz1-basez)^2)
distsz2<-sqrt((distancesz2-basez)^2)
#
distsstart<-sqrt(distsx1^2+distsy1^2+distsz1^2)
distsstop<-sqrt(distsx2^2+distsy2^2+distsz2^2)
dists<-c(distsstart,distsstop)
#
nozeros<-ifelse(dists==0,(max(dists)*10),dists)
dists<-ifelse(dists==0,(min(nozeros)/2),dists)
inv<-1/(dists^expo)
inv[anchor]<-0
ll<-length(inv)/2
inv[(anchor+ll)]<-0
sums<-sum(inv)
ones<-rep(1,length(inv))
if (anchor==0){probs<-ones} else
{probs<-inv/sums}
checker<-c(newsp,newsp)
probs<-ifelse(is.na(checker),0,probs)
probmatrix<-rbind(probmatrix,probs)
}
probs<-apply(probmatrix,2,prod)
sums<-sum(probs)
probs<-probs/sums
return(probs)
}
