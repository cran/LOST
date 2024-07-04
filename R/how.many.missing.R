how.many.missing <-
function (x) 
{
  if(length(dim(x))==3){x<-two.d.array(x)}
    totaldata <- nrow(x) * ncol(x)
    excluded <- ifelse(is.na(x), 1, 0)
    sums <- sum(excluded)
    percent <- sums/totaldata
    return(percent)
}
