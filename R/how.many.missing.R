how.many.missing <-
function (x) 
{
  if(class(x)=="array"){x<-two.d.array(x)}
    totaldata <- nrow(x) * ncol(x)
    excluded <- ifelse(is.na(x), 1, 0)
    sums <- sum(excluded)
    percent <- sums/totaldata
    return(percent)
}
