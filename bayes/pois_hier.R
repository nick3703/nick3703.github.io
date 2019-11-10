marg <- function(y,a,b) {
  l <- a*log(b) - lgamma(a) - lfactorial(y) + lgamma(a + y) - (a + y)*log(b + 1)
  return(exp(l))
}

joint <- function(y1,y2,a,b) {
  l <- a*log(b) - lgamma(a) - lfactorial(y1) - lfactorial(y2) + 
       lgamma(a + y1 + y2) - (a + y1 + y2)*log(b + 2)
  return(exp(l))
}

cond <- function(y1,y2,a,b) {
  l <-  lgamma(a + y1 + y2) + (a + y1)*log(b + 1) - lfactorial(y2) + 
        - (a + y1 + y2)*log(b + 2) - lgamma(a + y1) 
  return(exp(l))
}

plot(0:40,cond(0:40,3,5,1),type="h",xlab=expression(Y[1]),ylab=expression(paste("P(",Y[2],"=3|",Y[1],")")))

