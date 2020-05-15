library(LearnBayes)

lpost <- function(parm,data) {
  mu <- parm[1]
  sig2 <- parm[2]
  n <- length(data)

  z <- (-n-2)*0.5*log(sig2) - (0.5/(sig2))*((n-1)*var(data) + n*(mean(data) - mu)**2)
  return(z)
}

y <- rnorm(50,0,1)
mycontour(lpost,c(-1,1,0.01,3),y)

s <- rigamma(1000,(length(y)-1)/2,((length(y)-1)/2)*var(y))
m <- rnorm(1000,mean(y),sqrt(s/length(y)))
points(m,s)