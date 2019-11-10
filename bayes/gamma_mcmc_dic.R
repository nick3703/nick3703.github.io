library(coda)
library(LearnBayes)

y <- c(0.46162842,0.32546280,0.18797008,0.76845166,0.04301092,0.41355833,
       0.36185887,0.54708582,0.23919662,0.28029385,0.45574949,0.08003376,
       0.42898375,0.95687202,0.46743537,0.03854264,0.33286173,0.22253679,
       0.22063744,0.44068476,0.21077251,0.46399354,0.27928481,0.10385025,
       0.29417748,0.51466444,0.99218757,0.31253086,0.22571566,0.45682085)

lpost <- function(parm,data) {
  a <- parm[1]
  b <- parm[2] 
#  a1 <- 2
#  a2 <- 1
#  b1 <- 5
#  b2 <- 1
  a1 <- 12
  a2 <- 1
  b1 <- 15
  b2 <- 1
#  a1 <- 0.001
#  a2 <- 0.001
#  b1 <- 0.001
#  b2 <- 0.001
  z <- sum(dgamma(data,a,b,log=T)) + dgamma(a,a1,a2,log=T) + dgamma(b,b1,b2,log=T)
  return(z)
}

mycontour(lpost,c(0.01,5,0.01,12),y,xlab=expression(alpha),ylab=expression(beta))
z <- simcontour(lpost,c(0.01,5,0.01,12),y,10000)
lap <- laplace(lpost,c(2,1),y)

s <- 50000
a <- rep(0,s)
b <- rep(0,s)
acca <- 0
accb <- 0
a[1] <- 2
b[1] <- 5

for (i in 2:s) {
  # Sample from beta
  b[i] <- b[i-1]
  bstar <- rnorm(1,b[i-1],1.5)
  if (bstar > 0) {
    lnew <- lpost(c(a[i-1],bstar),y)
    lold <- lpost(c(a[i-1],b[i-1]),y)
    if (lnew - lold > log(runif(1))) {
      b[i] <- bstar
      accb <- accb + 1
    }
  }

  #Sample from alpha 
  a[i] <- a[i-1]
  astar <- rnorm(1,a[i-1],0.5)
  if (astar > 0) {
    lnew <- lpost(c(astar,b[i]),y)
    lold <- lpost(c(a[i-1],b[i]),y)
    if (lnew - lold > log(runif(1))) {
      a[i] <- astar
      acca <- acca + 1
    }
  }
}

d <- rep(0,49901)
for (i in 100:50000) {
  d[i-99] <- -2*sum(dgamma(y,a[i],b[i],log=T))
}
davg <- mean(d)

dthetahat <- -2*sum(dgamma(y,mean(a[100:50000]),mean(b[100:50000]),log=T))

pd <- davg - dthetahat

dic <- davg + pd


