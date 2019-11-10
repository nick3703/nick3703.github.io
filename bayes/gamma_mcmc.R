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
  z <- sum(dgamma(data,a,b,log=T)) + dgamma(a,2,1,log=T) + dgamma(b,5,1,log=T)
  return(z)
}

mycontour(lpost,c(0.01,5,0.01,12),y,xlab=expression(alpha),ylab=expression(beta))
z <- simcontour(lpost,c(0.01,5,0.01,12),y,10000)
lap <- laplace(lpost,c(2,1),y)

s <- 5000
a <- rep(0,s)
b <- rep(0,s)
acc <- 0
a[1] <- 5
b[1] <- 12

for (i in 2:s) {
  # Sample from beta
  b[i] <- rgamma(1,length(y)*a[i-1] + 5,1 + sum(y))

  #Sample from alpha 
  a[i] <- a[i-1]
  astar <- rnorm(1,a[i-1],0.5)
  if (astar > 0) {
    lnew <- lpost(c(astar,b[i]),y)
    lold <- lpost(c(a[i-1],b[i]),y)
    if (lnew - lold > log(runif(1))) {
      a[i] <- astar
      acc <- acc + 1
    }
  }
}

print(acc/s)

plot(as.mcmc(a))
plot(as.mcmc(b))








