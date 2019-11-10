library(LearnBayes)

nd <- scan("newcomb_data.txt")

n <- length(nd)
s2 <- var(nd)
ybar <- mean(nd)
sig2 <- rigamma(10000,0.5*(n-1),0.5*(n-1)*s2)
mu <- rnorm(10000,ybar,sig2/n)

par(mfrow=c(1,1))
mobs <- rep(0,10000)
pred <- rep(0,10000)
rep <- matrix(0,ncol=n,nrow=10000)
nd6 <- sort(nd)[6]
nd61 <- sort(nd)[61]
bp1 <- 0
bp2 <- 0
bp3 <- 0

for (i in 1:10000) {
  rep[i,] <- rnorm(n,mu[i],sqrt(sig2[i]))
  mobs[i] <- min(rep[i,])

  sr <- sort(rep[i,])
  trep1 <- abs(sr[61]-mu[i]) - abs(sr[6] - mu[i])
  ty1 <- abs(nd61 - mu[i]) - abs(nd6 - mu[i])
  bp1 <- bp1 + (trep1 >= ty1)

  trep2 <- sum((rep[i,] - mu[i])^2/sig2[i])
  ty2 <- sum((nd - mu[i])^2/sig2[i])
  bp2 <- bp2 + (trep2 >= ty2)

  trep3 <- -2*sum(dnorm(rep[i,],mu[i],sqrt(sig2[i]),log=T))
  ty3 <- -2*sum(dnorm(nd,mu[i],sqrt(sig2[i]),log=T))
  bp3 <- bp3 + (trep3 >= ty3)
}

print(c(bp1/10000,bp2/10000,bp3/10000))
