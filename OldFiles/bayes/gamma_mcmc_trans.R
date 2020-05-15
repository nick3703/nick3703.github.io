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
chk <- simcontour(lpost,c(0.01,5,0.01,12),y,10000)
lap <- laplace(lpost,c(2,1),y)

k <- t(chol(lap$var))
kinv <- solve(k)
k1 <- kinv[1,1]
k2 <- kinv[1,2]
k3 <- kinv[2,1]
k4 <- kinv[2,2]
plot(k1*(chk$x - lap$mode[1]) + k2*(chk$y - lap$mode[2]),k3*(chk$x - lap$mode[1]) + k4*(chk$y - lap$mode[2]))
plot(k1*chk$x + k2*chk$y,k3*chk$x + k4*chk$y)

lpost2 <- function(parm,data) {
  eta1 <- parm[1]  
  eta2 <- parm[2]
#  a <- (k4*eta1 - k2*eta2)/(k1*k4 - k2*k3) + lap$mode[1]
#  b <- (k3*eta1 - k1*eta2)/(k2*k3 - k1*k4) + lap$mode[2]
   a <- (k4*eta1 - k2*eta2)/(k1*k4 - k2*k3)
   b <- (k3*eta1 - k1*eta2)/(k2*k3 - k1*k4)

  z <- sum(dgamma(data,a,b,log=T)) + dgamma(a,2,1,log=T) + dgamma(b,5,1,log=T)
  return(z)
}
#mycontour(lpost2,c(-3,5,-3.2,5),y,xlab=expression(eta[1]),ylab=expression(eta[2]))
mycontour(lpost2,c(1.9,10,-3,5),y,xlab=expression(eta[1]),ylab=expression(eta[2]))

points(log(chk$x),log(chk$y))
lpost3 <- function(parm,data) {
  theta1 <- parm[1]
  theta2 <- parm[2]
  a <- exp(theta1)
  b <- exp(theta2)

  z <- sum(dgamma(data,a,b,log=T)) + dgamma(a,2,1,log=T) + dgamma(b,5,1,log=T) + theta1 + theta2
  return(z)
}
mycontour(lpost3,c(-0.5,1.5,0.5,2.5),y,xlab=expression(theta[1]),ylab=expression(theta[2]))



