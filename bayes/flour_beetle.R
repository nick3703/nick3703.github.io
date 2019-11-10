library(LearnBayes)
library(coda)

#Raw Data
w <- c(1.6907,1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)
y <- c(6,13,18,28,52,53,61,60)
n <- c(59,60,62,56,63,59,62,60)

#Initialization
niter <- 60000
nstart <- 5
k <- 2.0
mu <- matrix(0,niter,nstart)
theta2 <- matrix(0,niter,nstart)
theta3 <- matrix(0,niter,nstart)
sigma <- matrix(0,niter,nstart)
m1 <- matrix(0,niter,nstart)

logpost <- function(theta,data) {
  mu <- theta[1]
  theta2 <- theta[2]
  theta3 <- theta[3]

  sigma <- exp(theta2)
  m1 <- exp(theta3)
  a0 <- 0.25
  b0 <- 4.0
  y <- data[,1]
  w <- data[,2]
  n <- data[,3]
  
  z <- 0
  for (i in 1:8) {
    x <- (w[i] - mu)/exp(theta2)
    p <- (exp(x)/(1 + exp(x)))**(exp(theta3))
    if (p == 1) return(-Inf)
    z <- z + y[i]*log(p) + (n[i] - y[i])*log(1 - p)
  }
  z <- z - b0*exp(theta3) + a0*theta3
  return(z)
}

#Normal approximation and starting values
lap <- laplace(logpost,c(0,1,1),cbind(y,w,n))
propvar <- k*lap$var
start <- matrix(rmt(nstart,lap$mode,2.0*lap$var,4),ncol=3)

x <- seq(1.78,1.86,length=1000)
plot(x,dnorm(x,lap$mode[1],sqrt(lap$var[1,1])),type="l",xlab=expression(mu),ylab="")

for (i in 1:nstart) {
  #Initialize starting values and acceptance count
  mu[1,i] <- start[i,1]
  theta2[1,i] <- start[i,2]
  theta3[1,i] <- start[i,3]
  ar <- 0

  for (j in 2:niter) {
    old <- c(mu[j-1,i],theta2[j-1,i],theta3[j-1,i])
    mu[j,i] <- old[1]
    theta2[j,i] <- old[2]
    theta3[j,i] <- old[3]
    sigma[j,i] <- exp(theta2[j,i])
    m1[j,i] <- exp(theta3[j,i])

    cand <- rmnorm(1,c(mu[j-1,i],theta2[j-1,i],theta3[j-1,i]),propvar)
    lpo <- logpost(old,cbind(y,w,n))
    lpn <- logpost(as.vector(cand),cbind(y,w,n))

    if (lpn - lpo > log(runif(1))) {
       mu[j,i] <- cand[1,1]
       theta2[j,i] <- cand[1,2]
       theta3[j,i] <- cand[1,3]
       sigma[j,i] <- exp(theta2[j,i])
       m1[j,i] <- exp(theta3[j,i])

       ar <- ar + 1
    }
  }
  print(ar/niter)
}

autocorr.diag(as.mcmc(mu[,1]))
autocorr.plot(as.mcmc(mu[,1]))
autocorr.diag(as.mcmc(theta2[,1]))
autocorr.plot(as.mcmc(theta2[,1]))
autocorr.diag(as.mcmc(theta3[,1]))
autocorr.plot(as.mcmc(theta3[,1]))
summary(as.mcmc(mu[,1]))
summary(as.mcmc(theta2[,1]))
summary(as.mcmc(theta3[,1]))
raftery.diag(as.mcmc(mu[,1]))
raftery.diag(as.mcmc(theta2[,1]))
raftery.diag(as.mcmc(theta3[,1]))

draws <- cbind(c(mu[101:60000,1],mu[101:60000,2],mu[101:60000,3],mu[101:60000,4],mu[101:60000,5]),
  c(exp(theta2[101:60000,1]),exp(theta2[101:60000,1]),exp(theta2[101:60000,1]),exp(theta2[101:60000,1]),exp(theta2[101:60000,1])),
  c(exp(theta3[101:60000,1]),exp(theta3[101:60000,2]),exp(theta3[101:60000,3]),exp(theta3[101:60000,4]),exp(theta3[101:60000,5])))

dose <- seq(1.60,1.90,length=100)
p <- matrix(0,nrow=299500,ncol=100)
for (i in 1:100) {
  x <- (dose[i] - draws[,1])/draws[,2]
  p1 <- exp(x)/(1 + exp(x))
  p[,i] <- p1**draws[,3]
}

plot(dose,apply(p,2,median),xlab="Dose",ylab="P(Death)",type="l",ylim=c(0,1),xlim=c(1.60,1.90))
lines(dose,apply(p,2,quantile,0.05),lty=2)
lines(dose,apply(p,2,quantile,0.95),lty=2)

w50 <- draws[,1] + draws[,2]*log((0.5**(1/draws[,3]))/(1 - 0.5**(1/draws[,3])))
plot(density(w50),xlab="Dose where 50% of flour beetles die after five hours",main="",ylab="Posterior Distribution",xlim=c(1.6,1.9))
quantile(w50,c(0.05,0.5,0.95))