library(LearnBayes)

marga <- function(a) {
  n <- 30
  y <- c(0.46162842,0.32546280,0.18797008,0.76845166,0.04301092,0.41355833,
         0.36185887,0.54708582,0.23919662,0.28029385,0.45574949,0.08003376,
         0.42898375,0.95687202,0.46743537,0.03854264,0.33286173,0.22253679,
         0.22063744,0.44068476,0.21077251,0.46399354,0.27928481,0.10385025,
         0.29417748,0.51466444,0.99218757,0.31253086,0.22571566,0.45682085)

  z <- log(a) - a - n*lgamma(a) + lgamma(n*a + 5) - (n*a + 5)*log(sum(y)+1) + (a - 1)*log(prod(y))
  return(exp(z)/189.4804)
}


lmarga <- function(a) {
  return(log(marga(a)))
}

a1 <- rgamma(1000,2,1)
b1 <- rgamma(1000,5,1)
AoverB <- a1/b1
hist(AoverB)
p1 <- rgamma(1000,a1,b1)
hist(p1)

x <- seq(0.05,4.5,length=1000)
plot(x,marga(x),type="l",xlab=expression(alpha),ylab="Marginal Distribution")

y1 <- sample(x,100000,replace=T,prob=marga(x)/sum(marga(x)))
hist(y1,freq=F,xlab=expression(alpha),ylab="Marginal Distribution",main="Brute Force")
lines(x,marga(x))

lap <- laplace(lmarga,2)
hist(y1,freq=F,xlab=expression(alpha),ylab="Marginal Distribution",main="Normal Approximation")
lines(x,dnorm(x,lap$mode,sqrt(lap$var)))

# t distribution variance is sigma^2 df/(df - 2)
plot(x,marga(x),type="l",xlab=expression(alpha),ylab="Marginal Distribution",main="Rejection Sampling")

dtns <- function(x,mu,s2) {
  return(dt((x - mu)/sqrt(s2),4))
}

lines(x,dtns(x,lap$mode,2*lap$var),lty=2)
lines(x,2.45*dtns(x,lap$mode,2*lap$var),lty=2,col="red")

y2a <- sqrt(2*lap$var)*rt(25000,4) + lap$mode
ui <- runif(25000)
keep <- TRUE
keep[y2a < 0] <- FALSE
keep[y2a > 0] <- (ui[y2a > 0] <= (marga(y2a[y2a > 0]))/(2.45*dtns(y2a[y2a > 0],lap$mode,2*lap$var)))
sum(keep)
hist(y2a[keep],freq=F,xlab=expression(alpha),ylab="Marginal Distribution",main="Rejection Sampling")
lines(x,marga(x))

wi <- rep(0,25000)
wi[y2a > 0] <- marga(y2a[y2a > 0])/dtns(y2a[y2a > 0],lap$mode,2*lap$var)
qi <- wi/sum(wi)
hist(qi)
y3 <- sample(y2a,25000,replace=T,prob=qi)

hist(y3,freq=F,xlab=expression(alpha),ylab="Marginal Distribution",main="SIR")
lines(x,marga(x))

y4 <- rep(0,25000)
acc <- 0
y4[1] <- rnorm(1,lap$mode,sqrt(lap$var))
for (i in 2:25000) {
  y4[i] <- y4[i-1]
  ystar <- rnorm(1,y4[i-1],1)
  if (ystar > 0) {
    lnew <- marga(ystar)
    lold <- marga(y4[i-1])
    if (lnew/lold > runif(1)) {
      y4[i] <- ystar
      acc <- acc + 1
    }
  }
}

hist(y4[10001:25000],freq=F,xlab=expression(alpha),ylab="Marginal Distribution",main="Random Walk Metropolis")
lines(x,marga(x))

y5 <- rep(0,25000)
acc <- 0
y5[1] <- rnorm(1,lap$mode,sqrt(lap$var))
for (i in 2:25000) {
  y5[i] <- y5[i-1]
  ystar <- sqrt(2*lap$var)*rt(1,4) + lap$mode
  if (ystar > 0) {
    lnew <- marga(ystar)
    lold <- marga(y5[i-1])
    j1 <- dtns(ystar,lap$mode,2*lap$var)
    j2 <- dtns(y5[i-1],lap$mode,2*lap$var)
    r <- (lnew/lold)*(j2/j1)
    if (r > runif(1)) {
      y5[i] <- ystar
      acc <- acc + 1
    }
  }
}


hist(y5[10001:25000],freq=F,xlab=expression(alpha),ylab="Marginal Distribution",main="Independence Sampler Metropolis")
lines(x,marga(x))



    

