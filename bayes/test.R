s <- scan("test_scores.dat")

lpost <- function(parm,data) {
  llik <- sum(dnorm(data,parm,9,log=T))
#  lpri <- dnorm(parm,80,5,log=T)
   lpri <- dgamma(parm,6400/25,80/25,log=T)
  return(llik + lpri + 100)
}

lpost2 <- function(data,parm) {
  sapply(parm,lpost,data)
}

mu <- seq(50,110,length=1000)

plot(mu,lpost2(s,mu),type="l")

post <- function(parm,data) {
  llik <- sum(dnorm(data,parm,9,log=T))
  lpri <- dnorm(parm,80,5,log=T)
  return(exp(llik + lpri + 100)/0.001635394)
}

post2 <- function(parm,data) {
  sapply(parm,post,data)
}

plot(mu,post2(mu,s),type="l")
integrate(post2,60,110,s)
