library(LearnBayes)

p <- seq(0,1,length=10000)

plot(p,histprior(p,c(0.125,0.375,0.625,0.875),c(1/7,2/7,0,4/7)),xlab="p",ylab="Prior",type="l")

plot(p,dbeta(p,6,2),type="l",xlab="p",ylab="Likelihood")

post <- function(p) {
  histprior(p,c(0.125,0.375,0.625,0.875),c(1/7,2/7,0,4/7))*dbeta(p,6,2)
}

plot(p,post(p))
abline(v=c(0.25,0.5,0.75),lty=2)
integrate(post,0,1)

post <- function(p) {
  histprior(p,c(0.125,0.375,0.625,0.875),c(1/7,2/7,0,4/7))*dbeta(p,6,2)
}
plot(p,post(p))
abline(v=c(0.25,0.5,0.75),lty=2)
integrate(post,0,1)

probpost <- post(p)
y <- sample(p,100000,replace=T,prob=probpost)
hist(y,freq=F,main="100,000 Posterior Samples",xlab="",ylab="")
lines(p,post(p)/integrate(post,0,1)$value)
quantile(y,seq(0.05,0.95,by=0.05))