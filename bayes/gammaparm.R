parametersolver = function(qu,p,init) {
# I think that the p-th quantile occurs at qu
# parametersolver(c(3.0,4.0),c(.5,.75),c(6,2))
  qu <- qu
  p <- p
  gammaoptim = function(param) { 
  q1 <- qu[1]
  q2 <- qu[2]
  p1 <- p[1]
  p2 <- p[2]
  (pgamma(q1,param[1],param[2])-p1)^2 + (pgamma(q2,param[1],param[2])-p2)^2
}

r = optim(init,gammaoptim)
v = unlist(r)
t = c(v[1],v[2])
print(t)
}

