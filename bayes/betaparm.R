parametersolver = function(qu,p,init) {
# qu is the value of the quantile: median = qu[1]
# p is the quantile: 0.5, 0.25 (median, Q1)
# init is an initial value
  qu <- qu
  p <- p
  betaoptim = function(param) { 
  q1 <- qu[1]
  q2 <- qu[2]
  p1 <- p[1]
  p2 <- p[2]
  (pbeta(q1,param[1],param[2])-p1)^2 + (pbeta(q2,param[1],param[2])-p2)^2
}

r = optim(init,betaoptim)
v = unlist(r)
t = c(v[1],v[2])
print(t)
}

parametersolver(c(.05,.15),c(.5,.9),c(1,1))