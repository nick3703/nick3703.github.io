data beetle ;
 input w y n ;
 cards ;
 1.6907 6 59
 1.7242 13 60
 1.7552 18 62
 1.7842 28 56
 1.8113 52 63
 1.8369 53 59
 1.8610 61 62
 1.8839 60 60
 ;
 run ;

 proc mcmc data =  beetle outpost = beetpost nmc = 500000 propcov=quanew ;
   parms mu 0 theta2 0 theta3 0 ;

   x = (w - mu)/exp(theta2) ;
   p = (exp(x)/(1 + exp(x)))**exp(theta3) ;
   sigma = exp(theta2) ;
   m1 = exp(theta3) ;

   model y ~ binomial(n,p) ;  
   prior mu ~ general(0) ;
   prior theta2 ~ general(0) ;
   prior theta3 ~ expgamma(0.25,is=0.25) ;
run ;