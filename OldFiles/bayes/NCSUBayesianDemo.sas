*****DEMO 1 SLR****************************************************;

title 'Simple Linear Regression';
 
data Class;
   input Name $ Height Weight @@;
   datalines;
Alfred  69.0 112.5   Alice  56.5  84.0   Barbara 65.3  98.0
Carol   62.8 102.5   Henry  63.5 102.5   James   57.3  83.0
Jane    59.8  84.5   Janet  62.5 112.5   Jeffrey 62.5  84.0
John    59.0  99.5   Joyce  51.3  50.5   Judy    64.3  90.0
Louise  56.3  77.0   Mary   66.5 112.0   Philip  72.0 150.0
Robert  64.8 128.0   Ronald 67.0 133.0   Thomas  57.5  85.0
William 66.5 112.0
;

run;quit;

proc mcmc data=class outpost=classout nmc=10000 thin=2 seed=246810
   mchistory=detailed;
   parms beta0 0 beta1 0;
   parms sigma2 1;
   prior beta0 beta1 ~ normal(mean = 0, var = 1e6);
   prior sigma2 ~ igamma(shape = 3/10, scale = 10/3);
   mu = beta0 + beta1*height;
   model weight ~ n(mu, var = sigma2);
run;quit;


*For Comparison to Frequentist Approach;
title 'SLR - Comparison to PROC REG results';
proc reg data=class;
	model weight=height;
run;quit;




****************************************************************************;
*****DEMO 2 Random-Effects Model;

title 'Random-Effects Model';
 
data heights;
   input Family G$ Height @@;
   datalines;
1 F 67   1 F 66   1 F 64   1 M 71   1 M 72   2 F 63
2 F 63   2 F 67   2 M 69   2 M 68   2 M 70   3 F 63
3 M 64   4 F 67   4 F 66   4 M 67   4 M 67   4 M 69
;
run;quit;

*PROC MCMC does not currently support class statements -- user defines design variables;
data input;
   set heights;
   if g eq 'F' then gf = 1;
   else gf = 0;
   drop g;
run;quit;

proc mcmc data=input outpost=postout nmc=50000 thin=5 seed=7893;
   ods select Parameters REparameters PostSummaries PostIntervals 
      tadpanel;
   parms b0 0 b1 0 s2 1 s2g 1;
 
   prior b: ~ normal(0, var = 10000);
   prior s: ~ igamma(0.01, scale = 0.01); 
   random gamma ~ normal(0, var = s2g) subject=family monitor=(gamma);
   mu = b0 + b1 * gf + gamma;
   model height ~ normal(mu, var = s2);
run;quit;


********************************************************************************************;
*****DEMO 3 Behrens-Fisher Problem;

*Two Independent Samples from two different normal distributions;
*Testing is mu1 = mu2 with the added fun that n1 != n2;

title 'The Behrens-Fisher Problem';
 
data behrens;
   input y ind @@;
   datalines;
121 1  94 1 119 1 122 1 142 1 168 1 116 1
172 1 155 1 107 1 180 1 119 1 157 1 101 1
145 1 148 1 120 1 147 1 125 1 126 2 125 2
130 2 130 2 122 2 118 2 118 2 111 2 123 2
126 2 127 2 111 2 112 2 121 2
;
run;quit;


proc mcmc data=behrens outpost=postout seed=123
          nmc=40000 thin=10 monitor=(_parms_ mudif)
          statistics(alpha=0.01)=(summary interval);
   *ods select PostSummaries PostIntervals;
   parm mu1 0 mu2 0;
   parm sig21 1;
   parm sig22 1;
   prior mu: ~ general(0);
   prior sig21 ~ general(-log(sig21), lower=0);
   prior sig22 ~ general(-log(sig22), lower=0);
   mudif = mu1 - mu2;
   if ind = 1 then
      llike = lpdfnorm(y, mu1, sqrt(sig21));
   else
      llike = lpdfnorm(y, mu2, sqrt(sig22));
   model y ~ general(llike);
run;quit;

proc format;
   value diffmt   low-0 = 'mu1 - mu2 <= 0'
                0<-high = 'mu1 - mu2 > 0';
run;quit;

proc freq data = postout;
   tables mudif /nocum;
   format mudif diffmt.;
run;quit;


***********************************************************************************;
*****DEMO 4 Gelman-Rubin Diagnostics***********************************************;

*Returning to the SLR problem from DEMO 1;

title 'SLR, Gelman-Rubin Diagnostics';

data init;
   input Chain beta0 beta1 sigma2;
   datalines;
   1   10  -5   1
   2  -15  10  20
   3    0   0  50
;
run;quit;

/* define constants */
%let nchain = 3;
%let nparm = 3;
%let nsim = 50000;
%let var = beta0 beta1 sigma2;
 
%macro gmcmc;
   %do i=1 %to &nchain;
      data _null_;
         set init;
         if Chain=&i;
         %do j = 1 %to &nparm;
            call symputx("init&j", %scan(&var, &j));
         %end;
         stop;
      run;
 
      proc mcmc data=class outpost=out&i init=reinit nbi=0 nmc=&nsim
                stats=none seed=7;
         parms beta0 &init1 beta1 &init2;
         parms sigma2 &init3 / n;
         prior beta0 beta1 ~ normal(0, var = 1e6);
         prior sigma2 ~ igamma(3/10, scale = 10/3);
         mu = beta0 + beta1*height;
         model weight ~ normal(mu, var = sigma2);
      run;
   %end;
%mend;


ods html close;
%gmcmc;
ods html;


data all;
   set out1(in=in1) out2(in=in2) out3(in=in3);
   if in1 then Chain=1;
   if in2 then Chain=2;
   if in3 then Chain=3;
run;
 


/* plot the trace plots of three Markov chains. */
%macro trace;
   %do i = 1 %to &nparm;
        proc sgplot data=all cycleattrs;
           series x=Iteration y=%scan(&var, &i) / group=Chain;
        run;
   %end;
%mend;
%trace;


******************************************************************************;
*****DEMO 5 Zero-Inflated Poisson;

title 'Zero-Inflated Poisson';

data roots;
	input photo bap roots n;
	label bap='BAP concentration'
		  photo='Photoperiod'
	      roots='No. of roots'
		  ;
datalines;
8   2.2  0  0
8   4.4  0  0
8   8.8  0  0
8  17.6  0  2
16  2.2  0 15
16  4.4  0 16
16  8.8  0 12
16 17.6  0 19
8   2.2  1  3
8   4.4  1  0
8   8.8  1  0
8  17.6  1  0
16  2.2  1  0
16  4.4  1  2
16  8.8  1  3
16 17.6  1  2
8   2.2  2  2
8   4.4  2  3
8   8.8  2  1
8  17.6  2  0
16  2.2  2  2
16  4.4  2  1
16  8.8  2  2
16 17.6  2  2
8   2.2  3  3
8   4.4  3  0
8   8.8  3  2
8  17.6  3  2
16  2.2  3  2
16  4.4  3  1
16  8.8  3  1
16 17.6  3  4
8   2.2  4  6
8   4.4  4  1
8   8.8  4  4
8  17.6  4  2
16  2.2  4  1
16  4.4  4  2
16  8.8  4  2
16 17.6  4  3
8   2.2  5  3
8   4.4  5  0
8   8.8  5  4
8  17.6  5  5
16  2.2  5  2
16  4.4  5  1
16  8.8  5  2
16 17.6  5  1
8   2.2  6  2
8   4.4  6  3
8   8.8  6  4
8  17.6  6  5
16  2.2  6  1
16  4.4  6  2
16  8.8  6  3
16 17.6  6  4
8   2.2  7  2
8   4.4  7  7
8   8.8  7  4
8  17.6  7  4
16  2.2  7  0
16  4.4  7  0
16  8.8  7  1
16 17.6  7  3
8   2.2  8  3
8   4.4  8  3
8   8.8  8  7
8  17.6  8  8
16  2.2  8  1
16  4.4  8  1
16  8.8  8  0
16 17.6  8  0
8   2.2  9  1
8   4.4  9  5
8   8.8  9  5
8  17.6  9  3
16  2.2  9  3
16  4.4  9  0
16  8.8  9  2
16 17.6  9  2
8   2.2 10  2
8   4.4 10  3
8   8.8 10  4
8  17.6 10  4
16  2.2 10  1
16  4.4 10  3
16  8.8 10  0
16 17.6 10  0
8   2.2 11  1
8   4.4 11  4
8   8.8 11  1
8  17.6 11  4
16  2.2 11  1
16  4.4 11  0
16  8.8 11  1
16 17.6 11  0
8   2.2 12  0
8   4.4 12  0
8   8.8 12  2
8  17.6 12  0
16  2.2 12  1
16  4.4 12  1
16  8.8 12  1
16 17.6 12  0
8   2.2 13  1
8   2.2 17  1
8   4.4 13  1
8   8.8 14  1
8   8.8 14  1
8  17.6 14  1
;
run;

proc sort data=roots;
	by photo bap;
run;

data roots(drop=i n);
	set roots;
	do i = 1 to n;
		output;
	end;
run;

data roots;
   set roots;
   photo_bap = photo*bap;
run;

proc mcmc data=roots diag=all dic propcov=quanew 
      nbi=5000 ntu=5000 nmc=500000 thin=10 plots(smooth)=all seed=27513;
   parms (beta0 beta1 beta2 beta3) 0;
   parms (gamma0 gamma1) 0;
   prior beta: ~ normal(0,var=1000);
   prior gamma: ~ normal(0,var=10);
   mu= exp(beta0 + beta1*photo + beta2*bap + beta3*photo_bap);
   p0= logistic(gamma0 + gamma1*photo);
   llike=log(p0*(roots eq 0) + (1-p0)*pdf("poisson",roots,mu));
   model dgeneral(llike);
   title "Bayesian Analysis of Roots Data Set";
run;quit;

proc mcmc data=roots diag=all dic propcov=quanew 
      ntu=5000 nmc=250000 thin=10 plots(smooth)=all seed=27513;
   parms (beta0 beta1 beta2 beta3) 0;
   prior beta: ~ normal(0,var=1000);
   mu= exp(beta0 + beta1*photo + beta2*bap + beta3*photo_bap);
   model roots~Poisson(mu);
   title "Bayesian Analysis of Roots Data Set";
run;quit;


***************************************************************************************;
*****DEMO 6 GENMOD with BAYES STATEMENT;

title 'PROC GENMOD with BAYES STATEMENT';

data birth;
  input low mother_age mother_wt socio alcohol 
        prev_pretrm hist_hyp uterine_irr phy_visit;
  label low='Indicator for Birth Weight'
        mother_age='Mother''s age'
        mother_wt='Weight at Last Menstrual Period'
        socio='Socio-Economic Status'
        alcohol='Did the mother drink during pregnancy?'
        hist_hyp='History of Hypertension'
        prev_pretrm='Previous Preterm Labors'
        uterine_irr='Uterine Irritability'
        phy_visit='Physician Visit in 1st Trimester';
cards;
 1 28  120  2  1  1  0  1  0
 1 29  130  3  0  0  0  1  1
 1 34  187  1  1  0  1  0  0
 1 25  105  2  0  1  1  0  0
 1 25   85  2  0  0  0  1  0
 1 27  150  2  0  0  0  0  0
 1 23   97  2  0  0  0  1  1
 1 24  128  1  0  1  0  0  1
 1 24  132  2  0  0  1  0  0
 1 21  165  3  1  0  1  0  1
 1 32  105  3  1  0  0  0  0
 1 19   91  3  1  1  0  1  0
 1 25  115  2  0  0  0  0  0
 1 16  130  2  0  0  0  0  1
 1 25   92  3  1  0  0  0  0
 1 20  150  3  1  0  0  0  1
 1 21  200  1  0  0  0  1  1
 1 24  155  3  1  1  0  0  0
 1 21  103  2  0  0  0  0  0
 1 20  125  2  0  0  0  1  0
 1 25   89  2  0  1  0  0  1
 1 19  102  3  0  0  0  0  1
 1 19  112  3  1  0  0  1  0
 1 26  117  3  1  1  0  0  0
 1 24  138  3  0  0  0  0  0
 1 17  130  2  1  1  0  1  0
 1 20  120  1  1  0  0  0  1
 1 22  130  3  1  1  0  1  1
 1 27  130  1  0  0  0  1  0
 1 20   80  2  1  0  0  1  0
 1 17  110  3  1  0  0  0  0
 1 25  105  2  0  1  0  0  1
 1 20  109  2  0  0  0  0  0
 1 18  148  2  0  0  0  0  0
 1 18  110  1  1  1  0  0  0
 1 20  121  3  1  1  0  1  0
 1 21  100  2  0  1  0  0  1
 1 26   96  2  0  0  0  0  0
 1 31  102  3  1  1  0  0  1
 1 15  110  3  0  0  0  0  0
 1 23  187  1  1  0  0  0  1
 1 20  122  1  1  0  0  0  0
 1 24  105  1  1  0  0  0  0
 1 15  115  2  0  0  0  1  0
 1 23  120  2  0  0  0  0  0
 1 30  142  3  1  1  0  0  0
 1 22  130  3  1  0  0  0  1
 1 17  120  3  1  0  0  0  1
 1 23  110  3  1  1  0  0  0
 1 17  120  1  0  0  0  0  1
 1 26  154  2  0  1  1  0  1
 1 20  105  2  0  0  0  0  1
 1 26  190  3  1  0  0  0  0
 1 14  101  2  1  1  0  0  0
 1 28   95  3  1  0  0  0  1
 1 14  100  2  0  0  0  0  1
 1 23   94  2  1  0  0  0  0
 1 17  142  1  0  0  1  0  0
 1 21  130  3  1  0  1  0  1
 0 19  182  1  0  0  0  1  0
 0 33  155  2  0  0  0  0  1
 0 20  105  3  1  0  0  0  1
 0 21  108  3  1  0  0  1  1
 0 18  107  3  1  0  0  1  0
 0 21  124  2  0  0  0  0  0
 0 22  118  3  0  0  0  0  1
 0 17  103  2  0  0  0  0  1
 0 29  123  3  1  0  0  0  1
 0 26  113  3  1  0  0  0  0
 0 19   95  2  0  0  0  0  0
 0 19  150  2  0  0  0  0  1
 0 22   95  2  0  0  1  0  0
 0 30  107  2  0  1  0  1  1
 0 18  100  3  1  0  0  0  0
 0 18  100  3  1  0  0  0  0
 0 15   95  1  0  0  0  0  0
 0 25  118  3  1  0  0  0  1
 0 20  120  2  0  0  0  1  0
 0 28  120  3  1  0  0  0  1
 0 32  121  2  0  0  0  0  1
 0 31  100  3  0  0  0  1  1
 0 36  202  3  0  0  0  0  1
 0 28  120  2  0  0  0  0  0
 0 25  120  2  0  0  0  1  1
 0 28  167  3  0  0  0  0  0
 0 17  122  3  1  0  0  0  0
 0 29  150  3  0  0  0  0  1
 0 26  168  1  1  0  0  0  0
 0 17  113  1  0  0  0  0  1
 0 17  113  1  0  0  0  0  1
 0 24   90  3  1  1  0  0  1
 0 35  121  1  1  1  0  0  1
 0 25  155  3  0  0  0  0  1
 0 25  125  1  0  0  0  0  0
 0 29  140  3  1  0  0  0  1
 0 19  138  3  1  0  0  0  1
 0 27  124  3  1  0  0  0  0
 0 31  215  3  1  0  0  0  1
 0 33  109  3  1  0  0  0  1
 0 21  185  1  1  0  0  0  1
 0 19  189  3  0  0  0  0  1
 0 23  130  1  0  0  0  0  1
 0 21  160  3  0  0  0  0  0
 0 18   90  3  1  0  0  1  0
 0 18   90  3  1  0  0  1  0
 0 32  132  3  0  0  0  0  1
 0 19  132  2  0  0  0  0  0
 0 24  115  3  0  0  0  0  1
 0 22   95  2  1  0  0  0  0
 0 22  120  3  0  0  1  0  1
 0 23  128  2  0  0  0  0  0
 0 22  130  3  1  0  0  0  0
 0 30   95  3  1  0  0  0  1
 0 19  115  2  0  0  0  0  0
 0 16  110  2  0  0  0  0  0
 0 21  110  2  1  0  0  1  0
 0 30  153  2  0  0  0  0  0
 0 20  103  2  0  0  0  0  0
 0 17  119  2  0  0  0  0  0
 0 17  119  2  0  0  0  0  0
 0 23  119  2  0  0  0  0  1
 0 24  110  2  0  0  0  0  0
 0 28  140  3  0  0  0  0  0
 0 26  133  2  1  1  0  0  0
 0 20  169  2  0  1  0  1  1
 0 24  115  2  0  0  0  0  1
 0 28  250  2  1  0  0  0  1
 0 20  141  3  0  1  0  1  1
 0 22  158  1  0  1  0  0  1
 0 22  112  3  1  1  0  0  0
 0 31  150  2  1  0  0  0  1
 0 23  115  2  1  0  0  0  1
 0 16  112  1  0  0  0  0  0
 0 16  135  3  1  0  0  0  0
 0 18  229  1  0  0  0  0  0
 0 25  140  3  0  0  0  0  1
 0 32  134  3  1  1  0  0  1
 0 20  121  1  1  0  0  0  0
 0 23  190  3  0  0  0  0  0
 0 22  131  3  0  0  0  0  1
 0 32  170  3  0  0  0  0  0
 0 30  110  2  0  0  0  0  0
 0 20  127  2  0  0  0  0  0
 0 23  123  2  0  0  0  0  0
 0 17  120  2  1  0  0  0  0
 0 19  105  2  0  0  0  0  0
 0 23  130  3  0  0  0  0  0
 0 36  175  3  0  0  0  0  0
 0 22  125  3  0  0  0  0  1
 0 24  133  3  0  0  0  0  0
 0 21  134  2  0  0  0  0  1
 0 19  235  3  1  0  1  0  0
 0 25   95  3  1  1  0  1  0
 0 16  135  3  1  0  0  0  0
 0 29  135  3  0  0  0  0  1
 0 29  154  3  0  0  0  0  1
 0 19  147  3  1  0  0  0  0
 0 19  147  3  1  0  0  0  0
 0 30  137  3  0  0  0  0  1
 0 24  110  3  0  0  0  0  1
 0 19  184  3  1  0  1  0  0
 0 24  110  2  0  1  0  0  0
 0 23  110  3  0  0  0  0  1
 0 20  120  2  0  0  0  0  0
 0 25  241  1  0  0  1  0  0
 0 30  112  3  0  0  0  0  1
 0 22  169  3  0  0  0  0  0
 0 18  120  3  1  0  0  0  1
 0 16  170  1  0  0  0  0  1
 0 32  186  3  0  0  0  0  1
 0 18  120  2  0  0  0  0  1
 0 29  130  3  1  0  0  0  1
 0 33  117  3  0  0  0  1  1
 0 20  170  3  1  0  0  0  0
 0 28  134  2  0  0  0  0  1
 0 14  135  3  0  0  0  0  0
 0 28  130  2  0  0  0  0  0
 0 25  120  3  0  0  0  0  1
 0 16   95  2  0  0  0  0  1
 0 20  158  3  0  0  0  0  1
 0 26  160  2  0  0  0  0  0
 0 21  115  3  0  0  0  0  1
 0 22  129  3  0  0  0  0  0
 0 25  130  3  0  0  0  0  1
 0 31  120  3  0  0  0  0  1
 0 35  170  3  0  1  0  0  1
 0 19  120  3  1  0  0  0  0
 0 24  116  3  0  0  0  0  1
 0 45  123  3  0  0  0  0  1
;
run;quit;


proc genmod data=birth desc;
    model low = alcohol hist_hyp mother_wt prev_pretrm
                / dist=binomial link=logit;
	bayes seed=27513;
    title 'Bayesian Analysis of Low Birth Weight Model';
run;quit;


data Prior;
   input _TYPE_ $ alcohol1 alcohol0 hist_hyp mother_wt prev_pretrm;
datalines;
Mean 1.0986 0 0 0 0
Var 0.00116 1e6 1e6 1e6 1e6 
;
run;

proc genmod data=birth desc;
    class alcohol(desc);
    model low = alcohol hist_hyp mother_wt prev_pretrm
                / dist=binomial link=logit;
	lsmeans alcohol / diff oddsratio plots=all cl;
    bayes seed=27513 coeffprior=normal(input=Prior)  
          outpost=bayes_prob1 plots(smooth)=all diag=all;
    title 'Bayesian Analysis of Low Birth Weight Model';
run;


**Extend number of iterations to improve convergence statistics;
proc genmod data=birth desc;
    class alcohol(desc);
    model low = alcohol hist_hyp mother_wt prev_pretrm
                / dist=binomial link=logit;
	lsmeans alcohol / diff oddsratio plots=all cl;
    bayes seed=27513 coeffprior=normal(input=Prior) nmc=25000  
          outpost=bayes_prob1 plots(smooth)=all diagnostics=all;
    title 'Bayesian Analysis of Low Birth Weight Model';
run;
