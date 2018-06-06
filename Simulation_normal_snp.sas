
/*	Purpose: Simulation study  									*/
/*	Paper: Mediation analysis for longitudinal data using 		*/
/*         regression calibration when the mediator is measured */
/*         with error  											*/
/*	By: John Ssenkusu 								    		*/
/*  Date last modified: February 13, 2017			    		*/

* Notes;
* 1. Modify the 'home' libname to point to where you want the simulated results stored ;
* 2. Modify the path of 'slmm.mac' to point to the directory where you've stored the SNP code ;

*libname home 'C:\..............\Simulation' filelockwait=15; /* Reference the directory where results will be saved */

%macro simulate(SimNum);

options nonotes;
%do i=1 %to &SimNum;

%let a_N_norm=;
%let b_N_norm=;
%let ADE_N_norm=;

%let a_SN_norm=;
%let b_SN_norm=;
%let ADE_SN_norm=;

%let a_MN_norm=;
%let b_MN_norm=;
%let ADE_MN_norm=;

ods graphics off;
ods exclude all;
ods noresults;

%let TotalSubjects = 600;  							/* size of sample*/
%let NumberOfVisits = 6; 
 
data Datasim_1(keep=x_i bi_N bi_MN bi_SN id time bi_Y );
array prob [2] _temporary_ (0.65 0.35); 			/* mixing probabilities for mixture of normals*/
p = 1/2; 											/* proportion exposed */
id = 0;
do i = 1 to &TotalSubjects; 
	id = id + 1;

	x_i = rand("Bernoulli", p);
	bi_N = rand("Normal",0 , 2);       				/* random numbers from N(0,1) */
    
	Type = rand("Table", of prob[*]); 				/* returns 1, or 2 */
	if Type=1 then bi_MN = rand("Normal", 0, 1);	/* Generating a mixture of normals*/
	else bi_MN = rand("Normal", 4, 1);
	
    bi_SN = rand('chisquare', 4);
	bi_Y = rand("normal", 0, 1);					/* Random intercepts for y*/

do k = 1 to &NumberOfVisits;
	time = k;
	output;
	end;
	end;
run;	

/* Generate errors for the mediator and outcome from a multivariate normal distribution */

proc iml;
obs_num = &NumberOfVisits*&TotalSubjects;
norm_mean = {0, 0};
norm_cov = {1 0, 
            0 1.5};
x = RandNormal(obs_num, norm_mean, norm_cov );
id = colvec(repeat(T(1:&TotalSubjects), 1, &NumberOfVisits));
Z = ID || x;
create errors from Z[c={"ID" "e_y" "e_m"}];
append from Z;
close errors;
quit;

/* Creating a design matrix for the mediator model */

data Temp/View=Temp; set Datasim_1;
   FakeY = 0;
run;
proc logistic data=Temp outdesign=RefDesign(drop=FakeY) outdesignonly;
   class x_i(ref="0") time(ref="1") / param=reference; 
   model FakeY = x_i time;
run;

proc iml;
use Refdesign;
read all;
X = Intercept || x_i1 || time2 || time3 || time4 || time5 || time6 ;
beta = {4, -1.7, 3.3, 3.8, 4.3, 4.8, 5.3};  
xb = X * beta;
close Refdesign;
create xb from xb[colname='xb']; 
append from xb;
quit;

/* Creating a design matrix for the outcome model */

proc iml;
use Refdesign;
read all;
X = Intercept || x_i1 || time2 || time3 || time4 || time5 || time6 ;
gamma = {-1, -0.6, 0.15, 0.25, 0.4, 0.55, 0.7};
yb = X * gamma;
close Refdesign;
create yb from yb[colname='yb']; 
append from yb;
quit;

/* Creating the simulated data set */

data xb; set xb; fakeid = _N_; run;
data yb; set yb; fakeid = _N_; run;
data errors; set errors; fakeid = _N_; run;
data datasim_1; set datasim_1; fakeid = _N_; run;

proc sort data=datasim_1; by fakeid time; run;
proc sort data=xb; by fakeid; run;
proc sort data=yb; by fakeid; run;
proc sort data=errors; by fakeid; run;

data datasim_2; 
merge  xb yb errors datasim_1;
by fakeid;
run;

%let gamma_m = 0.8;

data datasim; set datasim_2;
drop fakeid xb yb;
M_N_error = xb + bi_N + e_m;
M_SN_error = xb + bi_SN + e_m;
M_MN_error = xb + bi_MN + e_m;

M_N_true = xb + bi_N;
M_SN_true = xb + bi_SN;
M_MN_true = xb + bi_MN;

Y_N = yb + &gamma_m*M_N_true + bi_Y + e_Y;
Y_SN = yb + &gamma_m*M_SN_true + bi_Y + e_Y;
Y_MN = yb + &gamma_m*M_MN_true + bi_Y + e_Y;

if time = 2 then time2 = 1; else time2 = 0;
if time = 3 then time3 = 1; else time3 = 0;
if time = 4 then time4 = 1; else time4 = 0;
if time = 5 then time5 = 1; else time5 = 0;
if time = 6 then time6 = 1; else time6 = 0;
run;

proc datasets library=work nolist;
save datasim;
quit;
run; 

/*** Calibration model: Mediator with normal random effects (RE) ***/

proc mixed data=datasim method=reml;
class time(ref="1") x_i(ref="0");
model M_N_error = x_i time / solution outpred=pred_N_rand outpredm=pred_N_fixed;
random intercept / subject=id type=un;
run;

data M_N_calib_norm; set pred_N_rand; 	/* Get the calibrated mediator under normal assumption*/
keep id time pred;
rename pred = M_N_calib_norm;
run;

*** Normal RE: Computing the MSE between estimated RE from normal calibration and simulated normal RE   ***;

data fixed_rand_N; 
	set pred_N_rand (keep = id time pred bi_N bi_SN bi_MN);
	rename pred = fixed_rand_N;
run;
	
data fixed_N; 
	set pred_N_fixed (keep = id time pred);
	rename pred = fixed_N;
run;

data RE_N;
	merge fixed_rand_N fixed_N;
	by id time;
run;

data RE_N2;
	set RE_N;
	sqerr_NRE_N = (bi_N - (fixed_rand_N - fixed_N))**2;
run;

/*** Calibration model: Mediator with skewed normal RE ***/

proc mixed data=datasim method=reml;
class time(ref="1") x_i(ref="0");
model M_SN_error = x_i time / solution outpred=pred_SN_rand outpredm=pred_SN_fixed;
random intercept / subject=id type=un;
run;

data M_SN_calib_norm; set pred_SN_rand;
keep id time pred;
rename pred = M_SN_calib_norm;

*** SN RE: Computing the MSE between estimated RE from normal calibration and simulated normal RE   ***;

data fixed_rand_SN; 
	set pred_SN_rand(keep = id time pred bi_SN);
	rename pred = fixed_rand_SN;
run;
	
data fixed_SN; 
	set pred_SN_fixed (keep = id time pred);
	rename pred = fixed_SN;
run;

data RE_SN;
	merge fixed_rand_SN fixed_SN;
	by id time;
run;

data RE_SN2;
	set RE_SN;
	sqerr_SNRE_N = (bi_SN - 4 - (fixed_rand_SN - fixed_SN))**2; /* We subtract 4 (the mean), to centre the generated bi_SN at 0 for comparison with other square errors */
run;

/*** Calibration model: Mediator with mixture of normal RE ***/

proc mixed data=datasim method=reml;
class time(ref="1") x_i(ref="0");
model M_MN_error = x_i time / solution outpred=pred_MN_rand outpredm=pred_MN_fixed;
random intercept / subject=id type=un;
run;

data M_MN_calib_norm; set pred_MN_rand;
keep id time pred;
rename pred = M_MN_calib_norm;
run;

*** MN RE: Computng the MSE between estimated RE from normal calibration and simulated normal RE   ***;

data fixed_rand_MN; 
	set pred_MN_rand(keep = id time pred bi_MN);
	rename pred = fixed_rand_MN;
run;
	
data fixed_MN; 
	set pred_MN_fixed (keep = id time pred);
	rename pred = fixed_MN;
run;

data RE_MN;
	merge fixed_rand_MN fixed_MN;
	by id time;
run;

data RE_MN2;
	set RE_MN;
	sqerr_MNRE_N = (bi_MN - 1.4 - (fixed_rand_MN - fixed_MN))**2; /* We subtract 1.4 (the mean), to centre the generated bi_SN at 0 for comparison with other square errors */
run;

/*** Combining the calibrated mediators with the simulated data ***/

proc sort data=datasim; by id time; run;
proc sort data=M_N_calib_norm; by id time; run; 
proc sort data=M_SN_calib_norm; by id time; run; 
proc sort data=M_MN_calib_norm; by id time; run; 

data datasim_final; merge M_MN_calib_norm M_SN_calib_norm datasim M_N_calib_norm; by id time; run;

/* Prepare datasets for ACME and ADE computation: Code adapted from Bauer, though you could as well 
/* fit models separately and then extract the estimates as necessary 										*/

/* Bauer, D.J., Preacher, K.J. & Gil, K.M. (2006). Conceptualizing and testing random indirect effects     	*/
/* and moderated mediation in multilevel models: new procedures and recommendations.. Psychological Methods, 11, 142-163. 	*/

**** NORMAL : Mediator calibrated assuming normality ****;

data dt_N_norm; set datasim_final;
keep id Y_N x_i time M_N_calib_norm M_N_error;
rename Y_N = Y;
rename x_i = X;
rename M_N_error = M;
run;

data dt_N_norm; set dt_N_norm;     
 Z = Y;            *first assigning Z the value of Y;
 Sy = 1;           *setting Sy selection variable to 1 to indicate Z is Y;
 Sm = 0;           *setting Sm selection variable to 0 to indicate Z is not M;
 dv = 'Y';         *creating variable dv also to differentiate Y from M;
 output;           *outputting the first record for the observation;
 Z = M;            *now assigning Z the value of M;
 Sy = 0;           *setting Sy selection variable to 0 to indicate Z is not Y;
 Sm = 1;           *setting Sm selection variable to 1 to indicate Z is M;
 dv = 'M';         *creating variable dv also to differentiate Y from M;
 output;           *outputting the second record for the observation;
run;

/* Fitting the outcome and mediator models simultaneously */

proc mixed data=dt_N_norm noclprint method=reml;
 class dv time(ref="1");
 model Z = Sm Sm*X Sy Sy*M_N_calib_norm Sy*X dv*time /noint solution ddfm=kr; 
 repeated / group=dv subject=id;
 random intercept / group=dv subject=id;
 ods output solutionf=estfix_N_norm;
run;

data _null_; set Estfix_N_norm;
    if Effect="Sm*X" then call symput("a_N_norm", estimate);
    if Effect="Sy*M_N_calib_norm" then call symput("b_N_norm", estimate);
    if Effect="X*Sy" then call symput("ADE_N_norm", estimate);
	if Effect="X*Sy" then call symput("ADE_N_SE_norm", StdErr);
run;


**** CHISQUARE RE: Mediator calibrated assuming normality ****;

data dt_SN_norm; set datasim_final;
keep id Y_SN x_i time M_SN_calib_norm M_SN_error;
rename Y_SN = Y;
rename x_i = X;
rename M_SN_error = M;
run;

data dt_SN_norm; set dt_SN_norm;     
 Z = Y;            
 Sy = 1;           
 Sm = 0;           
 dv = 'Y';         
 output;           
 Z = M;            
 Sy = 0;           
 Sm = 1;           
 dv = 'M';         
 output;           
run;

proc mixed data=dt_SN_norm noclprint method=reml;
 class dv time(ref="1");
 model Z = Sm Sm*X Sy Sy*M_SN_calib_norm Sy*X dv*time /noint solution ddfm=kr;
 repeated / group=dv subject=id;
 random intercept / group=dv subject=id;
 ods output solutionf=estfix_SN_norm ;
run;

data _null_; set Estfix_SN_norm;
    if Effect="Sm*X" then call symput("a_SN_norm", estimate);
    if Effect="Sy*M_SN_calib_norm" then call symput("b_SN_norm", estimate);
    if Effect="X*Sy" then call symput("ADE_SN_norm", estimate);
	if Effect="X*Sy" then call symput("ADE_SN_SE_norm", StdErr);
run;

**** MIXTURE OF NORMALS RE: Mediator calibrated assuming normality ****;

data dt_MN_norm; set datasim_final;
keep id Y_MN x_i time M_MN_calib_norm M_MN_error;
rename Y_MN = Y;
rename x_i = X;
rename M_MN_error = M;
run;

data dt_MN_norm; set dt_MN_norm;     
 Z = Y;            
 Sy = 1;           
 Sm = 0;           
 dv = 'Y';         
 output;           
 Z = M;            
 Sy = 0;           
 Sm = 1;           
 dv = 'M';         
 output;           
run;

proc mixed data=dt_MN_norm noclprint method=reml;
 class dv time(ref="1");
 model Z = Sm Sm*X Sy Sy*M_MN_calib_norm Sy*X dv*time /noint solution ddfm=kr;
 repeated / group=dv subject=id;
 random intercept / group=dv subject=id;
 ods output solutionf=estfix_MN_norm;
run;

data _null_; set Estfix_MN_norm;
    if Effect="Sm*X" then call symput("a_MN_norm", estimate);
    if Effect="Sy*M_MN_calib_norm" then call symput("b_MN_norm", estimate);
    if Effect="X*Sy" then call symput("ADE_MN_norm", estimate);
	if Effect="X*Sy" then call symput("ADE_MN_SE_norm", StdErr);
run;

**** NORMAL : Mediator un-calibrated ****;

data dt_N_error; set datasim_final;
keep id Y_N x_i time M_N_error;
rename Y_N = Y;
rename x_i = X;
rename M_N_error = M;
run;

data dt_N_error; set dt_N_error;     
 Z = Y;            *first assigning Z the value of Y;
 Sy = 1;           *setting Sy selection variable to 1 to indicate Z is Y;
 Sm = 0;           *setting Sm selection variable to 0 to indicate Z is not M;
 dv = 'Y';         *creating variable dv also to differentiate Y from M;
 output;           *outputting the first record for the observation;
 Z = M;            *now assigning Z the value of M;
 Sy = 0;           *setting Sy selection variable to 0 to indicate Z is not Y;
 Sm = 1;           *setting Sm selection variable to 1 to indicate Z is M;
 dv = 'M';         *creating variable dv also to differentiate Y from M;
 output;           *outputting the second record for the observation;
run;

proc mixed data=dt_N_error noclprint method=reml ;
 class dv time(ref="1");
 model Z = Sm Sm*X Sy Sy*M Sy*X dv*time /noint solution ddfm=kr;
 repeated / group=dv subject=id;
 random intercept / group=dv subject=id;
 ods output solutionf=estfix_N_error;
run;

data _null_; set Estfix_N_error;
    if Effect="Sm*X" then call symput("a_N_error", estimate);
    if Effect="Sy*M" then call symput("b_N_error", estimate);
    if Effect="X*Sy" then call symput("ADE_N_error", estimate);
	if Effect="X*Sy" then call symput("ADE_N_SE_error", StdErr);
run;

**** CHISQUARE RE: Mediator uncalibrated ****;

data dt_SN_error; set datasim_final;
keep id Y_SN x_i time M_SN_error;
rename Y_SN = Y;
rename x_i = X;
rename M_SN_error = M;
run;

data dt_SN_error; set dt_SN_error;     
 Z = Y;            
 Sy = 1;           
 Sm = 0;           
 dv = 'Y';         
 output;           
 Z = M;            
 Sy = 0;           
 Sm = 1;           
 dv = 'M';         
 output;           
run;

proc mixed data=dt_SN_error noclprint method=reml ;
 class dv time(ref="1");
 model Z = Sm Sm*X Sy Sy*M Sy*X dv*time /noint solution ddfm=kr;
 repeated / group=dv subject=id;
 random intercept / group=dv subject=id;
 ods output solutionf=estfix_SN_error;
run;

data _null_; set Estfix_SN_error;
    if Effect="Sm*X" then call symput("a_SN_error", estimate);
    if Effect="Sy*M" then call symput("b_SN_error", estimate);
    if Effect="X*Sy" then call symput("ADE_SN_error", estimate);
	if Effect="X*Sy" then call symput("ADE_SN_SE_error", StdErr);
run;

**** MIXTURE OF NORMALS RE: Mediator uncalibrated  ****;

data dt_MN_error; set datasim_final;
keep id Y_MN x_i time M_MN_error;
rename Y_MN = Y;
rename x_i = X;
rename M_MN_error = M;
run;

data dt_MN_error; set dt_MN_error;     
 Z = Y;            
 Sy = 1;           
 Sm = 0;           
 dv = 'Y';         
 output;           
 Z = M;            
 Sy = 0;           
 Sm = 1;           
 dv = 'M';         
 output;           
run;

proc mixed data=dt_MN_error noclprint method=reml ;
 class dv time(ref="1");
 model Z = Sm Sm*X Sy Sy*M Sy*X dv*time/noint solution ddfm=kr;
 repeated / group=dv subject=id;
 random intercept / group=dv subject=id;
 ods output solutionf=estfix_MN_error;
run;

data _null_; set Estfix_MN_error;
    if Effect="Sm*X" then call symput("a_MN_error", estimate);
    if Effect="Sy*M" then call symput("b_MN_error", estimate);
    if Effect="X*Sy" then call symput("ADE_MN_error", estimate);
	if Effect="X*Sy" then call symput("ADE_MN_SE_error", StdErr);
run;

/* Computing and storing ADE and ACME estimates after assuming normal RE */

data estimates_N(keep=ACME_N_norm ADE_N_norm  
                   ACME_SN_norm ADE_SN_norm  
                   ACME_MN_norm ADE_MN_norm  

				   ACME_N_error ADE_N_error 
                   ACME_SN_error ADE_SN_error 
                   ACME_MN_error ADE_MN_error);

	ACME_N_norm = &a_N_norm*&b_N_norm;
	ADE_N_norm = &ADE_N_norm;

	ACME_SN_norm = &a_SN_norm*&b_SN_norm ;
	ADE_SN_norm = &ADE_SN_norm;

	ACME_MN_norm = &a_MN_norm*&b_MN_norm ;
	ADE_MN_norm = &ADE_MN_norm;

	ACME_N_error = &a_N_error*&b_N_error;
	ADE_N_error = &ADE_N_error;

	ACME_SN_error = &a_SN_error*&b_SN_error;
	ADE_SN_error = &ADE_SN_error;

	ACME_MN_error = &a_MN_error*&b_MN_error;
	ADE_MN_error = &ADE_MN_error;

run;

proc datasets library=work nolist;
save estimates_N Datasim RE_N2 RE_SN2 RE_MN2;
quit;
run;

proc append base=home.est_N_uncorr data = estimates_N; run; 

/* Next, we use the SNP calibration for the mediator */

data Temp2/view=Temp2; set Datasim;
   FakeY = 0;
run;

proc logistic data=Temp2 outdesign=RefDesign2(drop=FakeY) outdesignonly;
   class x_i(ref="0") time(ref="1") / param=reference; 
   model FakeY = x_i time;
run;

*filename slmm 'C:\..............\Simulation\slmm.mac'; /* Specify the directory where you've stored the SNP code */
%include slmm;

%let a_N_SNP=;
%let b_N_SNP=;
%let ADE_N_SNP=;


%let a_SN_SNP=;
%let b_SN_SNP=;
%let ADE_SN_SNP=;


%let a_MN_SNP=;
%let b_MN_SNP=;
%let ADE_MN_SNP=;

/*******************************************************************************************************/

%slmm(data=datasim, dep=M_N_error, xvar= x_i time2 time3 time4 time5 time6, kz=0, id=id, print=Y, 
method=gs, gridpt=4, outparm = params_N_0); /* If excluded, somehow the 'params_N_0' will not exist */

%slmm(data=datasim, dep=M_N_error, xvar= x_i time2 time3 time4 time5 time6, kz=0, id=id, print=Y, 
method=gs, gridpt=4, outparm = params_N_0);
data param_N_0; set params_N_0; keep aic k; k = 0; run;

%slmm(data=datasim, dep=M_N_error, xvar= x_i time2 time3 time4 time5 time6, kz=1, id=id, print=Y, 
method=gs, gridpt=4, outparm = params_N_1);
data param_N_1; set params_N_1; keep aic k; k = 1; run;

%slmm(data=datasim, dep=M_N_error, xvar= x_i time2 time3 time4 time5 time6, kz=2, id=id, print=Y, 
method=gs, gridpt=4, outparm = params_N_2);
data param_N_2; set params_N_2; keep aic k; k = 2; run;

/*******************************************************************************************************/

/* Choosing a k with the minimum AIC */

data param_N; set param_N_0 param_N_1 param_N_2; run;

proc sql; /* SQL code that picks out the record with the minimum AIC */
create table N_minim_k as
select * from param_N
having aic=min(aic);
quit;

data _null_; set N_minim_k;
	call symput("N_k", k);
run;

%slmm(data=datasim, dep=M_N_error, xvar= x_i time2 time3 time4 time5 time6, kz=&N_k, id=id, print=Y, 
method=gs, gridpt=4, outrand=rand_N, outparm = param_N_SNP)

* Extract fixed effects from SNP model, create design matrix, and add the SNP RE;

data _null_; set param_N_SNP;
	call symput("xi_N", X_I);
	call symput("TIME2_N", TIME2);
	call symput("TIME3_N", TIME3);  
    call symput("TIME4_N", TIME4);
	call symput("TIME5_N", TIME5);
	call symput("TIME6_N", TIME6);
run;

proc iml;
use Refdesign2;
read all;
X_design = x_i1 || time2 || time3 || time4 || time5 || time6 ;    		
params_N = {&xi_N, &TIME2_N, &TIME3_N, &TIME4_N, &TIME5_N, &TIME6_N};
mu_N = X_design * params_N;
create mu_N from mu_N[colname='mu_N']; 
append from mu_N;
quit;

data mu_N; set mu_N; fakeid = _N_; run;
data datasim; set datasim; fakeid = _N_; run;

proc sort data=datasim; by fakeid time; run;
proc sort data=mu_N; by fakeid; run;

data datasim_N; 
merge  mu_N datasim;
by fakeid;
run;

proc sort data=rand_N; by id; run;
proc sort data=datasim_N; by id; run;
data pred_N; merge rand_N datasim_N; by id; run;

data M_N_calib_SNP; set pred_N; /* Get fixed part of model from SNP + RE estimated by SNP */
keep id time M_N_calib_SNP;
M_N_calib_SNP = mu_N + effect;
run;

** Normal RE: Computng the MSE between estimated RE using SNP calibration and simulated normal RE **;

data rand_N1; set rand_N(keep=id effect);
	rename effect = SNP_NRE_SNP;
run;

data rand_normal;
	merge RE_N2 rand_N1;
	by id;
run;

data RE_N3; set rand_normal (keep= id bi_N SNP_NRE_SNP sqerr_NRE_N);
	sqerr_NRE_SNP = (bi_N + 4 - SNP_NRE_SNP)**2;   /* SNP approach estimates are without intercept, so we add it back */
run;

/*******************************************************************************************************/

%slmm(data=datasim, dep=M_SN_error, xvar= x_i time2 time3 time4 time5 time6, kz=0, id=id, print=Y, 
method=gs, gridpt=4, outparm = params_SN_0); 

%slmm(data=datasim, dep=M_SN_error, xvar= x_i time2 time3 time4 time5 time6, kz=1, id=id, print=Y, 
method=gs, gridpt=4, outparm = params_SN_1);

%slmm(data=datasim, dep=M_SN_error, xvar= x_i time2 time3 time4 time5 time6, kz=2, id=id, print=Y, 
method=gs, gridpt=4, outparm = params_SN_2);

/*******************************************************************************************************/

/* Choosing a k with the minimum AIC */

data params_SN_0; set params_SN_0; keep aic k; k = 0; run;
data params_SN_1; set params_SN_1; keep aic k; k = 1; run;
data params_SN_2; set params_SN_2; keep aic k; k = 2; run;

data params_SN; set params_SN_0 params_SN_1 params_SN_2; run; 

proc sql; /* SQL code that picks out the record with the minimum aic */
create table SN_minim_k as
select * from params_SN
having aic=min(aic);
quit;

data _null_; set SN_minim_k;
	call symput("SN_k", k);
run;

%slmm(data=datasim, dep=M_SN_error, xvar= x_i time2 time3 time4 time5 time6, kz=&SN_k, id=id, print=Y, 
method=gs, gridpt=4, outrand=rand_SN, outparm = param_SN_SNP);


* Extract fixed effects from SNP model, create design matrix, and add the SNP RE;

data _null_; set param_SN_SNP;
	call symput("xi_SN", X_I);
	call symput("TIME2_SN", TIME2);
	call symput("TIME3_SN", TIME3);  
	call symput("TIME4_SN", TIME4);
	call symput("TIME5_SN", TIME5);
	call symput("TIME6_SN", TIME6);
run;

proc iml;
use Refdesign2;
read all;
X_design = x_i1 || time2 || time3 || time4 || time5 || time6 ;
params_SN = {&xi_SN, &TIME2_SN, &TIME3_SN, &TIME4_SN, &TIME5_SN, &TIME6_SN};
mu_SN = X_design * params_SN;
create mu_SN from mu_SN[colname='mu_SN']; 
append from mu_SN;
quit;

data mu_SN; set mu_SN; fakeid = _N_; run;
data datasim; set datasim; fakeid = _N_; run;

proc sort data=datasim; by fakeid time; run;
proc sort data=mu_SN; by fakeid; run;

data datasim_SN; 
merge  mu_SN datasim;
by fakeid;
run;

proc sort data=rand_SN; by id; run;
proc sort data=datasim_SN; by id; run;
data pred_SN; merge rand_SN datasim_SN; by id; run;

data M_SN_calib_SNP; set pred_SN; /* Get fixed part of model from SNP + RE estimated by SNP*/
keep id time M_SN_calib_SNP;
M_SN_calib_SNP = mu_SN + effect;
run;

** SN RE: Computng the MSE between estimated RE using SNP calibration and simulated normal RE **;

data rand_SN1; set rand_SN(keep=id effect);
	rename effect = SNP_SNRE_SNP;
run;

data rand_skewnormal;
	merge RE_SN2 rand_SN1;
	by id;
run;

data RE_SN3; set rand_skewnormal(keep=id bi_SN SNP_SNRE_SNP sqerr_SNRE_N);
	sqerr_SNRE_SNP = (bi_SN + 4 - SNP_SNRE_SNP)**2;   /* SNP approach estimates are without intercept, so we add it back */
run;

/*******************************************************************************************************/

%slmm(data=datasim, dep=M_MN_error, xvar= x_i time2 time3 time4 time5 time6, kz=0, id=id, print=Y, 
method=gs, gridpt=4, outparm = params_MN_0); 

%slmm(data=datasim, dep=M_MN_error, xvar= x_i time2 time3 time4 time5 time6, kz=1, id=id, print=Y, 
method=gs, gridpt=4, outparm = params_MN_1);

%slmm(data=datasim, dep=M_MN_error, xvar= x_i time2 time3 time4 time5 time6, kz=2, id=id, print=Y, 
method=gs, gridpt=4, outparm = params_MN_2);

/*******************************************************************************************************/

/* Choosing a k with the minimum AIC */

data params_MN_0; set params_MN_0; keep aic k; k = 0; run;
data params_MN_1; set params_MN_1; keep aic k; k = 1; run;
data params_MN_2; set params_MN_2; keep aic k; k = 2; run;

data params_MN; set params_MN_0 params_MN_1 params_MN_2; run;

proc sql; /* SQL code that picks out the record with the minimum aic */
create table MN_minim_k as
select * from params_MN
having aic=min(aic);
quit;

data _null_; set MN_minim_k;
	call symput("MN_k", k);
run;

%slmm(data=datasim, dep=M_MN_error, xvar= x_i time2 time3 time4 time5 time6, kz=&MN_k, id=id, print=Y, 
method=gs, gridpt=4, outrand=rand_MN, outparm = param_MN_SNP);


* Extract fixed effects from SNP model, create design matrix, and add the SNP RE;

data _null_; set param_MN_SNP;
	call symput("xi_MN", X_I);
	call symput("TIME2_MN", TIME2);
	call symput("TIME3_MN", TIME3);  
	call symput("TIME4_MN", TIME4);
	call symput("TIME5_MN", TIME5);
	call symput("TIME6_MN", TIME6);
run;


proc iml;
use Refdesign2;
read all;
X_design = x_i1 || time2 || time3 || time4 || time5 || time6 ;
params_MN = {&xi_MN, &TIME2_MN, &TIME3_MN, &TIME4_MN, &TIME5_MN, &TIME6_MN};
mu_MN = X_design * params_MN;
create mu_MN from mu_MN[colname='mu_MN']; 
append from mu_MN;
quit;

data mu_MN; set mu_MN; fakeid = _N_; run;
data datasim; set datasim; fakeid = _N_; run;

proc sort data=datasim; by fakeid time; run;
proc sort data=mu_MN; by fakeid; run;

data datasim_MN; 
	merge  mu_MN datasim;
	by fakeid;
run;

proc sort data=rand_MN; by id; run;
proc sort data=datasim_MN; by id; run;
data pred_MN; merge rand_MN datasim_MN; by id; run;

data M_MN_calib_SNP; set pred_MN;     /* Get fixed part of model + RE estimated by SNP*/
keep id time M_MN_calib_SNP;
M_MN_calib_SNP = mu_MN + effect;
run;

data datasim_SNP;
	merge datasim M_MN_calib_SNP M_N_calib_SNP M_SN_calib_SNP;
	by id time;
run;

** MN RE: Computng the MSE between estimated RE using SNP calibration and simulated normal RE **;

data rand_MN1; set rand_MN(keep=id effect);
	rename effect = SNP_MNRE_SNP;
run;

data rand_MN2;
	merge RE_MN2 rand_MN1;
	by id;
run;

data RE_MN4; set rand_MN2(keep=id bi_MN SNP_MNRE_SNP sqerr_MNRE_N);
	sqerr_MNRE_SNP = (bi_MN + 4 - SNP_MNRE_SNP)**2;  /* SNP approach estimates are without intercept, so we add it back */
run;

******* Storing the MSE (comparing simulated RE to those predicted from Normal and SNP models)*********;

data est_re;
	merge RE_N3 RE_SN3 RE_MN4;
	by id;
run;

proc means data=est_re; 
	var sqerr_NRE_N sqerr_NRE_SNP sqerr_SNRE_N sqerr_SNRE_SNP sqerr_MNRE_N sqerr_MNRE_SNP;
	output out=out_re_mse mean=; 
run;

data est_re_mse2; 
	set out_re_mse(drop=_TYPE_ _FREQ_);
run;

proc append base=home.mse_re data = est_re_mse2; run;  /* Appending MSE results */

**** NORMAL RE: Mediator calibrated using SNP ****;

data dt_N_SNP; set datasim_SNP;
keep id Y_N x_i time M_N_calib_SNP M_N_error;
rename Y_N = Y;
rename x_i = X;
rename M_N_error = M;
run;

data dt_N_SNP; set dt_N_SNP;     
 Z = Y;            *first assigning Z the value of Y;
 Sy = 1;           *setting Sy selection variable to 1 to indicate Z is Y;
 Sm = 0;           *setting Sm selection variable to 0 to indicate Z is not M;
 dv = 'Y';         *creating variable dv also to differentiate Y from M;
 output;           *outputting the first record for the observation;
 Z = M;            *now assigning Z the value of M;
 Sy = 0;           *setting Sy selection variable to 0 to indicate Z is not Y;
 Sm = 1;           *setting Sm selection variable to 1 to indicate Z is M;
 dv = 'M';         *creating variable dv also to differentiate Y from M;
 output;           *outputting the second record for the observation;
run;

proc mixed data=dt_N_SNP noclprint method=reml;
 class dv time(ref="1");
 model Z = Sm Sm*X Sy Sy*M_N_calib_SNP Sy*X dv*time /noint solution ddfm=kr;
 repeated / group=dv subject=id;
 random intercept / group=dv subject=id;
 ods output covb=acovfix_N_SNP solutionf=estfix_N_SNP;
run;

data _null_; set Estfix_N_SNP;
    if Effect="Sm*X" then call symput("a_N_SNP", estimate);
    if Effect="Sy*M_N_calib_SNP" then call symput("b_N_SNP", estimate);
    if Effect="X*Sy" then call symput("ADE_N_SNP", estimate);
	if Effect="X*Sy" then call symput("ADE_N_SE_SNP", StdErr);
run;

**** CHISQUARE RE: Mediator calibrted using SNP ****;

data dt_SN_SNP; set datasim_SNP;
keep id Y_SN x_i time M_SN_calib_SNP M_SN_error;
rename Y_SN = Y;
rename x_i = X;
rename M_SN_error = M;
run;

data dt_SN_SNP; set dt_SN_SNP;     
 Z = Y;            
 Sy = 1;           
 Sm = 0;           
 dv = 'Y';         
 output;           
 Z = M;            
 Sy = 0;           
 Sm = 1;           
 dv = 'M';         
 output;           
run;

proc mixed data=dt_SN_SNP noclprint method=reml;
 class dv time(ref="1");
 model Z = Sm Sm*X Sy Sy*M_SN_calib_SNP Sy*X dv*time /noint solution ddfm=kr;
 repeated / group=dv subject=id;
 random intercept / group=dv subject=id;
 ods output solutionf=estfix_SN_SNP;
run;

data _null_; set Estfix_SN_SNP;
    if Effect="Sm*X" then call symput("a_SN_SNP", estimate);
    if Effect="Sy*M_SN_calib_SNP" then call symput("b_SN_SNP", estimate);
    if Effect="X*Sy" then call symput("ADE_SN_SNP", estimate);
	if Effect="X*Sy" then call symput("ADE_SN_SE_SNP", StdErr);
run;

**** MIXTURE OF NORMALS : Mediator calibrated using SNP****;

data dt_MN_SNP; set datasim_SNP;
keep id Y_MN x_i time M_MN_calib_SNP M_MN_error;
rename Y_MN = Y;
rename x_i = X;
rename M_MN_error = M;
run;

data dt_MN_SNP; set dt_MN_SNP;     
 Z = Y;            
 Sy = 1;           
 Sm = 0;           
 dv = 'Y';         
 output;           
 Z = M;            
 Sy = 0;           
 Sm = 1;           
 dv = 'M';         
 output;           
run;

proc mixed data=dt_MN_SNP noclprint method=reml;
 class dv time(ref="1");
 model Z = Sm Sm*X Sy Sy*M_MN_calib_SNP Sy*X dv*time /noint solution ddfm=kr;
 repeated / group=dv subject=id;
 random intercept / group=dv subject=id;
 ods output solutionf=estfix_MN_SNP;
run;


data _null_; set Estfix_MN_SNP;
    if Effect="Sm*X" then call symput("a_MN_SNP", estimate);
    if Effect="Sy*M_MN_calib_SNP" then call symput("b_MN_SNP", estimate);
    if Effect="X*Sy" then call symput("ADE_MN_SNP", estimate);
	if Effect="X*Sy" then call symput("ADE_MN_SE_SNP", StdErr);
run;

/* Computing and storing ADE and ACME estimates after assuming SNP aaproach for the RE */

data estimates_SNP(keep= ACME_N_SNP ADE_N_SNP 
                     ACME_SN_SNP ADE_SN_SNP 
                     ACME_MN_SNP ADE_MN_SNP );

	ACME_N_SNP = &a_N_SNP*&b_N_SNP;
	ADE_N_SNP = &ADE_N_SNP;

	ACME_SN_SNP = &a_SN_SNP*&b_SN_SNP ;
	ADE_SN_SNP = &ADE_SN_SNP;

	ACME_MN_SNP = &a_MN_SNP*&b_MN_SNP;
	ADE_MN_SNP = &ADE_MN_SNP;

run;

proc datasets library=work nolist;
save estimates_SNP;
quit;
run;

proc append base=home.est_SNP data = estimates_SNP; run; 

%end;
%mend;


%Simulate(2000);

* Analyzing the results;
* True ACME = -1.360 and true ADE = -0.6;

ods graphics on;
ods exclude none;
ods results;

data estimates_n; set home.est_N_uncorr;
bias_ACME_N_norm = ACME_N_norm - (-1.36);
bias_ACME_SN_norm = ACME_SN_norm - (-1.36);
bias_ACME_MN_norm = ACME_MN_norm - (-1.36);

bias_ACME_N_error = ACME_N_error - (-1.36);
bias_ACME_SN_error = ACME_SN_error - (-1.36);
bias_ACME_MN_error = ACME_MN_error - (-1.36);

bias_ADE_N_norm = ADE_N_norm - (-0.6);
bias_ADE_SN_norm = ADE_SN_norm - (-0.6);
bias_ADE_MN_norm = ADE_MN_norm - (-0.6);

bias_ADE_N_error= ADE_N_error - (-0.6);
bias_ADE_SN_error = ADE_SN_error - (-0.6);
bias_ADE_MN_error = ADE_MN_error - (-0.6);

SE_ACME_N_norm = bias_ACME_N_norm**2;
SE_ACME_SN_norm = bias_ACME_SN_norm**2; 
SE_ACME_MN_norm = bias_ACME_MN_norm**2;

SE_ACME_N_error = bias_ACME_N_error**2;
SE_ACME_SN_error = bias_ACME_SN_error**2; 
SE_ACME_MN_error = bias_ACME_MN_error**2; 

SE_ADE_N_norm = bias_ADE_N_norm**2; 
SE_ADE_SN_norm = bias_ADE_SN_norm**2; 
SE_ADE_MN_norm = bias_ADE_MN_norm**2; 

SE_ADE_N_error = bias_ADE_N_error**2; 
SE_ADE_SN_error = bias_ADE_SN_error**2; 
SE_ADE_MN_error = bias_ADE_MN_error**2;

run;
 
data estimates_snp; set home.est_SNP;

bias_ACME_N_SNP = ACME_N_SNP - (-1.36);
bias_ACME_SN_SNP = ACME_SN_SNP - (-1.36);
bias_ACME_MN_SNP = ACME_MN_SNP - (-1.36);

bias_ADE_N_SNP = ADE_N_SNP - (-0.6);
bias_ADE_SN_SNP = ADE_SN_SNP - (-0.6);
bias_ADE_MN_SNP = ADE_MN_SNP - (-0.6);

SE_ACME_N_SNP = bias_ACME_N_SNP**2; 
SE_ACME_SN_SNP = bias_ACME_SN_SNP**2;
SE_ACME_MN_SNP = bias_ACME_MN_SNP**2; 

SE_ADE_N_SNP = bias_ADE_N_SNP**2; 
SE_ADE_SN_SNP = bias_ADE_SN_SNP**2; 
SE_ADE_MN_SNP = bias_ADE_MN_SNP**2;

run;


proc means data = estimates_n; run;
proc means data = estimates_snp; run;

