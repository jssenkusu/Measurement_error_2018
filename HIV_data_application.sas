
/*	Purpose: Apply mediation method to HIV-LIVE data	*/
/*	By: John Ssenkusu 								    */
/*  Date last modified: February 13, 2017		    	*/

* Notes ;
* 1. This data is available on request: https://www.bumc.bu.edu/care/research-studies/past-research-studies/hiv-live-study/ ;
* 2. Remember to change the libname to point to where you want your results stored ;
* 2. Modify the path of 'slmm.mac' to point to the directory where you've stored the SNP code ;

libname home 'C:\Data application';

PROC IMPORT OUT= WORK.HIV 
            DATAFILE= "C:\HIV_LIVE.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data hiv; set hiv; 
where pct3d_p ge 0 or pct3d_p = .;
	sqrt_cd4 = sqrt(cd4);
run;

proc sort data=HIV;
	by id tp;
run;

data hiv2; set hiv;
	where tp = 1;
run;

data hiv3; set hiv2 (keep=id race realm income niaahaz2 age);
	rename race = race2;
	rename realm = realm2;
	rename income = income2;
	rename niaahaz2 = niaahaz_BL;
	rename age = age_BL;
run;

data hiv_live; merge hiv hiv3; by id; run;

/* Calibration model: Mediator with normal RE */

/* Variables considered to include in model: realm, age_BL, homeless, hiv_qol */
/* Interraction terms between time and niaahaz_BL were not significant, thus excluded */

proc mixed data=hiv_live method=reml;
class TP(ref="1") race2(ref="1");
model pct3d_p = niaahaz_BL age_BL hiv_qol realm2 TP TP*niaahaz_BL/ solution outpred=pred_N_rand outpredm=pred_N_fixed;
random intercept / subject=id type=un ;
ods output solutionf=Est_adh_N_model(keep=effect estimate StdErr);
run;

data _null_; set Est_adh_N_model;
    if Effect="niaahaz_BL" then call symput("a", estimate);
	if Effect="niaahaz_BL" then call symput("a_SE", StdErr);
run;

data adh_Norm_calib; set pred_N_rand; 
keep id tp pred;
if pct3d_p = . then pred = .;
rename pred = adh_Norm_calib;
run;

data hiv_live; set hiv_live;
	where pct3d_p ne . and realm2 ne . and age_BL ne . and niaahaz_BL ne . and hiv_qol ne .;
	if tp = 2 then tp2 = 1; else tp2 = 0;
	if tp = 3 then tp3 = 1; else tp3 = 0;
	if tp = 4 then tp4 = 1; else tp4 = 0;
	if tp = 5 then tp5 = 1; else tp5 = 0;
	if tp = 6 then tp6 = 1; else tp6 = 0;
	if tp = 7 then tp7 = 1; else tp7 = 0;
	if tp = 8 then tp8 = 1; else tp8 = 0;

	if race = 2 then race_2 = 1; else race_2 = 0;
	if race = 3 then race_3 = 1; else race_3 = 0;
	if race = 4 then race_4 = 1; else race_4 = 0;

	tp2_BL_alcoh = tp2 * niaahaz_BL;
	tp3_BL_alcoh = tp3 * niaahaz_BL;
	tp4_BL_alcoh = tp4 * niaahaz_BL;
	tp5_BL_alcoh = tp5 * niaahaz_BL;
	tp6_BL_alcoh = tp6 * niaahaz_BL;
	tp7_BL_alcoh = tp7 * niaahaz_BL;
	tp8_BL_alcoh = tp8 * niaahaz_BL;
run;

filename slmm 'C:\slmm.mac';
%include slmm;

* Notes: Used AIC to determine the appropriate value for kz for each 'heavy alcohol' use group;
*		 The lowest for all: AIC, BIC and HQ is when kz = 2 ;

%slmm(data=hiv_live, dep=pct3d_p, xvar=niaahaz_BL age_BL hiv_qol realm2 tp2 tp3 tp4 tp5 tp6 tp7 tp8
tp2_BL_alcoh tp3_BL_alcoh tp4_BL_alcoh tp5_BL_alcoh tp6_BL_alcoh tp7_BL_alcoh tp8_BL_alcoh, 
kz=2, id=id, print=Y, method=gs, gridpt=4, outrand=rand_SNP, outparm = param_SNP);


* Storing fixed effects from SNP approach;

data _null_; set param_snp;
	call symput("TP2", TP2);
	call symput("TP3", TP3);
	call symput("TP4", TP4);
	call symput("TP5", TP5);
	call symput("TP6", TP6);
	call symput("TP7", TP7);
	call symput("TP8", TP8); 

	call symput("TP2_BL_alcoh", TP2_BL_ALCOH);
	call symput("TP3_BL_alcoh", TP3_BL_ALCOH);
	call symput("TP4_BL_alcoh", TP4_BL_ALCOH);
	call symput("TP5_BL_alcoh", TP5_BL_ALCOH);
	call symput("TP6_BL_alcoh", TP6_BL_ALCOH);
	call symput("TP7_BL_alcoh", TP7_BL_ALCOH);
	call symput("TP8_BL_alcoh", TP8_BL_ALCOH);
	
	call symput("niaahaz_BL", niaahaz_bl); 
	call symput("REALM2", realm2); 
	call symput("AGE_BL", age_BL); 
	call symput("HIV_QOL", hiv_qol); 
run;

* Creating a design matrix;

data Temp / view=Temp;
   set hiv_live;
   FakeY = 0;
run;

proc logistic data=Temp outdesign=Design(drop=FakeY) outdesignonly;
   class niaahaz_BL(ref="0") TP(ref="1") / param=reference; 
   model FakeY = tp tp*niaahaz_BL realm2 age_BL hiv_qol niaahaz_BL;
run;

proc iml;
use Design;
read all;
X_design = tp2 || tp3 || tp4 || tp5 || tp6 || tp7 || tp8 || niaahaz_BL1TP2 || niaahaz_BL1TP3 || niaahaz_BL1TP4 ||
niaahaz_BL1TP5 || niaahaz_BL1TP6 || niaahaz_BL1TP7 || niaahaz_BL1TP8 || realm2 || age_BL || hiv_qol || niaahaz_BL1;    
params_SNP = {&TP2, &TP3, &TP4, &TP5, &TP6, &TP7, &TP8, &TP2_BL_alcoh, &TP3_BL_alcoh, &TP4_BL_alcoh, &TP5_BL_alcoh, 
&TP6_BL_alcoh, &TP7_BL_alcoh, &TP8_BL_alcoh, &realm2, &age_BL, &hiv_qol, &niaahaz_BL};
mu_SNP = X_design * params_SNP;
create mu_SNP from mu_SNP[colname='mu_SNP']; 
append from mu_SNP;
quit;

* Calibrated adherence for patients using SNP;

data mu_SNP; set mu_SNP; fakeid = _N_; run;
data hiv_live; set hiv_live; fakeid = _N_; run;
proc sort data=hiv_live; by fakeid; run;
proc sort data=mu_SNP; by fakeid; run;
data dt_SNP; merge  mu_SNP hiv_live; by fakeid; run;
proc sort data=rand_snp; by id; run;
proc sort data=dt_SNP; by id; run;
data pred; merge rand_snp dt_SNP; by id; run;

data adh_SNP_calib; set pred; /* Get fixed part of model from SNP + RE estimated by SNP*/
keep id tp adh_SNP_calib;
adh_SNP_calib = mu_SNP + effect;
run;


proc sort data=adh_SNP_calib; by id tp; run;
proc sort data=adh_Norm_calib; by id tp; run;
proc sort data=hiv_live; by id tp; run;

data hiv_live_final;
merge hiv_live adh_SNP_calib  adh_Norm_calib ; 
by id tp; 
drop fakeid;
run;

proc datasets library=work nolist;
save hiv_live_final;
quit;
run;

data home.hiv_live_final; set hiv_live_final; run;

* Outcome model: Mediator (adherence) with error ;

proc mixed data=hiv_live_final method=reml;
class TP(ref="1");
model sqrt_cd4 = niaahaz_BL pct3d_p age_BL homeless hiv_qol TP / solution;
random intercept / subject=id type=un ;
ods output solutionf=Est_adh_error;
run;

data _null_; set Est_adh_error;
    if Effect="pct3d_p" then call symput("b_error", estimate);
	if Effect="pct3d_p" then call symput("b_SE_error", StdErr);
    if Effect="niaahaz_BL" then call symput("ADE_error", estimate);
	if Effect="niaahaz_BL" then call symput("ADE_SE_error", StdErr);
run;

* Outcome model: Mediator (adherence) calibrated using SNP ;

proc mixed data=hiv_live_final method=reml;
class TP(ref="1");
model sqrt_cd4 = niaahaz_BL adh_SNP_calib age_BL homeless hiv_qol TP / solution;
random intercept / subject=id type=un ;
ods output solutionf=Est_adh_snp(keep=effect estimate StdErr);
run;

data _null_; set Est_adh_snp;
    if Effect="adh_SNP_calib" then call symput("b_snp", estimate);
	if Effect="adh_SNP_calib" then call symput("b_SE_snp", StdErr);
    if Effect="niaahaz_BL" then call symput("ADE_snp", estimate);
	if Effect="niaahaz_BL" then call symput("ADE_SE_snp", StdErr);
run;

* Outcome model: Mediator (adherence) calibrated assuming normal RE;

proc mixed data=hiv_live_final method=reml;
class TP(ref="1");
model sqrt_cd4 = niaahaz_BL adh_Norm_calib age_BL homeless hiv_qol TP / solution;
random intercept / subject=id type=un ;
ods output solutionf=Est_adh_norm(keep=effect estimate StdErr);
run;

data _null_; set Est_adh_norm;
    if Effect="adh_Norm_calib" then call symput("b_norm", estimate);
	if Effect="adh_Norm_calib" then call symput("b_SE_norm", StdErr);
    if Effect="niaahaz_BL" then call symput("ADE_norm", estimate);
	if Effect="niaahaz_BL" then call symput("ADE_SE_norm", StdErr);
run;

* Storing estimates ;

data home.data_hiv_estimates(keep= ACME_error ADE_error  
							  	   ADE_norm ACME_norm 
							       ACME_SNP ADE_SNP );
	ACME_norm = &a*&b_norm;
	ADE_norm = &ADE_norm;

	ACME_error = &a*&b_error;
	ADE_error = &ADE_error;

	ACME_SNP = &a*&b_SNP;
	ADE_SNP = &ADE_SNP;

run;


