#delimit cr
capture clear all 
capture log close 
set more off 

global projdir "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\" 
log using "${projdir}programs\_LOGS\17_create_table_S1.log", replace 

/******************************************************************
PROGRAM NAME: 17_create_table_S1.do
AUTHOR: GW, BP
PURPOSE: create table of toxins by nhpoverty quartile - wave 1
*******************************************************************/

/********
LOAD DATA
*********/
use "${projdir}data\eclsb\v04_eclsb_mi.dta"

/******************
CREATE MI META VARS
*******************/
rename minum _mj
sort _mj caseid
by _mj: gen _mi=_n

/************
CREATE MACROS
*************/
global chems ///
	co_2001 no2_2001 o3_2001 pm10_2001 pm25_2001 so2_2001 ///
	chem6_2001 chem10_2001 chem12_2001 chem17_2001 chem36_2001 chem39_2001 ///
	chem40_2001 chem70_2001 chem75_2001 chem108_2001 chem109_2001 chem114_2001 ///
	chem141_2001 chem163_2001 chem180_2001 chem217_2001 chem225_2001 ///
	chem236_2001 chem290_2001 chem293_2001 chem323_2001 chem332_2001 /// 
	chem346_2001 chem347_2001 chem351_2001 chem355_2001 chem356_2001 chem359_2001 /// 
	chem360_2001 chem362_2001 chem364_2001 chem381_2001 chem447_2001 /// 
	chem454_2001 chem474_2001 chem519_2001 chem522_2001 chem531_2001 /// 
	chem549_2001 chem565_2001 chem566_2001 chem567_2001 chem572_2001 /// 
	chem586_2001 chem592_2001 chem599_2001 

/***************
COMPUTE QUANTILES
****************/
quietly qreg nhpovrt01 [iw=fnlwt] if _mj!=0, q(0.25)
scalar p25=_b[_cons]
quietly qreg nhpovrt01 [iw=fnlwt] if _mj!=0, q(0.50)
scalar p50=_b[_cons]
quietly qreg nhpovrt01 [iw=fnlwt] if _mj!=0, q(0.75)
scalar p75=_b[_cons]
gen nhpovqrt=.
replace nhpovqrt=1 if nhpovrt01<=p25
replace nhpovqrt=2 if nhpovrt01<=p50 & nhpovrt01>p25
replace nhpovqrt=3 if nhpovrt01<=p75 & nhpovrt01>p50
replace nhpovqrt=4 if nhpovrt01>p75 
tab nhpovqrt if _mj!=0, missing
tab nhpovqrt [iw=fnlwt] if _mj!=0, missing

/**********************
CREATE OUTPUT FOR TABLE
***********************/
/*OVERALL*/
sum *_2001 [iw=fnlwt] if _mj!=0

/*BY NH POVERTY QUARTILE*/
mi import ice, auto clear

gen clust=(10*strat)+psu

foreach v in $chems {
	di " "
	di " "
	di "`v'"
	quietly mi estimate, post: reg `v' i.nhpovqrt [pw=sampwt], cluster(clust)
	quietly predict yhat
	tabstat yhat, stats(mean) by(nhpovqrt)
	quietly drop yhat
	mi test 2.nhpovqrt 3.nhpovqrt 4.nhpovqrt
	}

log close
