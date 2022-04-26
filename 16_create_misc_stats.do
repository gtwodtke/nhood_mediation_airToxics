#delimit cr
capture clear all 
capture log close 
set more off 

global projdir "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\" 
log using "${projdir}programs\_LOGS\16_create_misc_stats.log", replace 

/******************************************************************
PROGRAM NAME: 16_create_misc_stats.do
AUTHOR: GW, BP
PURPOSE: create misc stats reported in text
*******************************************************************/

/********
LOAD DATA
*********/
use "${projdir}data\eclsb\v04_eclsb_mi.dta"
keep if minum==0

/***********
CREATE STATS
************/

/*PROPORTION OF MISSING INFO*/
/**Note: variables with no missing values are excluded from the numerator**/

/*numerator*/
global allvars /*gender*/ race /*twinid*/ birthwt momage01 dadage01 wic01 ///
	foodst01 medicd01 tanf01 ///
	nhpovrt01 /*faminc01*/ pared01 parocc01 momemp01 dademp01 /*hhtotal01*/ ///
	/*biodad01*/ married01 house01 rbooks01 /*prmlang01*/ age01 age05 region01 ///
	/*urban01*/ nhpopden01 ///
	co_2001 no2_2001 o3_2001 pm10_2001 pm25_2001 so2_2001 chem6_2001 ///
	chem10_2001 chem12_2001 chem17_2001 chem36_2001 chem39_2001 chem40_2001 ///
	chem70_2001 chem75_2001 chem108_2001 chem109_2001 chem114_2001 ///
	chem141_2001 chem163_2001 chem180_2001 chem217_2001 chem225_2001 ///
	chem236_2001 chem290_2001 chem293_2001 chem323_2001 chem332_2001 ///
	chem346_2001 chem347_2001 chem351_2001 chem355_2001 chem356_2001 ///
	chem359_2001 chem360_2001 chem362_2001 chem364_2001 chem381_2001 ///
	chem447_2001 chem454_2001 chem474_2001 chem519_2001 chem522_2001 ///
	chem531_2001 chem549_2001 chem565_2001 chem566_2001 chem567_2001 ///
	chem572_2001 chem586_2001 chem592_2001 chem599_2001 ///
	maththeta05 readtheta05
local miss_count 0
foreach x in $allvars {
	misstable summarize `x'
	local miss_count = `miss_count' + r(N_eq_dot) + r(N_gt_dot)
	scalar miss_total = `miss_count'
}
di miss_total

/*denominator*/
describe
local denom = 80*`r(N)' 
scalar denom = `denom'
di "`denom'"

/*proportion of missing info*/
di miss_total/denom

/*SUMMARY OF COVARIATES W/O IMPUTATION*/
tab1 gender race region01 urban01 married01 biodad01 house01 ///
	prmlang01 rbooks01 dademp01 momemp01 pared01 ///
	wic01 foodst01 medicd01 tanf01 [iw=sampwt]
sum birthwt age01 age05 momage01 dadage01 faminc01 parocc01 hhtotal01 ///
	nhpovrt01 maththeta05 readtheta05 [iw=sampwt]
	
/*SCALING SENSITIVITY PARAMETERS*/
use "${projdir}data\eclsb\v04_eclsb_mi.dta", clear
keep if minum!=0

egen maththeta05_std = std(maththeta05)
egen readtheta05_std = std(readtheta05)

egen faminc01_std = std(faminc01)
egen parocc01_std = std(parocc01)
egen momage01_std = std(momage01)

reg maththeta05_std faminc01_std parocc01_std momage01_std i.gender i.race i.twinid birthwt ///
	dadage01 i.wic01 i.foodst01 i.medicd01 i.tanf01 ///
	nhpovrt01 i.pared01 i.momemp01 i.dademp01 hhtotal01 i.biodad01 ///
	i.married01 i.house01 rbooks01 i.prmlang01 age01 age05 i.region01 i.urban01 nhpopden01 ///
	[iw=fnlwt]
	
reg readtheta05_std faminc01_std parocc01_std momage01_std i.gender i.race i.twinid birthwt ///
	dadage01 i.wic01 i.foodst01 i.medicd01 i.tanf01 ///
	nhpovrt01 i.pared01 i.momemp01 i.dademp01 hhtotal01 i.biodad01 ///
	i.married01 i.house01 rbooks01 i.prmlang01 age01 age05 i.region01 i.urban01 nhpopden01 ///
	[iw=fnlwt]
	
reg faminc01_std parocc01_std momage01_std i.gender i.race i.twinid birthwt ///
	dadage01 i.wic01 i.foodst01 i.medicd01 i.tanf01 ///
	nhpovrt01 pared01 i.momemp01 i.dademp01 hhtotal01 i.biodad01 ///
	i.married01 i.house01 rbooks01 i.prmlang01 age01 age05 i.region01 i.urban01 nhpopden01 ///
	[iw=fnlwt]	
lincom _b[nhpovrt01]*(0.25-0.05)

reg parocc01_std faminc01_std momage01_std i.gender i.race i.twinid birthwt ///
	dadage01 i.wic01 i.foodst01 i.medicd01 i.tanf01 ///
	nhpovrt01 pared01 i.momemp01 i.dademp01 hhtotal01 i.biodad01 ///
	i.married01 i.house01 rbooks01 i.prmlang01 age01 age05 i.region01 i.urban01 nhpopden01 ///
	[iw=fnlwt]	
lincom _b[nhpovrt01]*(0.25-0.05)

reg momage01_std parocc01_std faminc01_std i.gender i.race i.twinid birthwt ///
	dadage01 i.wic01 i.foodst01 i.medicd01 i.tanf01 ///
	nhpovrt01 pared01 i.momemp01 i.dademp01 hhtotal01 i.biodad01 ///
	i.married01 i.house01 rbooks01 i.prmlang01 age01 age05 i.region01 i.urban01 nhpopden01 ///
	[iw=fnlwt]	
lincom _b[nhpovrt01]*(0.25-0.05)

log close
