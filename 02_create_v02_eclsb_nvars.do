#delimit cr
capture log close
capture clear all 

global data_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\data\"
global log_directory "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\programs\_LOGS\"
log using "${log_directory}02_create_v01_eclsb_nvars.log", replace 

/*******************************************************************
PROGRAM NAME: 02_create_v02_eclsb_nvars.do
AUTHOR: KW, GW
PURPOSE: clean raw data from ECLS-B and create new vars for analysis
********************************************************************/

use "${data_directory}eclsb\v01_eclsb.dta"

/***************
SUBSET VARIABLES
****************/
keep I_ID I_TWINID X4CHSEX Y1CHRACE W1R0 W1C0 W1F0 W1FC0 W2R0 W2C0 W2F0 W2FC0 ///
W2C2J0 W2C2P0 W22J0 W22P0 W22F0 W2C1F0 W3R0 W3C0 W3D0 W31R0 W31C0 W33D0 W312F0 ///
W33J0 W32J0 W33P0 W323P0 W4R0 W41R0 W4C0 W44J0 W44T0 W44JT0 W44S0 W4123D0 W4234J0 ///
W42P0 W43P0 W423P0 W5R0 W54R0 W545J0 W523J0 W55T0 W545T0 W55S0 WKR0 WK1R0 WK1C0 ///
WK45J0 WK23J0 WK45T0 WK45S0 WK123D0 WK12F0 WK2P0 WK3P0 WK23P0 X4ELIGBL XKSTATUS ///
XKWHENK XKKDATA X1CHPREM X1MBRTST X1ASMTMM X1ASMTYY X2ASMTMM X2ASMTYY X3ASMTMM ///
X3ASMTYY X4ASMTMM X4ASMTYY X5ASMTMM X5ASMTYY X3RTHR2 X3MTHR2 X4RTHR2 X4MTHR2 X5RTHR2 ///
X5MTHR2 X1RMTLS X2MTLSCL X1RMTRS X2MTRSCL X3FMBLCK X3FMFORM X4FMFORM X5FMFORM  ///
X2TASCLS X4NCESID X4CCDPSS X5GRDLVL X5RPTR X5NCESID X5CCDPSS ///
BCBRTHWT Y1HMRACE Y1MOMED Y1HMEMP Y1MOMSCR Y2MOMED Y2HMEMP Y2MOMSCR Y3MOMED Y3HMEMP ///
Y3MOMSCR X4MOMED X4HMEMP X4MOMSCR X5MOMED X5HMEMP X5MOMSCR Y1FTHTYP Y1HFRACE Y1FTHED ///
Y1HFEMP Y1FTHSCR Y2FTHTYP Y2FTHED Y2HFEMP Y2FTHSCR Y3FTHTYP Y3FTHED Y3HFEMP Y3FTHSCR ///
X4FTHTYP X4FTHED X4HFEMP X4FTHSCR X5FTHTYP X5FTHED X5HFEMP X5FTHSCR X1HMAG_B X1HFAG_B ///
X1MARSTA Y2MARSTA Y3MARSTA X4MARSTA X5MARSTA Y1NUMSIB Y2NUMSIB Y3NUMSIB X4NUMSIB /// 
X5NUMSIB Y1HTOTAL Y2HTOTAL Y3HTOTAL X4HTOTAL X5HTOTAL X1LANGST X1INCOME X2INCOME ///
X3INCOME X4INCOME X5INCOME X1HHREGN X2HHREGN X3HHREGN X4HHREGN X5HHREGN X1HHURBN ///
X2HHURBN X3HHURBN X4HHURBN X5HHURBN X3HHLOCL X4HHLOCL X5HHLOCL X3CCLOCL X4CCLOCL ///
X4WCLOCL X5WCLOCL X1HHZIP X2HHZIP X3HHZIP X4HHZIP X5HHZIP X2CCZIP X3CCZIP X4CCZIP ///
X4WCZIP X5WCZIP R1POSAFF R1NEGAFF R1ATNTSK R1SOCIAL P1PRMLNG P1READBO P1WICBFT P1FDSTMP ///
P1MEDICD P1WELFR P1HSTYPE P1HSSIT BCDOBYY BCDOBMM BCMOMMAR R2POSAFF R2NEGAFF R2ATNTSK ///
R2SOCIAL P2PRMLNG P2READBO P2CRMRST P2FDSTMP P2WELFR P2MEDICD P2WICBFT P2HSTYPE P2HSSIT ///
P3READBO P3ADHD P3WICBFT P3HSTYPE P3HSSIT P4READBO P4SCLEXP P4ADHD P4MARSTS P4HSTYPE ///
P4HSSIT P5READBO P5ADHD P5MARSTS W1RRPSU W2RPSU W31RPSU W41RPS W5RPS W1RSTR W2RSTR /// 
X1IFHINC X2IFHINC X3IFHINC X4IFHINC X5IFHINC ///
W31RSTR W41RST W5RST P1MARSTS P2CRMRST P3SMMRST P3MARSTS P4SMMRST P4MARSTS P5SMMRST P5MARSTS ///
X2PRMLNG X3PRMLNG X4PRMLNG X5PRMLNG

/**********************
ACHIEVEMENT TEST SCORES
***********************/
/*MATH*/
replace X3MTHR2 = .i if X3MTHR2 == -9 
replace X3MTHR2 = .p if X3MTHR2 == .
rename X3MTHR2 maththeta05
label variable maththeta05 "Math theta in '05"

/*READING*/
replace X3RTHR2 = .i if X3RTHR2 == -9
replace X3RTHR2 = .p if X3RTHR2 == .
rename X3RTHR2 readtheta05
label variable readtheta05 "Reading theta in '05"

/******************************
INTERMEDIATE COGNITIVE OUTCOMES
*******************************/
/*ABILITY TO PAY ATTENTION*/
replace R1ATNTSK = .i if R1ATNTSK == -9
replace R1ATNTSK = .p if R1ATNTSK == .
replace R2ATNTSK = .i if R2ATNTSK == -9
replace R2ATNTSK = .p if R2ATNTSK == .
replace R2ATNTSK = . if R2ATNTSK == -1

rename R1ATNTSK payattn01
rename R2ATNTSK payattn03 

label variable payattn01 "Child pays attention to tasks - observation '01"
label variable payattn03 "Child pays attention to tasks - observation '03"

label define attnt ///
	0 "Mostly does not pay attention" ///
	1 "Off task half the time" ///
	2 "Typically pays attention" ///
	3 "Constantly pays attention"

recode payattn01 (1=0)(2=0)(3=1)(4=2)(5=3)
label values payattn01 attnt 

recode payattn03 (1=0)(2=0)(3=1)(4=2)(5=3)
label values payattn03 attnt

/*BAYLEY MENTAL SKILL SCALE SCORE*/
replace X1RMTLS = .i if X1RMTLS == -9 
replace X2MTLSCL = .i if X2MTLSCL == -9
replace X1RMTLS = .p if X1RMTLS == . 
replace X2MTLSCL = .p if X2MTLSCL == .
rename X1RMTLS mentalsc01 
rename X2MTLSCL mentalsc03 
label variable mentalsc01 "Cognitive Ability - Bayley Short Form Mental Scale Score '01"
label variable mentalsc03 "Cognitive Ability - Bayley Short Form Mental Scale Score '03"

/*BAYLEY MOTOR SKILL SCALE*/ 
replace X1RMTRS = .i if X1RMTRS == -9 
replace X2MTRSCL = .i if X2MTRSCL == -9
replace X1RMTRS = .p if X1RMTRS == . 
replace X2MTRSCL = .p if X2MTRSCL == .
rename X1RMTRS motorsc01 
rename X2MTRSCL motorsc03 
label variable motorsc01 "Motor Ability - Bayley Short Form Motor Scale Score '01"
label variable motorsc03 "Motor Ability - Bayley Short Form Motor Scale Score '03"

/*********
COVARIATES
**********/
/*CHILD BIRTH WEIGHT*/
replace BCBRTHWT = .p if BCBRTHWT == . 
replace BCBRTHWT = .i if BCBRTHWT == -8

rename BCBRTHWT birthwt

label variable birthwt "Child birth weight in grams, from birth certificate"

/*FAMILY INCOME*/ 
replace X1INCOME = .p if X1INCOME == .
replace X1INCOME = . if X1IFHINC == -1
replace X1INCOME = 2500 if X1INCOME == 1 
replace X1INCOME = 7500 if X1INCOME == 2 
replace X1INCOME = 12500 if X1INCOME == 3
replace X1INCOME = 17500 if X1INCOME == 4
replace X1INCOME = 22500 if X1INCOME == 5
replace X1INCOME = 27500 if X1INCOME == 6
replace X1INCOME = 32500 if X1INCOME == 7
replace X1INCOME = 37500 if X1INCOME == 8
replace X1INCOME = 45000 if X1INCOME == 9
replace X1INCOME = 62500 if X1INCOME == 10
replace X1INCOME = 87500 if X1INCOME == 11
replace X1INCOME = 150000 if X1INCOME == 12
replace X1INCOME = 260000 if X1INCOME == 13

rename X1INCOME faminc01

label variable faminc01 "Household income in 2001 dollars at wave 1"

/*MOTHER EDUCATION*/ 
replace Y1MOMED = .p if Y1MOMED == .
replace Y1MOMED = .i if Y1MOMED == -9
replace Y1MOMED = . if Y1MOMED == -1

rename Y1MOMED momed01

label variable momed01 "Mother's highest education in '01"

label define education ///
	0 "Less than high school diploma" ///
	1 "High school diploma or equivalent" ///
	2 "Vocational/technical degree or Some college" ///
	3 "Bachelor's degree" ///
	4 "Graduate degree"

recode momed01 (1=0)(2=0)(3=1)(4=2)(5=2)(6=3)(7=3)(8=4)(9=4)
label values momed01 education

/*FATHER EDUCATION*/
replace Y1FTHED = .p if Y1FTHED == .
replace Y1FTHED = .i if Y1FTHED == -9
replace Y1FTHED = . if Y1FTHED == -1

rename Y1FTHED daded01

label variable daded01 "Resident father's highest education in '01"

recode daded01 (1=0)(2=0)(3=1)(4=2)(5=2)(6=3)(7=3)(8=4)(9=4)
label values daded01 education

/*HIGHEST PARENT EDUCATION LEVEL*/
gen pared01 = max(momed01, daded01)
label values pared01 education
label variable pared01 "Parent highest education in '01"

/*MOTHER OCCUPATIONAL PRESTIGE*/ 
replace Y1MOMSCR = .p if Y1MOMSCR == .
replace Y1MOMSCR = .i if Y1MOMSCR == -9
replace Y1MOMSCR = . if Y1MOMSCR == -1

rename Y1MOMSCR momocc01

label variable momocc01 "Mother's occupational prestige score in '01"

/*FATHER OCCUPATIONAL PRESTIGE*/ 
replace Y1FTHSCR = .p if Y1FTHSCR == .
replace Y1FTHSCR = .i if Y1FTHSCR == -9
replace Y1FTHSCR = . if Y1FTHSCR == -1

rename Y1FTHSCR dadocc01

label variable dadocc01 "Resident father's occupational prestige score in '01"

/*HIGHEST PARENT OCCUPATIONAL PRESTIGE SCORE*/
gen parocc01 = max(momocc01, dadocc01)
label variable parocc01 "Parent highest occupational prestige in '01"

/*MOTHER EMPLOYMENT STATUS*/
replace Y1HMEMP = .p if Y1HMEMP == .
replace Y1HMEMP = .i if Y1HMEMP == -9
replace Y1HMEMP = . if Y1HMEMP == -1

rename Y1HMEMP momemp01 

label define employment ///
	0 "Not in the labor force" ///
	1 "Less than 35 hours per week" ///
	2 "35 or more per week"

recode momemp01 (4=0)(3=0)(1=2)(2=1)
label values momemp01 employment 
label variable momemp01 "Mother's employment status in '01"

/*FATHER EMPLOYMENT STATUS*/
replace Y1HFEMP = .p if Y1HFEMP == .
replace Y1HFEMP = .i if Y1HFEMP == -9
replace Y1HFEMP = . if Y1HFEMP == -1

rename Y1HFEMP dademp01 

label define employment2 ///
	0 "Not in the labor force" ///
	1 "Currently employed"

recode dademp01 (4=0)(3=0)(1=1)(2=1)
label values dademp01 employment2 
label variable dademp01 "Father's employment status in '01"

/*HOUSEHOLD SIZE*/ 
replace Y1HTOTAL = .p if Y1HTOTAL == .

rename Y1HTOTAL hhtotal01

label variable hhtotal01 "Total number of household members in '01"

/*BIO FATHER IN HOUSEHOLD*/
label define biodad ///
	0 "Biological father not in household" ///
	1 "Biological father is in household"

replace Y1FTHTYP = .p if Y1FTHTYP == . 
recode Y1FTHTYP (2=0)(3=0)(4=0)(5=0)(6=0)(7=0)(9=0)
label values Y1FTHTYP biodad 

rename Y1FTHTYP biodad01

label variable biodad01 "Biological father is in household in '01"

/*MOTHER MARITAL STATUS*/ 
label define marry ///
	0 "Not currently married" ///
	1 "Currently married" 

replace X1MARSTA = .i if X1MARSTA == -9
replace X1MARSTA = . if X1MARSTA == 6
recode X1MARSTA (2=1)(3=0)(4=0)(5=0)
label values X1MARSTA marry 

rename X1MARSTA married01 

label variable married01 "Mother is currently married in '01"

/*MOTHER AGE AT BASELINE*/ 
replace X1HMAG_B = . if X1HMAG_B == -1
rename X1HMAG_B momage01 
label variable momage01 "Mother's age in years in '01"

/*FATHER AGE AT BASELINE */
replace X1HFAG_B = . if X1HFAG_B == -1 
replace X1HFAG_B = .i if X1HFAG_B ==-9
rename X1HFAG_B dadage01 
label variable dadage01 "Father's age in years in '01"

/*HOMEOWNERSHIP*/ 
label define house ///
	0 "Don't own house" ///
	1 "Own house"

replace P1HSSIT = .p if P1HSSIT == . 
replace P1HSSIT = .i if P1HSSIT == -7 
replace P1HSSIT = .i if P1HSSIT == -8 
replace P1HSSIT = .i if P1HSSIT == -9 

recode P1HSSIT (91=0)(7=0)(6=0)(5=0)(4=0)(3=0)(2=0)
label values P1HSSIT house 

rename P1HSSIT house01

label variable house01 "Family owns house in '01"

/*READING WITH CHILD FREQUENCY*/
label define book ///
	0 "Not at all" ///
	1 "Once or twice" ///
	2 "3 to 6 times" ///
	3 "Everyday"

replace P1READBO = .p if P1READBO == . 
replace P1READBO = .i if P1READBO == -9 
replace P1READBO = .i if P1READBO == -8
replace P1READBO = .i if P1READBO == -7

recode P1READBO (1=0)(2=1)(3=2)(4=3)
label values P1READBO book 

rename P1READBO rbooks01 

label variable rbooks01 "In typical week, amount of time mother reads books to child in '01"

/*URRBAN/RURAL RESIDENCE*/
label define urban ///
	1 "Urban, inside urbanized area" ///
	2 "Urban, inside urban cluster" ///
	3 "Rural"

replace X1HHURBN = .p if X1HHURBN == . 
replace X1HHURBN = .i if X1HHURBN ==-9
label values X1HHURBN urban 

rename X1HHURBN urban01 

label variable urban01 "Household urbanicity in '01"

/*CHILD GENDER*/
rename X4CHSEX gender
label variable gender "Child gender"
replace gender=0 if gender==2 

label define gender ///
	0 "Female" ///
	1 "Male"

label values gender gender

/*CHILD RACE/ETHNICITY*/ 
replace Y1CHRACE = .i if Y1CHRACE == -9
rename Y1CHRACE race 
label variable race "Child race/ethnicity"
recode race (1=0)(2=1)(3=2)(4=2)(5=3)(6=4)(7=4)(8=4)

label define race ///
	0 "White, non-Hispanic" ///
	1 "Black, non-Hispanic" ///
	2 "Hispanic" ///
	3 "Asian, non-Hispanic" ///
	4 "Other"

label values race race

/*PRIMARY LANGUAGE SPOKEN AT HOME ENGLISH */
rename X1LANGST prmlang01
label variable prmlang01 "Primary language in household English '01"
recode prmlang01 (1=0)(2=1)

label define english ///
	0 "Not English" ///
	1 "English"

label values prmlang01 english

/*RECEIVED WIC*/ 
label define yes ///
	0 "No" ///
	1 "Yes"

replace P1WICBFT = .i if P1WICBFT == -7 | P1WICBFT == -9
rename P1WICBFT wic01 
label variable wic01 "Mother received WIC in last 12 months"
replace wic01 = 0 if wic01 == 2
label values wic01 yes

/*RECEIVED FOODSTAMPS*/
replace P1FDSTMP = .i if P1FDSTMP == -7 | P1FDSTMP == -8 | P1FDSTMP == -9
rename P1FDSTMP foodst01 
label variable foodst01 "Since child was born, received foodstamps"
replace foodst01 = 0 if foodst01 == 2
label values foodst01 yes

/*RECEIVED MEDICAID*/
replace P1MEDICD = .i if P1MEDICD == -7 | P1MEDICD == -8 | P1MEDICD == -9
rename P1MEDICD medicd01
label variable medicd01 "Since child was born, received medicaid"
replace medicd01 = 0 if medicd01 == 2 
label values medicd01 yes

/*RECEIVED TANF*/
replace P1WELFR = .i if P1WELFR == -7 | P1WELFR == -8 | P1WELFR == -9
rename P1WELFR tanf01 
label variable tanf01 "Since child was born, received TANF"
replace tanf01 = 0 if tanf01 == 2 
label values tanf01 yes

/*CHILD AGE*/ 
gen bday = ym(BCDOBYY, BCDOBMM)
format bday %tm
label variable bday "Child's month and year of birth"

gen date01 = ym(X1ASMTYY, X1ASMTMM)
format date01 %tm
label variable date01 "Date at child assessment in '01"
gen age01 = date01 - bday 
label variable age01 "Child's age at baseline"

gen date05 = ym(X3ASMTYY, X3ASMTMM)
format date05 %tm
label variable date05 "Date at child assessment in '05"
gen age05 = date05 - bday
label variable age05 "Child's age at '05 assessment"

/******
GEODATA
*******/
/*HOUSEHOLD ZCTA*/ 
destring X1HHZIP, generate(zip1)
replace zip1 = .i if zip1 == -9
replace zip1 = .p if zip1 == . 

rename zip1 zip01 

label variable zip01 "Household zip code in '01"

/***************
DESIGN VARIABLES
****************/
/*NORMALIZED SAMPLING WEIGHT*/
sum W1R0
gen sampwt=W1R0/r(mean)
sum sampwt

/*STRATUM AND PSU*/
rename W1RSTR strat
rename W1RRPSU psu

/***********
ID VARIABLES
************/
/*CASE ID*/
gen caseid=I_ID

/*FAMILY AND TWIN ID*/
destring I_ID I_TWINID, replace
gen famid=I_ID if I_TWINID==.
gen temp1=_n
forval i=1/10668 {
	quietly gen temp2=I_TWINID if _n==`i'
	quietly egen temp3=max(temp2)
	quietly replace famid=I_ID if I_ID==temp3 & famid==.
	quietly replace famid=I_TWINID if I_TWINID==temp3 & famid==.
	quietly drop temp2 temp3 
	}
drop temp1
gen twinid=0 if I_TWINID==.
replace twinid=1 if I_TWINID!=.
codebook caseid famid twinid
drop I_ID I_TWINID

/************
TRIM AND SAVE
*************/
order caseid famid twinid strat psu sampwt gender race birthwt age01 age05 maththeta05 readtheta05 payattn01 payattn03 mentalsc01 mentalsc03 motorsc01 motorsc03 faminc01 momed01 daded01 pared01 momocc01 dadocc01 parocc01 momemp01 dademp01 hhtotal01 biodad01 married01 momage01 dadage01 house01 rbooks01 urban01 prmlang01 wic01 foodst01 medicd01 tanf01 zip01

keep caseid famid twinid strat psu sampwt gender race birthwt age01 age05 maththeta05 readtheta05 payattn01 payattn03 mentalsc01 mentalsc03 motorsc01 motorsc03 faminc01 momed01 daded01 pared01 momocc01 dadocc01 parocc01 momemp01 dademp01 hhtotal01 biodad01 married01 momage01 dadage01 house01 rbooks01 urban01 prmlang01 wic01 foodst01 medicd01 tanf01 zip01

save "${data_directory}eclsb\v02_eclsb_nvars.dta", replace 

/*********************
DESCRIPTIVE STATISTICS   
**********************/
codebook*

log close 
