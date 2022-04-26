capture log close
capture clear all
set more off

global projdir "C:\Users\wodtke\Desktop\projects\nhood_mediation_toxins\" 

log using "${projdir}programs\_LOGS\08_create_figure_2.log", replace

/*****************************************
DO-FILE NAME: 08_create_figure_2.do
AUTHOR: GW
PURPOSE: create chloropleths of nh poverty
and selected toxins for Cook County, IL
******************************************/

/*****************************
CREATE DTA FILE FROM SHP FILES
******************************/
capture erase "${projdir}data\_TEMP\ILcoord.dta"
capture erase "${projdir}data\_TEMP\ILdb.dta"

shp2dta using "${projdir}data\shpfiles\tr17_d00.shp", ///
	database("${projdir}\data\_TEMP\ILdb.dta") ///
	coordinates("${projdir}\data\_TEMP\ILcoord.dta") genid(id)

/***************************************
CREATE NCDB FILE IN 2000 FIPS BOUNDARIES
****************************************/
import delimited "${projdir}data\ncdb\csvFiles\ncdb_raw.csv", encoding(ISO-8859-1) 

rename trctpop0 trctpop_2000 
rename trctpop1 trctpop_2010 
rename povrat0d povratd_2000 
rename povrat0n povratn_2000 
rename povrat1ad povratd_2010 
rename povrat1an povratn_2010 

keep geo2010 state region division county trctpop* povratd* povratn* 

save "${projdir}data\_TEMP\ncdb_temp.dta", replace 
clear 

insheet using "${projdir}data\crosswalks\us2010trf.txt", clear delim(",")

keep v4 v13 v29 v30

rename v4 geo2000
rename v13 geo2010
rename v29 hupct00
rename v30 hupct10

label var geo2000 "Concatenation of 2000 state, county and tract for Census 2000 Tract"
label var geo2010 "Concatenation of state, county and tract for 2010 Census Tracts"
label var hupct00 "Percent of the HU00 this record contains"
label var hupct10 "Percent of the HU10 this record contains"

format geo2000 %15.0g
format geo2010 %15.0g

bysort geo2010: gen freq=_n
bysort geo2010: egen maxfreq=max(freq)

gen xwgt=hupct10/100 if maxfreq>1
replace xwgt=hupct00/100 if maxfreq==1

egen sumwt = sum(xwgt), by(geo2010)

keep geo* xwgt

merge m:1 geo2010 using "${projdir}data\_TEMP\ncdb_temp.dta"
drop if _merge == 1
drop _merge

sort geo2010
foreach var of varlist trctpop_2000-povratd_2010 {
	replace `var'= `var'*xwgt
	} 

collapse (sum) trctpop_2000-povratd_2010 (firstnm) state division county, by(geo2000)
sort geo2000

capture macro drop vars 
global vars povratd povratn trctpop
			  
foreach v of global vars { 
	forval t=2001/2009 { 
		gen `v'_`t'=`v'_2000+(`v'_2010-`v'_2000)*((`t'-2000)/(2010-2000)) 
		} 
	} 

reshape long povratd_ povratn_ trctpop_ , i(geo2000) j(year) 

foreach v in povratd povratn trctpop { 
			 rename `v'_ `v' 
			 } 

gen nhpovrt=povratn/povratd 
sum nhpovrt 

keep geo2000 state division county year nhpovrt 
keep if inrange(year,2001,2007)

reshape wide nhpovrt, i(geo2000) j(year) 

format %13.0g geo2000
tostring geo2000, replace format(%13.0g)
gen tract = ""
replace tract = substr((geo2000),6,6)
replace tract = substr((geo2000),6,4) if real(substr((geo2000),10,2))==0
label var tract "2000 FIPS TRACT ID"
destring tract, replace

keep geo2000 state division county tract *2001 

save "${projdir}data\_TEMP\ncdb_temp.dta", replace
clear 

/****************************************
CREATE CACES FILE IN 2000 FIPS BOUNDARIES
*****************************************/
import delimited "${projdir}data\caces\caces_raw.csv", encoding(ISO-8859-1)

format fips %15.0g
sort fips year pollutant

rename pred_wght poll_ 
rename fips geo2010

keep geo2010 pollutant year poll_ 
sort geo2010 year pollutant
order geo2010 year pollutant poll_  

reshape wide poll_, i(geo2010 year) j(pollutant) string

rename poll_co co_
rename poll_no2 no2_
rename poll_o3 o3_
rename poll_pm10 pm10_
rename poll_pm25 pm25_
rename poll_so2 so2_

reshape wide co_ no2_ o3_ pm10_ pm25_ so2_, i(geo2010) j(year)

save "${projdir}data\_TEMP\caces_temp.dta", replace 
clear 

insheet using "${projdir}data\crosswalks\us2010trf.txt", clear delim(",")

keep v4 v13 v22 v24

rename v4 geo2000
rename v13 geo2010
rename v22 landpct00
rename v24 landpct10

format geo2000 %15.0g
format geo2010 %15.0g

bysort geo2010: gen freq=_n
bysort geo2010: egen maxfreq=max(freq)

gen xwgt=landpct10/100 if maxfreq>1
replace xwgt=landpct00/100 if maxfreq==1

sort geo2000
egen sumwt = sum(xwgt), by(geo2000)
replace xwgt = xwgt/sumwt
drop sumwt

keep geo* xwgt

merge m:1 geo2010 using "${projdir}data\_TEMP\caces_temp.dta"
drop if _merge == 1 
drop _merge

sort geo2010
foreach var of varlist co_1980-so2_2015 {
	replace `var'= `var'*xwgt
	} 

collapse (sum) co_1980-so2_2015, by(geo2000)

keep geo2000 *_2001

save "${projdir}data\_TEMP\caces_temp.dta", replace 

/***************************************
CREATE RSEI FILE IN 2000 FIPS BOUNDARIES
****************************************/
use "${projdir}data\rsei\raw\2001neuro_2000tracts.dta", clear
gen year=.
replace year = 2001

rename *, lower

rename fips geo2000
replace geo2000 = " " if regexm(geo2000, "nodata")
gen g = substr(geo2000,1,2)
order g
gen h = substr(geo2000,4,3)
order g h 
gen i = substr(geo2000,8,6)
order g h i
gen geo = g + h + i
order geo
drop geo2000 g h i
rename geo geo2000

destring geo2000, replace
format geo2000 %13.0f

drop *tc
drop *score

rename (*c)(*)
rename (*)(*_)
rename geo2000_ geo2000
rename year_ year

unab vlist: chem*
di "`vlist'"
local j "year"
di "`j'"

levelsof `j', local(J) //
foreach var of varlist `vlist' {
	foreach j of local J {
		local newlist `newlist' `var'`j'
		local lablist "`lablist' `"`:variable label `var'' (`j')"'"
	}
}

reshape wide chem*, i(geo2000) j(year)

foreach new of local newlist {
	gettoken lab lablist : lablist
	quietly label var `new' "`lab'"
}

keep geo2000 *_2001

foreach i of var chem* {
	quietly replace `i'=`i'*1000
	}

save "${projdir}data\_TEMP\rsei_temp.dta", replace 

/*********************************************
MERGE NCDB, CACES, RSEI, AND CARTOGRAPHIC DATA
**********************************************/
use "${projdir}data\_TEMP\ncdb_temp.dta", clear
destring geo2000, replace

merge 1:1 geo2000 using "${projdir}data\_TEMP\caces_temp.dta" 
keep if _merge==1 | _merge==3
drop _merge 
merge 1:1 geo2000 using "${projdir}data\_TEMP\rsei_temp.dta"
keep if _merge==1 | _merge==3
drop _merge

keep if state==17 & county==31 //Cook county
drop state division county

save "${projdir}data\_TEMP\ncdb_caces_rsei_temp.dta", replace

use "${projdir}\data\_TEMP\ILdb.dta", clear

foreach x in STATE COUNTY TRACT NAME LSAD_TRANS {
	destring `x', replace 
	}

keep if CO == 31 //Cook county 
sort TRACT

rename TRACT tract

merge m:1 tract using "${projdir}\data\_TEMP\ncdb_caces_rsei_temp.dta"
keep if _merge == 3
drop _merge

sort id
foreach x in nhpovrt2001 pm10_2001 no2_2001 chem355_2001 {
	replace `x'=(`x'[_n-1]+`x'[_n+1])/2 if `x'==.
	replace `x'=`x'[_n-1] if `x'==.
	replace `x'=`x'[_n+1] if `x'==.
	}

/**************************************
CREATE CHLOROPLETHS FOR COOK COUNTY, IL
***************************************/
//Panel A: Poverty
replace nhpovrt2001=round(nhpovrt2001,0.01)
grmap nhpovrt2001 using "${projdir}\data\_TEMP\ILcoord.dta", ///
	id(id) fcolor(Greys2) clm(q) cln(10) ndf(bluishgray) mosize(vvthin) ///
	title("A. Poverty Rate", size(medsmall)) ///
	legtitle("Deciles")

quietly graph save "${projdir}\figures\_TEMP\fig2a.gph", replace

//Panel B: PM10
replace pm10_2001=round(pm10_2001,0.1)
grmap pm10_2001 using "${projdir}\data\_TEMP\ILcoord.dta", ///
	id(id) fcolor(Greys2) clm(q) cln(10) ndf(bluishgray) mosize(vvthin)  ///
	title("B. PM10 (ug/m^3)", size(medsmall)) ///
	legtitle("Deciles")

quietly graph save "${projdir}\figures\_TEMP\fig2b.gph", replace

//Panel C: Nitrogren Dioxide
replace no2_2001=round(no2_2001,0.1)
grmap no2_2001 using "${projdir}\data\_TEMP\ILcoord.dta", ///
	id(id) fcolor(Greys2) clm(q) cln(10) ndf(bluishgray) mosize(vvthin)  ///
	title("C. NO2 (ppb)", size(medsmall)) ///
	legtitle("Deciles")

quietly graph save "${projdir}\figures\_TEMP\fig2c.gph", replace

//Panel D: Manganese
replace chem355_2001=round(chem355_2001,0.1)
grmap chem355_2001 using "${projdir}\data\_TEMP\ILcoord.dta", ///
	id(id) fcolor(Greys2) clm(q) cln(10) ndf(bluishgray) mosize(vvthin)  ///
	title("D. Mn (ng/m^3)", size(medsmall)) ///
	legtitle("Deciles")

quietly graph save "${projdir}\figures\_TEMP\fig2d.gph", replace

//Combine Panels
graph combine ///
	"${projdir}\figures\_TEMP\fig2a.gph" ///
	"${projdir}\figures\_TEMP\fig2b.gph" ///
	"${projdir}\figures\_TEMP\fig2c.gph" ///
	"${projdir}\figures\_TEMP\fig2d.gph", ///
	col(2) row(2) ysize(6.5) xsize(9) scheme(s2mono) imargin(medsmall)

graph export "${projdir}\figures\figure_2.emf", as(emf) replace

erase "${projdir}\data\_TEMP\ILcoord.dta"
erase "${projdir}\data\_TEMP\ILdb.dta"
erase "${projdir}data\_TEMP\ncdb_temp.dta"
erase "${projdir}data\_TEMP\caces_temp.dta"
erase "${projdir}data\_TEMP\rsei_temp.dta"
erase "${projdir}data\_TEMP\ncdb_caces_rsei_temp.dta"
erase "${projdir}\figures\_TEMP\fig2a.gph"
erase "${projdir}\figures\_TEMP\fig2b.gph"
erase "${projdir}\figures\_TEMP\fig2c.gph"
erase "${projdir}\figures\_TEMP\fig2d.gph"

log close
