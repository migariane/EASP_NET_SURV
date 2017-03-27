*********************************************************
*AGE-STANDARDISED NET SURVIVAL BY DEPRIVATION 
*PERIOD APPROACH
*EAPS-SEMINAR
*CANCER SURVIVAL ANALYSIS USING POPULATION-BASED DATA
*GRANADA, 28-29 MARCH 2017
*Miguel Angel Luque-Fernandez & Rhea Harewood
********************************************************

ssc install stns
which stns
//I Data consistency and setting time
//Setting your path and working directory
clear
set more off
cd "C:\Users\invitado-easp\Desktop\ASNS"
capture log using "Period.log"

//Loading the data

use breast_stns.dta, clear
set more off
//browse

//Describing the data
describe
summarize

//Calendar year at diagnosis
gen year = year(diagmdy)
labe var year "calendar year at diagnosis"
tab year

/***
Canlendar year at exit (last follow-up)
Note that in Stata the 1st of January of 1960 = 0
How you will check the consistency of the data?
***/
gen eyear = year(finmdy)
labe var eyear "calendar year last follow-up"
tab eyear dead
sum dead finmdy if dead==1 & (finmdy>15705 & finmdy<=16070)

/***
Five-band age groups needed for standardisation
***/
sum agediag, det
egen agegr =cut(agediag), at(0 45(10)75 100) icodes
recode agegr 0=1 1=2 2=3 3=4 4=5
tabstat agediag, statistics(min max) by(agegr)
label var agegr "5-band age groups for standardisation"

/***
SETTING time for the five years calendar PERIOD time 1981-1985 (note the ORIGIN, ENTRY and EXIT options)
Understanding your approach: note _t0 is not 0 as before in the cohort analysis
ORIGIN(time diagmdy): you are setting the date of cancer diagnosis as the ANALYSIS TIME
ENTER(time mdy(1,1,1980): the ONSET RISK start a the date specified by the user
EXIT(time mdy(12,31,1984): the exact date when subjects exit from the analysis STOP COUNTING person-time at risk
***/
stset finmdy, fail(dead==1) origin(time diagmdy) enter(time mdy(1, 1, 1981)) exit(time mdy(12, 31, 1985))
sort year
list diagmdy finmdy _t0 _t _d _st in 1/10 
list diagmdy finmdy _t0 _t _d _st in 130100/130110 
list diagmdy finmdy _t0 _t _d _st in 180100/180110 

/***
II) Net Survival Estimation: LIFE TABLE and STNS
The STRUCTURE OF THE LIFE TABLE (understanding it is really IMPORTANT)
Note: the unit of _age, stns needs _age in days
***/
preserve
clear
use Lifetable_stns.dta
list in 1/10
display 11*365.25
restore 

/***
In the breast cancer incident dataset we have to generate age at diagnosis in days to merge it with the life table data (rate in days)
***/
list diagmdy finmdy _t0 _t _d _st agediag in 1/10
gen agediagindays = agediag*365.25
label var agediagindays  "Age at diagnosis in days for Net Survival estimation"
list diagmdy finmdy _t0 _t _d _st agediag agediagindays in 1/10

/***
Understanding your approach:
Note _t0 is not 0 as before in the cohort analysis
Our focus in the analysis is just five-years calendar period
***/
preserve
keep if _st==1
list diagmdy finmdy _t0 _t _d _st year in 1/10 if year==1971
display 5046-3220
display %2.0f 365.24*5
restore

/***
Net Survival Estimates at 1, 2, 3, 4, 5, and 10 years after diagnosis
Q: Can you remark and enumerate the difference between the PERIOD and COHORT approaches?
A: Yes </br>
1. The setting of time was different for the PERIOD approach as we specified the study ENTRY and EXIT dates
2. The opotion IF after calling the life table data is not specified in the PERIOD approach
3. The option <font color="blue">begintime</font> was not specified before as in the cohort approach the ONSET RISK was the DATE of CANCER DIAGNOSIS
4. The begintime option is needed to compute the inverse-probability censoring weights (IPCW)
5. With the IPCW our Net Survival estimate is independent from other causes of death accounting for CENSORING
***/
stns list using lifetable_stns.dta, ///
	age(agediagindays=_age) period(diagmdy=yearindays) ///
	begintime(origin) strata(dep) rate(rate_day) ///
	at(.5 1 5(5)10, scalefactor(365.25) unit(year)) end_followup(3625.5)	///
	saving(period, replace)

/***
Comparing COHORT and PERIOD appraches
Q: Can you explain the differences?
***/
quietly stset finmdy, fail(dead==1) origin(time diagmdy) enter(time diagmdy)
stns list using lifetable_stns.dta if year(diagmdy)==1971, ///
	age(agediagindays=_age) period(diagmdy=yearindays) ///
	begintime(origin) strata(sex dep) rate(rate_day) ///
	at(.5 1 5(5)10, scalefactor(365.25) unit(year)) end_followup(3625.5) ///
 	saving(cohort, replace)
use cohort, clear
gen str11 approach="Cohort 1971"
append using period
replace approach="Period 1981-85" if approach==""
gen year=time/365.25
twoway (connected survival year if approach=="Cohort 1971", sort msymbol(none)) ///
	(connected survival year if approach=="Period 1981-85", sort msymbol(none)), ///
	yscale(range(0 1)) ylabel(0(.2)1, labels angle(horizontal) format(%9.1g)) ///
	ytitle("Net survival") xscale(range(0 10)) xlabel(, val angle(horizontal)) ///
	xtitle("Years since diagnosis") ///
	legend(order(1 "Cohort 1971" 2 "Period 1981-85") row(1))
rm cohort.dta
rm period.dta

/***
Q: Can you explain the differences?
A: Advances in Medical Treatment!
***/

/***
Five-years Net Survival Estimates at 1, 2, 3, 4, 5, and 10 years after diagnosis for the Period 1.1.1981 31.12.1985
Q: Which other options I have to know to understand STNS?
1. age and period (specification of the linking variables to merge 1:1 the life table and the cancer incident cases)
2. strata and rate (stratified Net Survial by deprivation)
3. at (years of follo-up to compute the Net Survival)
4. scale factor and unit
5. For display end_followup (in days) and by(agegr dep)
6. saving option (really important!). Please use a meaningful name
7. Please: read carefully the Stata stns help file and the Stata "http://www.stata-journal.com/article.html?article=st0326"
***/

use breast_stns, clear
egen agegr =cut(agediag), at(0 45(10)75 100) icodes
recode agegr 0=1 1=2 2=3 3=4 4=5
gen agediagindays = agediag*365.25
quietly stset finmdy, fail(dead==1) origin(time diagmdy) enter(time mdy(1, 1, 1981)) exit(time mdy(12, 31, 1985))
stns list using lifetable_stns.dta, ///
	begingtime(origin) age(agediagindays=_age) period(diagmdy=yearindays) ///
	strata(sex dep) rate(rate_day) ///
	at(1 2 3 4 5, scalefactor(365.25) unit(year)) end_followup(5480) by(agegr dep) ///
	saving(ASNetPeriod_1981-1985, replace)
clear 
use ASNetPeriod_1981-1985.dta	
drop if time > 1826.25
bysort dep agegr (time): keep if _n == _N

/***
III) Dealing with ICSS WEIGHTS
Weights from International Cancer Survival Standards (ICSS). Corazziari et al. European Journal of Cancer. 2004
Q: Where do you can get the information of the weights for other cancer sites?
A: Check out the Stata do file provide for the exercise and the article referred above
***/
gen weight=. 
*Standard cancer population one (stomach, colon, rectum, liver, lung, breast, ovary, leukaemia)
replace weight=.07 if agegr==1
replace weight=.12 if agegr==2
replace weight=.23 if agegr==3
replace weight=.29 if agegr==4
replace weight=.29 if agegr==5

/*Standard cancer population for Prostate
replace weight=.19 if age5STD==1 
replace weight=.23 if age5STD==2 
replace weight=.29 if age5STD==3 
replace weight=.23 if age5STD==4 
replace weight=.06 if age5STD==5 

*Standard cancer population for Cervix
replace weight=.28 if age5STD==1 
replace weight=.17 if age5STD==2 
replace weight=.21 if age5STD==3 
replace weight=.20 if age5STD==4 
replace weight=.14 if age5STD==5 

*For childhood leukaemia will use equal weights for each age group
replace weight=1/3 if age5STD==1 
replace weight=1/3 if age5STD==2 
replace weight=1/3 if age5STD==3 */

/***
IV) Age-Standardised Net Survival estimation
Weighted NET SURVIVAL estimate: from Corazziari et al. European Journal of Cancer. 2004
***/
bysort /*insert relevant variables*/ dep (agegr): gen s1ASN = (surv*weight)
list agegr dep surv weight s1ASN
display .70528846*.07
bysort /*insert relevant variables*/ dep (agegr): egen  ASNS = sum(s1ASN)
list agegr dep surv weight s1ASN ASNS 

/***
Remember: H(t) = -ln(S(t)) ==> exp(-H(t)) = S(t)
V) Age-Standardised Net Survival Statistical Inference
Weighted STANDARD ERROR of the net survival estimate. The formula for the the standard error for net survival (se_ns) is derived from the DELTA METHOD
based on Clayton and Hills. Statistical Models in Epidemiology, 1993
***/ 
gen ns=exp(-cum_hazard) //using the Delta method we ned the cummulative hazard H(t).
corr ns surv //checking consistency
gen se_ns=ns*std_err //where std_err is the standard error of the cumulative hazard and ns is the survival estimate from stns
bysort /*insert relevant variables*/ dep (agegr): gen seASN = sqrt(sum((se_ns*weight)^2))
//Keep age-standardise estimate by deprivation
bysort /*insert relevant variables*/ dep (agegr): keep if _n == _N

/***
95%CIs from Corazziari et al. European Journal of Cancer. 2004 and Clayton and Hills. Statistical Models in Epidemiology, 1993
***/
gen L95CI=(ASNS/exp(1.96*seASN/ASNS))
gen U95CI=(ASNS*exp(1.96*seASN/ASNS))

/***
VI) Age-standardised Five years Net Survival by Deprivation, for the Cohort 1971
***/
list dep survival ASNS L95CI U95CI 
eclplot ASNS L95CI U95CI dep, hori estopts(msize(vlarge)) ciopts(msize(vlarge)) yscale(range(1 6)) xline(0,lpattern(dot)) xtitle("Age-Standardised Net Survival")

log close
rm ASNetPeriod_1981-1985.dta

//Thank YOU! So, Beautiful is Granada, isn't it?
