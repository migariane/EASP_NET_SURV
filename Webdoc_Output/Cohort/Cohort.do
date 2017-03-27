webdoc init Cohort, replace logall header(width(800px) stscheme(classic))
/***
<html>
<head>
<title>AGE-STANDARDISED NET SURVIVAL BY DEPRIVATION COHORT APPROACH</title>
</head>
***/

/***
<body>
<h2>TUTORIAL: AGE-STANDARDISED NET SURVIVAL BY DEPRIVATION</br>
<p></p>
COHORT APPROACH</br>
<p></p>
EASP-COURSE, GRANADA, 28-29 MARCH 2017</br>
<p></p>
Miguel Angel Luque Fernandez, BSc, MA, MPH, MSc, GStat(UK), PhD</br>
Assistant Professor of Epidemiology</br>
Cancer Survival Group, LSHTM, London, UK</br>
<p></p>
Rhea Harewood, BSc, MSc</br>
Research Fellow in Epidemiology</br>
Cancer Survival Group, LSHTM, London, UK</br>
<p></p>
</h2>
***/

/***
<h2>Content</h4>
<h3>I) Data consistency and setting time</br>
<p></p>
II) Net Survival Estimation</br>
<p></p>
III) Dealing with ICSS WEIGHTS</br>
<p></p>
IV) Age-Standardised Five-year Net Survival Estimation, chohort 1971</br>
	 <p style="text-indent: 3em">-Life Table</p>
	 <p style="text-indent: 3em">-STNS Stata command</p>
V) Age-Standardised Five-year Net Survival Inference, chohort 1971</br>
<p></p>
VI) Age-Standardised Five-year Net Survival by Deprivation, cohort 1971</br></h3>
<p></p>
<p>Note: the data used in the tutorial has been modified from the original source for a teaching purpose and represent breast cancer incident cases between 1971 and 2001 in England</p>
<p>Note: the interpretation of the results is not applicable to the real-world setting</p>

***/

/***
<h3>I) Data consistency and setting time</h3>
***/

/***
<p>Setting your path and working directory</p>
***/
clear
set more off
cd "/Users/MALF/Desktop"

/***
<p>Loading the data</p>
***/ 
use breast_stns.dta, clear
set more off
//browse

/***
<p>Describing the data</p>
***/
describe
summarize

/***
<p>Calendar year at diagnosis</p>
***/
gen year = year(diagmdy)
labe var year "calendar year at diagnosis"
tab year

/***
<p>Canlendar year at exit (last known vital status)</p>
<p>Note that in Stata the 1st of January of 1960 = 0</p>
<p>How you will check the consistency of the data?</p>
***/
gen eyear = year(finmdy)
labe var eyear "calendar year last follow-up"
tab eyear dead
sum dead finmdy if dead==1 & (finmdy>15705 & finmdy<=16070)

/***
<p>Five age groups are needed for standardisation </p>
***/
sum agediag, det
egen agegr =cut(agediag), at(0 45(10)75 100) icodes
recode agegr 0=1 1=2 2=3 3=4 4=5
tabstat agediag, statistics( min max ) by(agegr)
label var agegr "Five age groups for standardisation"

/***
<p>Setting time (note origin and entry)</p>
***/
stset finmdy, fail(dead==1) origin(time diagmdy) enter(time diagmdy) 

/***
<p>Checking the assumptions of time (note _t0)</p>
***/
list diagmdy finmdy _t0 _t _d _st in 1/10
scalar  _t = (finmdy - diagmdy) in 1
display _t

/***
<h3>II) Net Survival Estimation: LIFE TABLE and STNS</h3>

<p>The STRUCTURE OF THE LIFE TABLE (understanding it is really IMPORTANT)</p>
<p>Note: the unit of _age. stns needs _age in days)</p>
***/
preserve
clear
use Lifetable_stns.dta
list in 1/10
display 11*365.25
restore 

/***
<p>In the breast cancer incident dataset we have to generate age at diagnosis in days to merge it with the life table data (rate in days)</p>
***/
list diagmdy finmdy _t0 _t _d _st agediag in 1/10
gen agediagindays = agediag*365.25
label var agediagindays  "Age at diagnosis in days for Net Survival estimation"
list diagmdy finmdy _t0 _t _d _st agediag agediagindays in 1/10

/***
<p>Net Survival Estimates at 1, 2, 3, 4, 5, and 10 years after diagnosis</p>
***/
stns list using lifetable_stns.dta if year(diagmdy)==1971, ///
	age(agediagindays=_age) period(diagmdy=yearindays) ///
	strata(dep) rate(rate_day) ///
	at(1/4 5(5)10, scalefactor(365.25) unit(year)) end_followup(3625.5)	///	
	saving(period, replace)
preserve
clear
use period
gen year=time/365.25
twoway (connected survival year, sort msymbol(none)), ///
	yscale(range(0 1)) ylabel(0(.2)1, labels angle(horizontal) format(%9.1g)) ///
	ytitle("Net survival") xscale(range(0 10)) xlabel(, val angle(horizontal)) ///
	xtitle("Years since diagnosis") 
rm period.dta
restore
webdoc graph, caption(Figure 1.Five Years Net Survival for the Cohort of Breast Cancer 1971) cabove ///
width(1000)	

/***
<p>Net Survival for the cohort of cancer incident cases diagnosed in 1971 and followed-up for 5 years</p>
<p>All what you have to know about STNS:</br>
1. if year(diagmdy) == 1971</br>
2. age and period (specification of the linking variables to merge 1:1 the life table and the incident cases)</br>
3. strata and rate (stratified Net Survial by deprivation)</br> 
4. at (years of follow up to compute the Net Survival)</br>
5. scale factor and unit</br>
6. For display end_follow up (in days) and by (agegr dep)</br>
7. saving option (really important!)</br>
8. Please: read carefully the Stata stns help file and the Stata <a href="http://www.stata-journal.com/article.html?article=st0326">Journal article</a></p>
***/
stns list using lifetable_stns.dta if year(diagmdy)==1971, ///
	age(agediagindays=_age) period(diagmdy=yearindays) ///
	strata(dep) rate(rate_day) ///
	at(1 2 3 4 5, scalefactor(365.25) unit(year)) end_followup(5480) by(agegr dep) ///
	saving(ASNetcohort_1971, replace)

/***
<p>Loading Net Survival estimates and checking consistency</p>
***/	
clear 
use ASNetcohort_1971.dta
list time if time==5480 & dep==1 & agegr==1 //checking consistency

/***
<p>Undesrtanding the results derived from STNS Stata command</p>
***/	
describe
list dep agegr time survival std_err lower_bound upper_bound in 1/10

/***
<p>Keep results for just five-year Net Survival estimates</p>
***/
/***
<p>Q: What do you have to make in order to get 10 years survival?</p>
<p>A: drop if time > ???. Remember: 365.25 x Years (display 365.25*10)
***/
display 365.25*10
display 365.25*5
drop if time > 1826.25
bysort dep agegr (time): keep if _n == _N

/***
<h3>III) Dealing with ICSS WEIGHTS</h3>
<p>International Cancer Survival Standard (ICSS) weights </br> 
(Corazziari I, Quinn M, Capocaccia R. Eur J Cancer. 2004; 40: 2307-16. Standard cancer patient population for age standardising survival ratios.)</p>
***/

/***
<p>Q: Where do you can get the information of the weights for other cancer sites?</p>
<p>A: Check out the Stata do file provide for the exercise and Corazziari et al. European Journal of Cancer. 2004</p>
***/
gen weight=. 
*Standard cancer population one (stomach, colon, rectum, liver, lung, breast, ovary, leukaemia)
replace weight=.07 if agegr==1
replace weight=.12 if agegr==2
replace weight=.23 if agegr==3
replace weight=.29 if agegr==4
replace weight=.29 if agegr==5

/***
<h3>IV) Age-Standardised Net Survival estimation</h3>
<p>Weighted NET SURVIVAL estimate: from Corazziari et al. European Journal of Cancer. 2004</p>
***/
bysort /*insert relevant variables*/ dep (agegr): gen s1ASN = (surv*weight)
list agegr dep surv weight s1ASN
display .70528846*.07
bysort /*insert relevant variables*/ dep (agegr): egen  ASNS = sum(s1ASN)
list agegr dep surv weight s1ASN ASNS 

/***
<h3>Remember: <font color="blue"><mark>H(t) = -ln(S(t)) ==> exp(-H(t)) = S(t)</mark></font></h3>
***/ 
/***
<h3>V) Age-Standardised Net Survival Statistical Inference</h3>
<p>Weighted STANDARD ERROR of the net survival estimate. The formula for the the standard error for net survival (se_ns) is derived from the DELTA METHOD
based on Clayton and Hills. Statistical Models in Epidemiology, 1993</p>
***/ 
gen ns=exp(-cum_hazard) //using the Delta method we ned the cummulative hazard H(t).
corr ns surv //checking consistency
gen se_ns=ns*std_err //where std_err is the standard error of the cumulative hazard and ns is the survival estimate from stns
bysort /*insert relevant variables*/ dep (agegr): gen seASN = sqrt(sum((se_ns*weight)^2))
//Keep age-standardise estimate by deprivation
bysort /*insert relevant variables*/ dep (agegr): keep if _n == _N

/***
<p>95%CIs from Corazziari et al. European Journal of Cancer. 2004 and Clayton and Hills. Statistical Models in Epidemiology, 1993</p>
***/
gen L95CI=(ASNS/exp(1.96*seASN/ASNS))
gen U95CI=(ASNS*exp(1.96*seASN/ASNS))

/***
<h3>VI) Age-standardised Five-year Net Survival by Deprivation, for the Cohort 1971</h3>
***/
list dep survival ASNS L95CI U95CI 
eclplot ASNS L95CI U95CI dep, hori estopts(msize(vlarge)) ciopts(msize(vlarge)) yscale(range(1 6)) xline(0,lpattern(dot)) xtitle("Age-Standardised Net Survival")
webdoc graph, caption(Figure 2. Age-Standardised Five-year Net Survival for the Breast Cancer Cohort of 1971) cabove ///
width(1000)
/***

<h3><font color="blue">THANK YOU FOR YOUR ATTENTION</font></h3>

</body>
</html>
***/
webdoc close
rm ASNetcohort_1971.dta
//webdoc do Webdoc_Cohort_Age_Sta_Net_Surv.do
