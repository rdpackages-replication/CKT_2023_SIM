/*###################--------------------------##############
# Replication files for "A Guide to Regression Discontinuity Designs in Medical Applications"
# ART application
# by Matias D. Cattaneo, Luke Keele, and Rocio Titiunik
# last modified: May 2022
###################--------------------------##############*/

/* To install necessary packages for RD analysis, execute the following commands 

net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace
net install rdlocrand, from(https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata) replace
net install rddensity, from(https://raw.githubusercontent.com/rdpackages/rddensity/master/stata) replace
net install lpdensity, from(https://raw.githubusercontent.com/nppackages/lpdensity/master/stata) replace

For more details, see https://rdpackages.github.io
*/
global covariates age1 age2  age3 age4 age5 age6 age7 age8 qtr1 qtr2 qtr3 qtr4 qtr5 qtr6 female clinic_a clinic_b clinic_c
clear all
use "CKT_2023_SIM--ART.dta", clear

/* --------------------------------
 Figure 1(a): Score Histogram
 ----------------------------------*/
histogram cd4, width(10) xline(350) frequency

 /* --------------------------------
 Figure 1(b): Outcome Plot
 ----------------------------------*/
rdplot visit_test_6_18 cd4, p(3) c(350) graph_options(yscale(r(0 1)) xscale(r(0 1000)))

/* --------------------------------
 Figure 3(a): Treatment Assignment
----------------------------------*/
rdplot below cd4, c(350) nbins(10 100) binselect(es) graph_options(yscale(r(0 1)) xscale(r(0 1000))) p=0 
 
/* --------------------------------
  Figure 3(b): Treatment Take-up
----------------------------------*/
rdplot art_6m cd4, c(350) p(3) graph_options(yscale(r(0 1)) xscale(r(0 1000)))

/* --------------------------------
  Table 1: Continuity-based RD Estimates with Robust Bias Corrected Inference
----------------------------------*/
rdrobust visit_test_6_18 cd4, c(350) all fuzzy(art_6m)

* save bandwidths for falsification test below
global h_original = e(h_l)
global b_original = e(b_l)

local h =  e(h_l)
local b =  e(b_l)

* Intention-to-treat effect in same bandwdiths
rdrobust visit_test_6_18 cd4, c(350) all h(`h') b(`b')

/* --------------------------------
  Table 2: Local Randomization RD Estimates and Inference
----------------------------------*/

** Window Selection  (not shown in paper)
rdwinselect cd4 age1 age2 age3 age4 age5 age6 age7 age8 qtr1 qtr2 ///
qtr3 qtr4 qtr5 qtr6 clinic_a clinic_b clinic_c, c(350) seed(987) reps(1000) wstep(1)

* In order to always use the same number of observations, keep only observations non-missing in covariates, Y, and D
keep if cd4!=. & art_6m !=. & visit_test_6_18 !=.
foreach var of global covariates {
	keep if `var'!=.
}
rdrandinf art_6m cd4, c(350) seed(5023) wl(346) wr(354)  ci(.05 -.5(.01).5) wmass
rdrandinf visit_test_6_18 cd4, c(350) seed(5023) wl(346) wr(354)  ci(.05 -.5(.01).5) wmass 
rdrandinf visit_test_6_18 cd4, c(350) seed(5023) wl(346) wr(354)  ci(.05 -2(.1)2) wmass fuzzy(art_6m tsls)
* Note: small discrepancies between R and Stata in the TSLS large-sample confidence intervals are due to differences in the implementation of dependencies

/* --------------------------------
   Figure 4(a): Histogram near cutoff: score in [345, 355]
----------------------------------*/
use "CKT_2023_SIM--ART.dta", clear
histogram cd4 if cd<=355 & cd4>=345, xline(350) frequency bin(11)

/* --------------------------------
    Figure 4(b): Histogram and local polynomial density estimate for score 
----------------------------------*/
rddensity cd4, c(350) vce(jackknife) plot

/* --------------------------------
   Table 3: Continuity-Based ITT RD Estimates for Predetermined Covariates
----------------------------------*/

foreach var of global covariates {
	rdrobust `var' cd4, c(350) bwselect(mserd)
}

/* --------------------------------
  Table 4: Continuity-Based Bandwidth Sensitivity Diagnostic
----------------------------------*/
local h =  $h_original + 10
local b =  $b_original + 10

rdrobust visit_test_6_18 cd4, c(350) all h(`h') b(`b') masspoints(adjust)
rdrobust visit_test_6_18 cd4, c(350) all h(`h') b(`b') fuzzy(art_6m) masspoints(adjust)

local h =  $h_original + 35
local b =  $b_original + 35

rdrobust visit_test_6_18 cd4, c(350) all h(`h') b(`b') masspoints(adjust)
rdrobust visit_test_6_18 cd4, c(350) all fuzzy(art_6m) masspoints(adjust)

/* --------------------------------
  Table 5:  Local Randomization Neighborhood Sensitivity Diagnostic
----------------------------------*/
* In order to always use the same number of observations, keep only observations non-missing in covariates, Y, and D
keep if cd4!=. & art_6m !=. & visit_test_6_18 !=.
foreach var of global covariates {
	keep if `var'!=.
}

** 340--360
rdrandinf art_6m cd4, c(350) seed(5023) wl(340) wr(360)  ci(.05 -.5(.01).5) wmass
rdrandinf visit_test_6_18 cd4, c(350) seed(5023) wl(340) wr(360)  ci(.05 -.5(.01).5) wmass 
rdrandinf visit_test_6_18 cd4, c(350) seed(5023) wl(340) wr(360)  ci(.05 -2(.1)2) wmass fuzzy(art_6m tsls)

** 335--365
rdrandinf art_6m cd4, c(350) seed(5023) wl(335) wr(365)  ci(.05 -.5(.01).5) wmass
rdrandinf visit_test_6_18 cd4, c(350) seed(5023) wl(335) wr(365)  ci(.05 -.5(.01).5) wmass 
rdrandinf visit_test_6_18 cd4, c(350) seed(5023) wl(335) wr(365)  ci(.05 -2(.1)2) wmass fuzzy(art_6m tsls)

** 330--370
rdrandinf art_6m cd4, c(350) seed(5023) wl(330) wr(370)  ci(.05 -.5(.01).5) wmass
rdrandinf visit_test_6_18 cd4, c(350) seed(5023) wl(330) wr(370)  ci(.05 -.5(.01).5) wmass 
rdrandinf visit_test_6_18 cd4, c(350) seed(5023) wl(330) wr(370)  ci(.05 -2(.1)2) wmass fuzzy(art_6m tsls)

/* --------------------------------
  Table 6:  Continuity-Based Donut Hole Diagnostic 
----------------------------------*/
use "CKT_2023_SIM--ART.dta", clear
* In order to always use the same number of observations, keep only observations non-missing in covariates, Y, and D
keep if cd4!=. & art_6m !=. & visit_test_6_18 !=.

rdrobust visit_test_6_18 cd4 if (cd4!=350 & cd4!=351 & cd4!=349), c(350) all fuzzy(art_6m)
local h =  e(h_l)
local b =  e(b_l)
rdrobust visit_test_6_18 cd4 if (cd4!=350 & cd4!=351 & cd4!=349), c(350) all h(`h') b(`b')

/* --------------------------------
 Table 7: Continuity-Based Placebo Cutoffs Diagnostic
----------------------------------*/
use "CKT_2023_SIM--ART.dta", clear

** Set to 300, keep only those below 350 (assigned to treated)
rdrobust art_6m cd4 if cd4<350, c(300) all masspoints(adjust)
rdrobust visit_test_6_18 cd4 if cd4<350, c(300) all masspoints(adjust)


** Set to 400, keep only those above 350 (assigned to control)
rdrobust art_6m cd4 if cd4>=350, c(400) all masspoints(adjust)
rdrobust visit_test_6_18 cd4 if cd4>=350, c(400) all masspoints(adjust)

/* --------------------------------
 Table 8: Continuity-Based Fuzzy RD Estimates for Predetermined Covariates
----------------------------------*/
foreach var of global covariates {
	rdrobust `var' cd4, c(350) fuzzy(art_6m) bwselect(mserd)
}


