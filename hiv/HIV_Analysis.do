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
clear all
use "./hiv/hiv_data.dta", clear

** Plots

** Figure 1a
rdplot art_6m cd4, c(350) binselect(qs) graph_options(yscale(r(0 1)) xscale(r(0 1000)))

** Figure 1b
rdplot visit_test_6_18 cd4, c(350) binselect(qs) graph_options(yscale(r(0 1)) xscale(r(0 1000)))

*####################################### Robustness Checks ####################################### 
* Figure 2a
histogram cd4, width(10) xline(350)


* Figure 2b (density test)
rddensity cd4, c(350) vce(jackknife)

*  RD Effects on Baseline Covariates
global covariates age1 age2  age3 age4 age5 age6 age7 age8 qtr1 qtr2 qtr3 qtr4 qtr5 qtr6 female clinic_a clinic_b clinic_c
** Table 1: Intention-to-treat RD effects on the covariates
foreach var of global covariates {
	rdrobust `var' cd4, c(350) bwselect(cerrd)
}


** Table 2: Fuzzy RD effects on the covariates
foreach var of global covariates {
	rdrobust `var' cd4, c(350) fuzzy(art_6m) bwselect(cerrd)
}

** Table 3
* Estimate RD Treatment Effect
* Local Linear w/ Quadratic - MSE Bandwdith - Bias Corrected CIs

* First-stage effect and Fuzzy effect
rdrobust visit_test_6_18 cd4, c(350) all fuzzy(art_6m)
local h =  e(h_l)
local b =  e(b_l)
* Intention-to-treat effect in same bandwdiths
rdrobust visit_test_6_18 cd4, c(350) all h(`h') b(`b')


** Table 4
** Placebo Cutoffs
** Set to 300, keep only those below 350 (assigned to treated)
rdrobust visit_test_6_18 cd4 if cd4<350, c(300) all 
rdrobust art_6m cd4 if cd4<350, c(300) all

** Set to 400, keep only those above 350 (assigned to control)
rdrobust visit_test_6_18 cd4 if cd4>=350, c(400) all 
rdrobust art_6m cd4 if cd4>=350, c(400) 

** Table 5
** Donut Analysis
rdrobust visit_test_6_18 cd4 if (cd4!=350 & cd4!=351 & cd4!=349), c(350) all fuzzy(art_6m)
local h =  e(h_l)
local b =  e(b_l)
rdrobust visit_test_6_18 cd4 if (cd4!=350 & cd4!=351 & cd4!=349), c(350) all h(`h') b(`b')

** Local Randomization
** Table 6
** Window Selection
rdwinselect cd4 age1 age2 age3 age4 age5 age6 age7 age8 qtr1 qtr2 ///
qtr3 qtr4 qtr5 qtr6 clinic_a clinic_b clinic_c, c(350) seed(5) reps(1000) wstep(1)


* In order to always use the same number of observations, keep only observations non-missing in covariates, Y, and D
keep if cd4!=. & art_6m !=. & visit_test_6_18 !=.
foreach var of global covariates {
	keep if `var'!=.
}
rdrandinf art_6m cd4, c(350) seed(5023) wl(346) wr(354)  ci(.05 -.5(.01).5) wmass
rdrandinf visit_test_6_18 cd4, c(350) seed(5023) wl(346) wr(354)  ci(.05 -.5(.01).5) wmass 
rdrandinf visit_test_6_18 cd4, c(350) seed(5023) wl(346) wr(354)  ci(.05 -2(.1)2) wmass fuzzy(art_6m tsls)


** Wider Windows
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
