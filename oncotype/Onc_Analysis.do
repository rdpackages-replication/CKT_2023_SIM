/*###################--------------------------##############
# Replication files for "A Guide to Regression Discontinuity Designs in Medical Applications"
# Chemotherapy application
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

use "./oncotype/onc_rdd_clean.dta"

keep if onc_score >= 15 & onc_score <=30

** Density Test
rddensity onc_score, c(25.5) vce(jackknife) h(5 5) bino_wstep(1 1)

** Balance tests for baseline covariates
global covariates age white afam other hisp  grade1 grade2 grade3 grade_miss tumor_sz ///
                  lv_inv estro_rec progest_rec under50 surg_t size_miss

*Table 7
foreach var of global covariates {
    rdrandinf `var' onc_score, c(25.9) seed(5023) wl(25) wr(26)
}

*Table 8
foreach var of global covariates {
    rdrandinf `var' onc_score, c(25.9) seed(5023) wl(24) wr(27)
}

** Table 9 
rdrandinf chemo onc_score, c(25.9) seed(5023) wl(25) wr(26)
reg chemo above if onc_score>=25 & onc_score<=26
rdrandinf chemo onc_score, c(25.9) seed(5023) wl(24) wr(27)
reg chemo above if onc_score>=24 & onc_score<=27

** Table 10
rdrandinf cancer2 onc_score, c(25.5) seed(5023) wl(25) wr(26)  ci(.05 -.5(.01).5) wmass
rdrandinf cancer2 onc_score, c(25.5) seed(5023) wl(24) wr(27)  ci(.05 -.5(.01).5) wmass 
