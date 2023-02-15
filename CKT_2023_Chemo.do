/*###################--------------------------##############
# Replication files for "A Guide to Regression Discontinuity Designs in Medical Applications"
# Chemotherapy application
# by Matias D. Cattaneo, Luke Keele, and Rocio Titiunik
###################--------------------------##############*/

/* To install necessary packages for RD analysis, execute the following commands 

net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace
net install rdlocrand, from(https://raw.githubusercontent.com/rdpackages/rdlocrand/master/stata) replace
net install rddensity, from(https://raw.githubusercontent.com/rdpackages/rddensity/master/stata) replace
net install lpdensity, from(https://raw.githubusercontent.com/nppackages/lpdensity/master/stata) replace

For more details, see https://rdpackages.github.io
*/

clear all
use "CKT_2023_Chemo.dta", clear
keep if onc_score >= 15 & onc_score <=30

/* --------------------------------
 Density Test, reported in text
 Check via Binomial Test
----------------------------------*/
rddensity onc_score, c(25.5) vce(jackknife) h(5 5) bino_wstep(1 1) q(2)

/* --------------------------------
 Pre-intervention Covariate Diagnostics
----------------------------------*/
** Balance tests for baseline covariates
global covariates age white afam other hisp  grade1 grade2 grade3 grade_miss tumor_sz ///
                  lv_inv estro_rec progest_rec under50 surg_t size_miss
* Window [25,26]
foreach var of global covariates {
    rdrandinf `var' onc_score, c(25.5) seed(5023) wl(25) wr(26)
}

* Window [24,27]
foreach var of global covariates {
    rdrandinf `var' onc_score, c(25.5) seed(5023) wl(24) wr(27)
}

* Window [23,28]
foreach var of global covariates {
    rdrandinf `var' onc_score, c(25.5) seed(5023) wl(23) wr(28)
}

* Window [22,29]
foreach var of global covariates {
    rdrandinf `var' onc_score, c(25.5) seed(5023) wl(22) wr(29)
}

* Window [21,30]
foreach var of global covariates {
    rdrandinf `var' onc_score, c(25.5) seed(5023) wl(21) wr(30)
}

