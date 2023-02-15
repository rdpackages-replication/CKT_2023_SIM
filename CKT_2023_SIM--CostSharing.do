
/*###################--------------------------##############
# Replication files for "A Guide to Regression Discontinuity Designs in Medical Applications"
# Cost-sharing application
# by Matias D. Cattaneo, Luke Keele, and Rocio Titiunik
# last modified: May 2022
###################--------------------------##############*/

* To install necessary packages for RD analysis, execute the following commands 

*install.packages('rdrobust')
*install.packages('rdlocrand')
*install.packages('rddensity')
*install.packages('lpdensity')

* For more details, see https://rdpackages.github.io


clear all

use "./taiwan/visit-data.dta", clear

** Figure 6
rdplot op_n_per sev_day, c(51.5) p(2)

** Figure 7 
keep if sev_day > 44 & sev_day<60
scatter op_n_per sev_day, xline(51) ytitle("Number of visits per 10,000") xtitle("Weeks from third birthday")

use "./taiwan/diagnostic-data.dta",clear
* Baseline Covariates Balance --- Window Selection
rdwinselect sev_day id_male h_inc_per tpe ///
            b_year, c(51.5) seed(50) reps(1000) wstep(1) wmin(1) plot
			
** Table 11		
rdrandinf id_male sev_day, c(51.5) seed(5023) wl(51) wr(52)
rdrandinf h_inc_per sev_day, c(51.5) seed(5023) wl(51) wr(52)
rdrandinf tpe sev_day, c(51.5) seed(5023) wl(51) wr(52)
rdrandinf b_year sev_day, c(51.5) seed(5023) wl(51) wr(52)

use "./taiwan/visit-data.dta", replace
* Table 12
rdrandinf op_n_per sev_day, c(51.5) seed(5023) wl(51) wr(52) wmass

* Continuity-based, reported in text
rdrobust op_n_per sev_day, c(51.5) all
rdrobust op_n_per age_diff, c(0.5) all


* Placebo analysis
** Table 13
use "./taiwan/placebo-data.dta", replace

rdrandinf op_n_per sev_day, c(51.5) seed(5023) wl(51) wr(52) wmass
