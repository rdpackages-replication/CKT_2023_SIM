/*###################--------------------------##############
# Replication files for "A Guide to Regression Discontinuity Designs in Medical Applications"
# Cost-sharing application
# by Matias D. Cattaneo, Luke Keele, and Rocio Titiunik
###################--------------------------##############*/

* To install necessary packages for RD analysis, execute the following commands 

*install.packages('rdrobust')
*install.packages('rdlocrand')
*install.packages('rddensity')
*install.packages('lpdensity')

* For more details, see https://rdpackages.github.io
clear all
use "CKT_2023_SIM--CostSharing.dta", clear

/* --------------------------------
 Figure 5(a): Scatterplot 
----------------------------------*/
scatter op_n_per week, ytitle("Number of visits per 10,000") xtitle("Weeks from Third Birthday") xline(0)

/* --------------------------------
 Figure 5(b): RD Plot 
----------------------------------*/
rdplot op_n_per week, c(0) p(1)  graph_options(ytitle("Number of visits per 10,000") xtitle("Weeks from Third Birthday"))

/* --------------------------------
 Table 9: Continuity-Based Methods 
----------------------------------*/
rdrobust op_n_per week, c(0) all

/* --------------------------------
 Window selection (in text)
----------------------------------*/
rdwinselect week sh_male h_inc_per sh_tpe b_year, c(0.5) seed(50) reps(1000) wstep(1) wmin(1) plot

/* --------------------------------
 Table 10: Distribution Pre-intervention covariates 
----------------------------------*/
foreach y in sh_male h_inc_per sh_tpe b_year {
		rdrandinf `y' week, c(-0.5) seed(5023) wl(-1) wr(0)
}

/* --------------------------------
 Table 11: Local Randomization Methods 
----------------------------------*/
rdrandinf op_n_per week, c(-0.5) seed(5023) wl(-1) wr(0)

/* --------------------------------
  Table 12: Local Randomization Methods Pre-intervention covariates
----------------------------------*/
rdrandinf pre_op_n_per week, c(-0.5) seed(5023) wl(-1) wr(0)