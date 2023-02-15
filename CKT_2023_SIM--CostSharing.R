###################--------------------------##############
# Replication files for "A Guide to Regression Discontinuity Designs in Medical Applications"
# Cost-sharing application
# by Matias D. Cattaneo, Luke Keele, and Rocio Titiunik
###################--------------------------##############

# To install necessary packages for RD analysis, execute the following commands 

#install.packages('rdrobust')
#install.packages('rdlocrand')
#install.packages('rddensity')
#install.packages('lpdensity')

#For more details, see https://rdpackages.github.io

rm(list=ls())
library(foreign)
library(rdrobust)
library(rddensity)
library(ggplot2)
library(xtable)
library(dplyr)
library(rdlocrand)

# Make sure to set working directory to "replication" folder
data <- read.dta("CKT_2023_SIM--CostSharing.dta")

####################################################################################
# Plots

# --------------------------------#
# Figure 5(a): Scatterplot 
# --------------------------------#
#data.zoom <- data %>% filter(sev_day > 44 & sev_day<60)                                 
ggplot(data, aes(x = week, y = op_n_per)) +
		geom_point() +
		geom_vline(xintercept = 0) +
		theme_bw() +
		ylab("Number of visits per 10,000") +
		xlab("Weeks from Third Birthday")

ggsave("output/CS_Fig5a.pdf", width = 7, height = (13/16)*7, units = "in")

# --------------------------------#
# Figure 5(b): RD Plot 
# --------------------------------#
rdplot(data$op_n_per, data$week, c=0, p=1,
      y.label="Number of visits per 10,000", 
      x.label="Weeks from Third Birthday", 
      title="")

ggsave("output/CS_Fig5b.pdf", width = 7, height = (13/16)*7, units = "in")

# --------------------------------#
# Table 9: Continuity-Based Methods 
# --------------------------------#
cont <- rdrobust(data$op_n_per, data$week, c=0, all=TRUE)
summary(cont)
rr = function(x) round(x, digits=2)
row.1 <- c(rr(cont$coef[1]), paste("[", rr(cont$ci[3,1]), ",", rr(cont$ci[3,2]), "]", sep=""), rr(cont$bws[1]), rr(cont$N_h[1]), rr(cont$N_h[2]))
row.1 <- matrix(row.1, 1, 5)

rownms = c("Number of Doctor Visits per 10,000")
names(row.1) = rownms
print(xtable(row.1), include.rownames=TRUE)


# --------------------------------#
# Window selection (in text)
# --------------------------------#
vars <- c("sh_male", "h_inc_per", "sh_tpe", "b_year")

Z <- data[vars]

out = rdwinselect(data$week, Z, c=-0.5, seed = 50, reps = 1000, wmin=1, wstep=1, level= 0.15)

# --------------------------------#
# Table 10: Distribution Pre-intervention covariates 
# --------------------------------#
n_runs <- length(vars)
out2 <- matrix(NA, n_runs,4)

for(i in 1:n_runs){
  temp <- rdrandinf(Z[,i], data$week, cutoff = -0.5, seed= 5023, wl=-1, wr=0)
  out2[i,1] <- temp$sumstats[3,1]
  out2[i,2] <- temp$sumstats[3,2]
  out2[i,3] <- temp$obs.stat
  out2[i,4] <- temp$p.value
}	

temp$sumstats[2,1] + temp$sumstats[2,2]                         

var_names <- c("Share of male", "Household Income per Capita",
               "Share of children born in Taipei", "Birth year") 

colnames(out2) <- c("Mean Below", "Mean Above", "Diff. in Means", "p-value")
rownames(out2) <- var_names

round(out2, 3)

# --------------------------------#
# Table 11: Local Randomization Methods 
# --------------------------------#
lr <- rdrandinf(data$op_n_per, data$week, cutoff=-0.5, seed=5023, 
                       wl=-1, wr=0)

# --------------------------------#
# Table 12: Local Randomization Methods Pre-intervention covariates
# --------------------------------#
placebo <- rdrandinf(data$pre_op_n_per, data$week, cutoff=-0.5, seed=5023, 
                     wl=-1, wr=0) 

