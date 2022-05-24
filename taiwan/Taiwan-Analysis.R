
###################--------------------------##############
# Replication files for "A Guide to Regression Discontinuity Designs in Medical Applications"
# Cost-sharing application
# by Matias D. Cattaneo, Luke Keele, and Rocio Titiunik
# last modified: May 2022
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
data <- read.dta("./taiwan/visit-data.dta")

####################################################################################
# Plots

## Figure 6: RD plot
rdplot(data$op_n_per, data$sev_day, c=51.5, p=2,
       y.label="Number of visits per 10,000", 
       x.label="Days from Third Birthday", 
       title="",
       ci=TRUE)

## Figure 7: Binning of daily observations be week above and below the cutoff of the third birthday.
data.zoom <- data %>% filter(sev_day > 44 & sev_day<60)                                 
ggplot(data.zoom, aes(x = sev_day, y = op_n_per)) +
		geom_point() +
		geom_vline(xintercept = 51) +
		theme_bw() +
		ylab("Number of visits per 10,000") +
		xlab("Weeks from third birthday")
	

## Diagnostics
## Baseline Covariates Balance --- Window Selection

rm(list=ls())
data <- read.dta("./taiwan/diagnostic-data.dta")

# Baseline Variables
vars <- c( "id_male", "h_inc_per", "tpe", 
            "b_year")

Z <- data[vars]

out = rdwinselect(data$sev_day, Z, c=51.5, seed = 50, reps = 1000, wmin=1, wstep=1, level= 0.15)


## Figure 8
## Minimum p-value plot
pvals <- c(out$results[1:10,1])
win1 <- seq(1:10)

df <- as.data.frame(cbind(pvals, win1))
df$window <- ordered(df$win1 , levels = c(1,2,3,4,5,6,7,8,9,10),
				labels = c("Window 1", "Window 2", "Window 3", "Window 4", "Window 5", "Window 6", "Window 7", "Window 8", "Window 9", "Window 10"))

ggplot(df, aes(x=window, y=pvals, group = 1)) + 
   geom_point() +
   ylab(" Minimum p-value") + 
   xlab("Window Length") +
   theme_bw() +
   theme(legend.title=element_blank())
  


## Table 11: Baseline Covariates Table
n_runs <- length(vars)
out2 <- matrix(NA, n_runs,4)

for(i in 1:n_runs){
  
  temp <- rdrandinf(Z[,i], data$sev_day, cutoff = 51.5, seed= 5023, wl=51, wr=52)
  out2[i,1] <- temp$sumstats[3,1]
  out2[i,2] <- temp$sumstats[3,2]
  out2[i,3] <- temp$obs.stat
  out2[i,4] <- temp$p.value
}	

temp$sumstats[2,1] + temp$sumstats[2,2]                         


var_names <- c("Share of male", "Household Income per Capita",
               "Share of children born in Taipei",
               "Birth year") 
               
               
colnames(out2) <- c("Mean Below", "Mean Above", "Diff. in Means", "p-value")
rownames(out2) <- var_names

# Table 11
round(out2, 3)
               
data <- read.dta("./taiwan/visit-data.dta")

# Outcome Analysis
# Local Randomization in selected window
## Table 12
lr.weekly <- rdrandinf(data$op_n_per, data$sev_day, 
                       cutoff=51.5, 
                       seed= 5023, 
                       wl=51, 
                       wr=52)

# Matches Continiuty-based results with larger window
lr.3weeks <- rdrandinf(data$op_n_per, data$sev_day, 
                       cutoff=51.5, 
                       seed= 5023, 
                       wl=49, 
                       wr=54)

# Continuity Based (reported in text)

cont.weekly <- rdrobust(data$op_n_per, data$sev_day, 
                        c=51.5, 
                        kernel="triangular", 
                        p=1, 
                        all=TRUE, 
                        masspoints="adjust")
summary(cont.weekly)

cont.daily <- rdrobust(data$op_n_per, data$age_diff, 
                        c=.5, 
                        kernel="triangular", 
                        p=1, 
                        all=TRUE, 
                        masspoints="adjust")
summary(cont.daily)



rm(list=ls())                                                                  
 ## Placebo Analysis  
data <- read.dta("./taiwan/placebo-data.dta")  

## Table 14                 
placebo1 <- rdrandinf(data$op_n_per, data$sev_day, 
                        cutoff=51.5, 
                        seed= 5023, 
                        wl=51, 
                        wr=52) 
                        
## One Month Window
placebo2 <- rdrandinf(data$op_n_per, data$sev_day,
                        cutoff=51.5, 
                        seed= 5023, 
                        wl=49, 
                        wr=54)    
                                                 
                                             