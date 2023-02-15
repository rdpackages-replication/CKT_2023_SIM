###################--------------------------##############
# Replication files for "A Guide to Regression Discontinuity Designs in Medical Applications"
# ART application
# by Matias D. Cattaneo, Luke Keele, and Rocio Titiunik
###################--------------------------##############

# To install necessary packages for RD analysis, execute the following commands 

#install.packages('rdrobust')
#install.packages('rdlocrand')
#install.packages('rddensity')
#install.packages('lpdensity')

#For more details, see https://rdpackages.github.io

rm(list=ls())
library(dplyr)
library(foreign)
library(rdrobust)
library(rddensity)
library(rdlocrand)
library(ggplot2)
library(xtable)

# Make sure to set working directory to "replication" folder
rr = function(x) round(x, digits=2)
data <- read.dta("CKT_2023_SIM--ART.dta")
head(data)

# Score: cd4 c = 350
# X: age1-age8, qtr1-qtr6, clinic_a, clinic_b, clinic_c, female
# Y: visit_test_6_18
# D: art_6m

# --------------------------------#
# Figure 1(a): Score Histogram
# --------------------------------#
ggplot(data, aes(x=cd4)) + 
  geom_histogram(binwidth = 10, color="black", fill="lightblue") +
  xlim(0,1000) +
  geom_vline(aes(xintercept = 350)) +
  labs(x="CD4 Count", y="Counts") + theme_bw() +
  scale_x_continuous(breaks=c(0, 250, 350, 500, 750, 1000))

ggsave("output/ART_Fig1a.pdf", width = 7, height = (13/16)*7, units = "in")


# --------------------------------#
# Figure 1(b): Outcome Plot
# --------------------------------#
p1 <- rdplot(data$visit_test_6_18, data$cd4, c=350, 
             y.lim=c(0,1),
             x.lim=c(0,1000), 
             y.label="Proportion Retained in Treatment", 
             x.label="CD4 Count", 
             title="", p=3)

p1$rdplot + scale_x_continuous(breaks=c(0, 250, 350, 500, 750, 1000))                                     

ggsave("output/ART_Fig1b.pdf", width = 7, height = (13/16)*7, units = "in")


# --------------------------------#
#  Figure 3(a): Treatment Assignment
# --------------------------------#
p3 <- rdplot(data$below, data$cd4, c=350, 
             y.lim=c(0,1),
             x.lim=c(0,1000),
             y.label="Assigned to ART Initiation", 
             x.label="CD4 Count", 
             title="", 
             binselect = "es",
             p=0, 
             nbins=c(10,100))

p3$rdplot + 
  geom_point(aes(x=458, y=0))+
  geom_point(aes(x=566, y=0))+
  geom_point(aes(x=674, y=0))+
  geom_point(aes(x=782, y=0)) +
  geom_point(aes(x=890, y=0)) +
  geom_point(aes(x=76.5, y=1)) +
  geom_point(aes(x=251., y=1)) +
  scale_x_continuous(breaks=c(0, 350, 500, 750, 1000)) 

ggsave("output/ART_Fig3a.pdf", width = 7, height = (13/16)*7, units = "in")

# --------------------------------#
#  Figure 3(b): Treatment Take-up
# --------------------------------#
p2 <- rdplot(data$art_6m, data$cd4, c=350, 
             y.lim=c(0,1),
             x.lim=c(0,1000), 
             y.label="Proportion with ART Initiation", 
             x.label="CD4 Count", 
             title="", p=3)

p2$rdplot + scale_x_continuous(breaks=c(0, 250, 350, 500, 750, 1000)) 

ggsave("output/ART_Fig3b.pdf", width = 7, height = (13/16)*7, units = "in")


# --------------------------------#
# Table 1: Continuity-based RD Estimates with Robust Bias Corrected Inference
# --------------------------------#

## To always use same number of observations, keep only observations non-missing in Y, and D
dd = data[complete.cases(data[,c("visit_test_6_18", "art_6m", "cd4")]),]

# Run fuzzy analysis first and save bandwidths to report first-stage and reduced form within same bandwidth
fzy <- rdrobust(dd$visit_test_6_18, dd$cd4, c=350, all=TRUE, fuzzy=dd$art_6m)
row.3 <- c(rr(fzy$coef[1]), paste("[", rr(fzy$ci[3,1]), ",", rr(fzy$ci[3,2]), "]", sep=""), rr(fzy$bws[1]), rr(fzy$N_h[1]), rr(fzy$N_h[2]))
h = fzy$bws[1]
b = fzy$bws[2]

# First Stage: effect of assignment on ART Initiation
itt1 <- rdrobust(dd$art_6m, dd$cd4, c=350, all=TRUE, h=h, b=b)
row.1 <- c(rr(itt1$coef[1]), paste("[", rr(itt1$ci[3,1]), ",", rr(itt1$ci[3,2]), "]", sep=""), rr(itt1$bws[1]), rr(itt1$N_h[1]), rr(itt1$N_h[2]))

# Reduced form: effect of assignment on program retention
itt2 <- rdrobust(dd$visit_test_6_18, dd$cd4, c=350, all=TRUE, h=h, b=b)
row.2 <-c(rr(itt2$coef[1]), paste("[", rr(itt2$ci[3,1]), ",", rr(itt2$ci[3,2]), "]", sep=""), rr(itt2$bws[1]), rr(itt2$N_h[1]), rr(itt2$N_h[2]))

# bias bandwidths
itt1$bws[2]
itt2$bws[2]
fzy$bws[2]

rownms = c("ITT Effect of ART Assignment on ART Initiation", 
           "ITT Effect of ART Assignment on Program Retention",
           "Fuzzy Effect of ART Initiation on Program Retention")

tab = rbind(row.1, row.2, row.3)
rownames(tab) = rownms
print(xtable(tab), include.rownames=TRUE)

# save bandwidths for falsification test below
bws_original = fzy$bws

# --------------------------------#
# Table 2: Local Randomization RD Estimates and Inference
# --------------------------------#
## Window Selection (not shown in paper)
vars <- c( "age1", "age2", "age3", "age4",
           "age5", "age6", "age7", "age8",
           "qtr1", "qtr2", "qtr3", "qtr4",
           "qtr5", "qtr6", "clinic_a",
           "clinic_b", "clinic_c")

X <- data[vars]

out = rdwinselect(data$cd4, X, c=350, seed = 50, reps = 1000, wmasspoints=FALSE, wstep=1)

## In order to always use the same number of observations, keep only observations non-missing in covariates, Y, and D
dd = data[complete.cases(data[,c("age1", "age2", "age3", "age4",
                                 "age5", "age6", "age7", "age8",
                                 "qtr1", "qtr2", "qtr3", "qtr4",
                                 "qtr5", "qtr6", "female", "clinic_a",
                                 "clinic_b", "clinic_c",
                                 "visit_test_6_18", "art_6m", "cd4")]),]
#### ITT - D Analysis
ci_vec = c(0.05, seq(from = -.5, to = .5, by = 0.01))
itt.d <- rdrandinf(dd$art_6m, dd$cd4, 
                   cutoff = 350, seed= 5023, wl=346, wr=354, ci=ci_vec, 
                   wmasspoints=TRUE)                                             
row.d <- c(rr(itt.d$obs.stat), paste("[", itt.d$ci[1], " , ", rr(itt.d$ci[2]), "]", sep=""),  itt.d$sumstats[2,1], itt.d$sumstats[2,2],
           itt.d$p.value, itt.d$asy.pvalue)


itt.y <- rdrandinf(dd$visit_test_6_18, dd$cd4, 
                   cutoff = 350, seed= 5023, wl=346, wr=354, ci=ci_vec, 
                   wmasspoints=TRUE)                       
row.y <- c(rr(itt.y$obs.stat),  paste("[", rr(itt.y$ci[1]), " , ", rr(itt.y$ci[2]), "]", sep=""), itt.y$sumstats[2,1], itt.y$sumstats[2,2],
           itt.y$p.value, itt.y$asy.pvalue)

### Fuzzy RD Results
ci_vec = c(0.05, seq(from = -1, to = 1, by = 0.01))
iv.w1 <- rdrandinf(dd$visit_test_6_18, dd$cd4, fuzzy=c(fuzzy.tr=dd$art_6m, fuzzy.stat="tsls"), 
                   cutoff = 350, seed= 5023, wl=346, wr=354, ci=ci_vec)
row1.iv.w1 <- c(rr(iv.w1$obs.stat),  paste("[", rr(iv.w1$ci[1]), " , ",rr(iv.w1$ci[2]), "]", sep=""), iv.w1$sumstats[2,1], iv.w1$sumstats[2,2], 
                iv.w1$asy.pvalue)

# print full table
row.1.tab <- row.d[-c(5:6)]
row.2.tab <- row.y[-c(5:6)]
row.3.tab <- row1.iv.w1[-5]

tab = rbind(row.1.tab, row.2.tab, row.3.tab)
rownames(tab) = rownms
print(xtable(tab), include.rownames=TRUE)

# --------------------------------#
#  Figure 4(a): Histogram near cutoff: score in [345, 355]
# --------------------------------#
data.p <- data %>% filter(between(cd4, 345, 355))

ggplot(data.p, aes(x=cd4)) + 
  geom_histogram(color="black", fill="lightblue", stat="count") +
  #xlim(0,1000) 
  geom_vline(aes(xintercept = 350)) +
  labs(x="CD4 Count", y="Counts") + theme_bw()

ggsave("output/ART_Fig4a.pdf", width = 7, height = (13/16)*7, units = "in")


# --------------------------------#
#  Figure 4(b): Histogram and local polynomial density estimate for score 
# --------------------------------#
mtest <- rddensity(X = data$cd4, vce="jackknife", c=350)
summary(mtest)

rdplotdensity(mtest,X = data$cd4)$Estplot + 
  geom_vline(xintercept=350, linetype="dashed") +
  scale_x_continuous(breaks=c(0, 200, 350, 400, 600)) +
  labs(x="CD4 Count", y="Density")

ggsave("output/ART_Fig4b.pdf", width = 7, height = (13/16)*7, units = "in")


# --------------------------------#
# Table 3: Continuity-Based ITT RD Estimates for Predetermined Covariates
# --------------------------------#
# Baseline Variables
vars <- c( "age1", "age2", "age3", "age4",
           "age5", "age6", "age7", "age8",
           "qtr1", "qtr2", "qtr3", "qtr4",
           "qtr5", "qtr6", "female", "clinic_a",
           "clinic_b", "clinic_c")

Z <- data[vars]

mean(vars %in% names(data))

n_runs <- length(vars)
out1 <- matrix(NA, n_runs, 7)

for(i in 1:n_runs){
  
  temp <- rdrobust(Z[,i], data$cd4, c=350, bwselect="mserd", masspoints="adjust")
  out1[i,1] <- temp$tau_cl[1]
  out1[i,2] <- temp$tau_cl[2]
  out1[i,3] <- temp$coef[1]
  out1[i,4] <- temp$pv[3]
  out1[i,5] <- temp$bws[1]
  out1[i,6] <- temp$N_h[1] 
  out1[i,7] <- temp$N_h[2]
}	

var_names <- c("Age 0-18", "Age 18 -25", "Age 25-30", "Age 30-35", "Age 35-40",
               "Age 40-45", "Age 45-55", "Age 55+", "2011 Qtr3", "2011 Qtr4",
               "2012 Qtr1", "2012 Qtr2", "2012 Qtr3", "2012 Qtr4", "Female",
               "Clinic A", "Clinic B", "Clinic C")                            

colnames(out1) <- c("Mean Below", "Mean Above", "Diff. in Means", "Robust p-value", "MSE-Optimal Bandwidth",  "N^-_h","N^+_h")
rownames(out1) <- var_names

round(out1, 2)
xtable(out1, display = c("f","f","f", "f", "f", "f", "d", "d"))


# --------------------------------#
# Table 4: Continuity-Based Bandwidth Sensitivity Diagnostic
# --------------------------------#

## To always use same number of observations, keep only observations non-missing in Y, and D
dd = data[complete.cases(data[,c("visit_test_6_18", "art_6m", "cd4")]),]

## Expanded Bandwidths
h = bws_original[1] + 10
b = bws_original[2] + 10

# Run fuzzy analysis first and save bandwidths to report first-stage and reduced form within same bandwidth
fzy <- rdrobust(dd$visit_test_6_18, dd$cd4, c=350, all=TRUE, fuzzy=dd$art_6m, h=h, b=b)
row.3 <- c(rr(fzy$coef[1]), paste("[", rr(fzy$ci[3,1]), ",", rr(fzy$ci[3,2]), "]", sep=""), rr(fzy$bws[1]), rr(fzy$N_h[1]), rr(fzy$N_h[2]))


# First Stage: effect of assignment on ART Initiation
itt1 <- rdrobust(dd$art_6m, dd$cd4, c=350, all=TRUE, h=h, b=b)
row.1 <- c(rr(itt1$coef[1]), paste("[", rr(itt1$ci[3,1]), ",", rr(itt1$ci[3,2]), "]", sep=""), rr(itt1$bws[1]), rr(itt1$N_h[1]), rr(itt1$N_h[2]))

# Reduced form: effect of assignment on program retention
itt2 <- rdrobust(dd$visit_test_6_18, dd$cd4, c=350, all=TRUE, h=h, b=b)
row.2 <-c(rr(itt2$coef[1]), paste("[", rr(itt2$ci[3,1]), ",", rr(itt2$ci[3,2]), "]", sep=""), rr(itt2$bws[1]), rr(itt2$N_h[1]), rr(itt2$N_h[2]))

rownms = c("ITT Effect of ART Assignment on ART Initiation", 
           "ITT Effect of ART Assignment on Program Retention",
           "Fuzzy Effect of ART Initiation on Program Retention")

tab = rbind(row.1, row.2, row.3)
rownames(tab) = rownms
print(xtable(tab), include.rownames=TRUE)

h = bws_original[1] + 35
b = bws_original[2] + 35

# Run fuzzy analysis first and save bandwidths to report first-stage and reduced form within same bandwidth
fzy <- rdrobust(dd$visit_test_6_18, dd$cd4, c=350, all=TRUE, fuzzy=dd$art_6m, h=h, b=b)
row.3 <- c(rr(fzy$coef[1]), paste("[", rr(fzy$ci[3,1]), ",", rr(fzy$ci[3,2]), "]", sep=""), rr(fzy$bws[1]), rr(fzy$N_h[1]), rr(fzy$N_h[2]))


# First Stage: effect of assignment on ART Initiation
itt1 <- rdrobust(dd$art_6m, dd$cd4, c=350, all=TRUE, h=h, b=b)
row.1 <- c(rr(itt1$coef[1]), paste("[", rr(itt1$ci[3,1]), ",", rr(itt1$ci[3,2]), "]", sep=""), rr(itt1$bws[1]), rr(itt1$N_h[1]), rr(itt1$N_h[2]))

# Reduced form: effect of assignment on program retention
itt2 <- rdrobust(dd$visit_test_6_18, dd$cd4, c=350, all=TRUE, h=h, b=b)
row.2 <-c(rr(itt2$coef[1]), paste("[", rr(itt2$ci[3,1]), ",", rr(itt2$ci[3,2]), "]", sep=""), rr(itt2$bws[1]), rr(itt2$N_h[1]), rr(itt2$N_h[2]))

rownms = c("ITT Effect of ART Assignment on ART Initiation", 
           "ITT Effect of ART Assignment on Program Retention",
           "Fuzzy Effect of ART Initiation on Program Retention")

tab = rbind(row.1, row.2, row.3)
rownames(tab) = rownms
print(xtable(tab), include.rownames=TRUE)


# Check 
itt2$coef[1]/itt1$coef[1]
fzy$coef[1]


# --------------------------------#
# Table 5: Local Randomization Neighborhood Sensitivity Diagnostic
# --------------------------------#

## In order to always use the same number of observations, keep only observations non-missing in covariates, Y, and D
dd = data[complete.cases(data[,c("age1", "age2", "age3", "age4",
                                 "age5", "age6", "age7", "age8",
                                 "qtr1", "qtr2", "qtr3", "qtr4",
                                 "qtr5", "qtr6", "female", "clinic_a",
                                 "clinic_b", "clinic_c",
                                 "visit_test_6_18", "art_6m", "cd4")]),]

# Window 1
l.win <- 340
r.win <- 360

#### ITT - D Analysis
ci_vec = c(0.05, seq(from = -.5, to = .5, by = 0.01))

itt.d <- rdrandinf(dd$art_6m, dd$cd4, 
                   cutoff = 350, seed= 5023, wl=l.win, wr=r.win, ci=ci_vec, 
                   wmasspoints=TRUE)                                             
row.d <- c(rr(itt.d$obs.stat),  paste("[", rr(itt.d$ci[1]), " , ",rr(itt.d$ci[2]), "]", sep=""), itt.d$sumstats[2,1], itt.d$sumstats[2,2],  itt.d$p.value, itt.d$asy.pvalue)


itt.y <- rdrandinf(dd$visit_test_6_18, dd$cd4, 
                   cutoff = 350, seed= 5023, wl=l.win, wr=r.win, ci=ci_vec, 
                   wmasspoints=TRUE)                       
row.y <- c(rr(itt.y$obs.stat),  paste("[", rr(itt.y$ci[1]), " , ",rr(itt.y$ci[2]), "]", sep=""), itt.y$sumstats[2,1], itt.y$sumstats[2,2], itt.y$p.value, itt.y$asy.pvalue)

### Fuzzy RD Results

ci_vec = c(0.05, seq(from = -1, to = 1, by = 0.01))                
iv.w2 <- rdrandinf(dd$visit_test_6_18, dd$cd4, fuzzy=c(fuzzy.tr=dd$art_6m, fuzzy.stat="tsls"), 
                   cutoff = 350, seed= 5023, wl=l.win, wr=r.win, ci=ci_vec)                                                  
row1.iv.w2 <- c(rr(iv.w2$obs.stat),  paste("[", rr(iv.w2$ci[1]), " , ",rr(iv.w2$ci[2]), "]", sep=""), iv.w2$sumstats[2,1], iv.w2$sumstats[2,2], iv.w2$asy.pvalue)

### Print Results
row.1.tab <- row.d[-c(5:6)]
row.2.tab <- row.y[-c(5:6)]
row.3.tab <- row1.iv.w2[-5]

tab =rbind(row.1.tab, row.2.tab, row.3.tab)
rownames(tab) = rownms
print(xtable(tab), include.rownames=TRUE)

# Window 2
l.win <- 335
r.win <- 365

#### ITT - D Analysis
ci_vec = c(0.05, seq(from = -.5, to = .5, by = 0.01))

itt.d <- rdrandinf(dd$art_6m, dd$cd4, 
                   cutoff = 350, seed= 5023, wl=l.win, wr=r.win, ci=ci_vec, 
                   wmasspoints=TRUE)                                             
row.d <- c(rr(itt.d$obs.stat),  paste("[", rr(itt.d$ci[1]), " , ",rr(itt.d$ci[2]), "]", sep=""), itt.d$sumstats[2,1], itt.d$sumstats[2,2],  itt.d$p.value, itt.d$asy.pvalue)

itt.y <- rdrandinf(dd$visit_test_6_18, dd$cd4, 
                   cutoff = 350, seed= 5023, wl=l.win, wr=r.win, ci=ci_vec, 
                   wmasspoints=TRUE)                       
row.y <- c(rr(itt.y$obs.stat), paste("[", rr(itt.y$ci[1]), " , ",rr(itt.y$ci[2]), "]", sep=""), itt.y$sumstats[2,1], itt.y$sumstats[2,2],itt.y$p.value, itt.y$asy.pvalue)

### Fuzzy RD Results

ci_vec = c(0.05, seq(from = -1, to = 1, by = 0.01))

iv.w2 <- rdrandinf(dd$visit_test_6_18, dd$cd4, fuzzy=c(fuzzy.tr=dd$art_6m, fuzzy.stat="tsls"), 
                   cutoff = 350, seed= 5023, wl=l.win, wr=r.win, ci=ci_vec)                                                  
row1.iv.w2 <- c(rr(iv.w2$obs.stat), paste("[", rr(iv.w2$ci[1]), " , ",rr(iv.w2$ci[2]), "]", sep=""), iv.w2$sumstats[2,1], iv.w2$sumstats[2,2],iv.w2$asy.pvalue)

### Print Results
row.1.tab <- row.d[-c(5:6)]
row.2.tab <- row.y[-c(5:6)]
row.3.tab <- row1.iv.w2[-5]

tab =rbind(row.1.tab, row.2.tab, row.3.tab)
rownames(tab) = rownms
print(xtable(tab), include.rownames=TRUE)

# Window 3
l.win <- 330
r.win <- 370

#### ITT - D Analysis
ci_vec = c(0.05, seq(from = -.5, to = .5, by = 0.01))

itt.d <- rdrandinf(dd$art_6m, dd$cd4, 
                   cutoff = 350, seed= 5023, wl=l.win, wr=r.win, ci=ci_vec, 
                   wmasspoints=TRUE)                                             
row.d <- c(rr(itt.d$obs.stat), paste("[", rr(itt.d$ci[1]), " , ",rr(itt.d$ci[2]), "]", sep=""), itt.d$sumstats[2,1], itt.d$sumstats[2,2],  itt.d$p.value, itt.d$asy.pvalue)


itt.y <- rdrandinf(dd$visit_test_6_18, dd$cd4, 
                   cutoff = 350, seed= 5023, wl=l.win, wr=r.win, ci=ci_vec, 
                   wmasspoints=TRUE)                       
row.y <- c(rr(itt.y$obs.stat), paste("[", rr(itt.y$ci[1]), " , ",rr(itt.y$ci[2]), "]", sep=""), itt.y$sumstats[2,1], itt.y$sumstats[2,2], itt.y$p.value, itt.y$asy.pvalue)

### Fuzzy RD Results

ci_vec = c(0.05, seq(from = -1, to = 1, by = 0.01))


iv.w2 <- rdrandinf(dd$visit_test_6_18, dd$cd4, fuzzy=c(fuzzy.tr=dd$art_6m, fuzzy.stat="tsls"), 
                   cutoff = 350, seed= 5023, wl=l.win, wr=r.win, ci=ci_vec)                                                  
row1.iv.w2 <- c(rr(iv.w2$obs.stat), paste("[", rr(iv.w2$ci[1]), " , ",rr(iv.w2$ci[2]), "]", sep=""), iv.w2$sumstats[2,1], iv.w2$sumstats[2,2], iv.w2$asy.pvalue)

### Print Results
row.1.tab <- row.d[-c(5:6)]
row.2.tab <- row.y[-c(5:6)]
row.3.tab <- row1.iv.w2[-5]

tab =rbind(row.1.tab, row.2.tab, row.3.tab)
rownames(tab) = rownms
print(xtable(tab), include.rownames=TRUE)

# --------------------------------#
# Table 6:  Continuity-Based Donut Hole Diagnostic 
# --------------------------------#
## To always use same number of observations, keep only observations non-missing in Y, and D
dd = data[complete.cases(data[,c("visit_test_6_18", "art_6m", "cd4")]),]

data.donut <- dd %>% filter(!((cd4==350) | (cd4==351) | (cd4==349)))

# Run fuzzy analysis first and save bandwidths to report first-stage and reduced form within same bandwidth
fzy <- rdrobust(data.donut$visit_test_6_18, data.donut$cd4, c=350, all=TRUE, fuzzy=data.donut$art_6m)
row.3 <- c(rr(fzy$coef[1]), paste("[", rr(fzy$ci[3,1]), ",", rr(fzy$ci[3,2]), "]", sep=""), rr(fzy$bws[1]), rr(fzy$N_h[1]), rr(fzy$N_h[2]))
h = fzy$bws[1]
b = fzy$bws[2]

itt1 <- rdrobust(data.donut$art_6m, data.donut$cd4, c=350, all=TRUE, h=h, b=b)
row.1 <-c(rr(itt1$coef[1]),  paste("[", rr(itt1$ci[3,1]), ",", rr(itt1$ci[3,2]), "]", sep=""),rr(itt1$bws[1]), rr(itt1$N_h[1]), rr(itt1$N_h[2]))

itt2 <- rdrobust(data.donut$visit_test_6_18, data.donut$cd4, c=350, all=TRUE, h=h, b=b)
row.2 <- c(rr(itt2$coef[1]), paste("[", rr(itt2$ci[3,1]), ",", rr(itt2$ci[3,2]), "]", sep=""), rr(itt2$bws[1]), rr(itt2$N_h[1]), rr(itt2$N_h[2]))

# bias bandwidths
itt1$bws[2]
itt2$bws[2]
fzy$bws[2]

tab = rbind(row.1, row.2, row.3)
rownames(tab) = rownms
print(xtable(tab), include.rownames=TRUE)

# --------------------------------#
# Table 7: Continuity-Based Placebo Cutoffs Diagnostic
# --------------------------------#
# Set to 300, first keep only those below 350 (assigned to treated)
datas = data[data$cd4 < 350,]
itt1 <- rdrobust(datas$art_6m, datas$cd4, c=300, all=TRUE)
row.1 <- c(rr(itt1$coef[1]),  paste("[", rr(itt1$ci[3,1]), ",", rr(itt1$ci[3,2]), "]", sep=""),rr(itt1$bws[1]), rr(itt1$N_h[1]), rr(itt1$N_h[2]))

itt2 <- rdrobust(datas$visit_test_6_18, datas$cd4, c=300, all=TRUE)
row.2 <- c(rr(itt2$coef[1]), paste("[", rr(itt2$ci[3,1]), ",", rr(itt2$ci[3,2]), "]", sep=""), rr(itt2$bws[1]), rr(itt2$N_h[1]), rr(itt2$N_h[2]))

# bias bandwidths
itt1$bws[2]
itt2$bws[2]

tab = rbind(row.1, row.2)
rownames(tab) = rownms[1:2]
print(xtable(tab), include.rownames=TRUE)

# Placebo Cutoffs  
# Set to 400, first keep only those above 350 (assigned to control)
datas = data[data$cd4 >= 350,]
itt1 <- rdrobust(datas$art_6m, datas$cd4, c=400, all=TRUE)
row.1 <- c(rr(itt1$coef[1]),  paste("[", rr(itt1$ci[3,1]), ",", rr(itt1$ci[3,2]), "]", sep=""),rr(itt1$bws[1]), rr(itt1$N_h[1]), rr(itt1$N_h[2]))

itt2 <- rdrobust(datas$visit_test_6_18, datas$cd4, c=400, all=TRUE)
row.2 <- c(rr(itt2$coef[1]), paste("[", rr(itt2$ci[3,1]), ",", rr(itt2$ci[3,2]), "]", sep=""), rr(itt2$bws[1]), rr(itt2$N_h[1]), rr(itt2$N_h[2]))

# bias bandwidths
itt1$bws[2]
itt2$bws[2]

tab = rbind(row.1, row.2)
rownames(tab) = rownms[1:2]
print(xtable(tab), include.rownames=TRUE)


# --------------------------------#
# Table 8: Continuity-Based Fuzzy RD Estimates for Predetermined Covariates
# --------------------------------#
out2 <- matrix(NA, n_runs, 5)

for(i in 1:n_runs){
  
  temp1 <- rdrobust(Z[,i], data$cd4, c=350, all=TRUE, fuzzy=data$art_6m)
  out2[i,1] <- temp1$coef[1]
  out2[i,2] <- temp1$pv[3]
  out2[i,3] <- temp1$bws[1]
  out2[i,4] <- temp1$N_h[1] 
  out2[i,5] <- temp1$N_h[2]
}	

colnames(out2) <- c("tau_frd", "Robust p-value", "MSE-Optimal Bandwidth", "N^+_h","N^+_h")
rownames(out2) <- var_names

round(out2, 2)
xtable(out2, display = c("f", "f", "f", "f", "d", "d"))

summary(rdrobust(data$art_6m, data$cd4, c=350, all=TRUE))


# --------------------------------#
# Additional results reported in text
# --------------------------------#

#### Weak IV Test
bw.out <- rdbwselect(data$art_6m, data$cd4, c=350, all=TRUE)
summary(bw.out)

data.weak.iv <- data %>% filter(cd4 > 350-bw.out$bws[1,1] & cd4 < 350+bw.out$bws[1,1]) 
nrow(data.weak.iv)

mod1 <- lm(art_6m ~ 1, data=data.weak.iv)            
mod2 <- lm(art_6m ~ below, data=data.weak.iv)

anova(mod1, mod2)


