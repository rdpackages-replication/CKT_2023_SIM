
###################--------------------------##############
# Replication files for "A Guide to Regression Discontinuity Designs in Medical Applications"
# Chemotherapy application
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
library(dplyr)
library(ggplot2)

# Make sure to set working directory to "replication" folder
data <- read.dta("./oncotype/onc_rdd_clean.dta", warn.missing.labels = FALSE)
data <- data %>% filter(onc_score >= 15, onc_score <=30)

predict.data <- data.frame(expand.grid(list(score = sort(unique(data$onc_score)))))
plot_data <- data %>% group_by(onc_score) %>% summarise(comp = mean(chemo, na.rm=TRUE), 
                                                       canc = mean(cancer2, na.rm=TRUE),
                                                       sd.ch = sd(chemo, na.rm=TRUE),
                                                       sd.ca = sd(cancer2, na.rm=TRUE),
                                                       n = length(chemo))

plot_data <- as.data.frame(cbind(predict.data, plot_data))


plot_data$se.chemo <- plot_data$sd.ch/sqrt(plot_data$n)
plot_data$se.cancer <- plot_data$sd.ca/sqrt(plot_data$n)
plot_data$ch.ci.1 <- plot_data$comp + qnorm(0.025, lower.tail=FALSE)*plot_data$se.chemo
plot_data$ch.ci.2 <- plot_data$comp - qnorm(0.025, lower.tail=FALSE)*plot_data$se.chemo
plot_data$ca.ci.1 <- plot_data$canc + qnorm(0.025, lower.tail=FALSE)*plot_data$se.cancer
plot_data$ca.ci.2 <- plot_data$canc - qnorm(0.025, lower.tail=FALSE)*plot_data$se.cancer

## Figure 1
ggplot(plot_data, aes(x=onc_score, y=comp)) +
      geom_line(aes(y=comp, color="orangered2"), show.legend=FALSE) +
      geom_point(aes(y=comp), show.legend=FALSE) +
      xlab("Oncotype Score") +
      ylab("Proprotion Receiving Chemotherapy") +
      theme(legend.position="none") +
      geom_vline(aes(xintercept=26), linetype =2) +
      geom_errorbar(aes(ymax = ch.ci.1, ymin = ch.ci.2), 
      position=position_dodge(0.05), width = .1) +
      theme_bw()


## Figure 2
ggplot(plot_data, aes(x=onc_score, y=canc)) +
  geom_line(aes(y=canc, color="orangered2"), show.legend=FALSE) +
  geom_point(aes(y=canc), show.legend=FALSE) +
  xlab("Oncotype Score") +
  ylab("Proportion with Cancer Recurrence") +
  theme(legend.position="none") +
  geom_vline(aes(xintercept=26), linetype =2) +
  geom_errorbar(aes(ymax = ca.ci.1, ymin = ca.ci.2), 
  position=position_dodge(0.05), width = .1) +
  theme_bw()

## Density Test, reported in text
# Check via Binomial Test
library(rddensity)
mtest <- rddensity(X = data$onc_score, vce="jackknife", c=25.5, binoWStep=1, q=2, h=c(5,5))

### Table of Output
cbind(mtest$bino$LeftN, mtest$bino$RightN, round(mtest$bino$pval, 3))


## Balance Tests for Baseline Covariates
library(rdlocrand)

vars <- c( "age", "white","afam","other","hisp", 
             "grade1","grade2", "grade3", "grade_miss",
             "tumor_sz","lv_inv", "estro_rec",
             "progest_rec", "under50", "surg_t", "size_miss")

Z <- data[vars]
n_runs <- length(vars)

##### Balance for Window 25 -- 26
out1 <- matrix(NA, n_runs,4)

for(i in 1:n_runs){
  temp <- rdrandinf(Z[,i], data$onc_score, cutoff = 25.9, seed= 5023, wl=25, wr=26)
  out1[i,1] <- temp$sumstats[3,1]
  out1[i,2] <- temp$sumstats[3,2]
  out1[i,3] <- temp$obs.stat
  out1[i,4] <- temp$p.value
}	

temp$sumstats[2,1] + temp$sumstats[2,2]

var_names <- c("Age", "White", "Af-American", "Other Race", "Hispanic", "Lo Grade Tumor",
               "Int. Grade Tumor", "High Grade Tumor",  "Grade Miss",             
               "Tumor Size", "Lymphovascular Invasion", "Estrogen Receptor", "Proges. Receptor",
               "Under 50", "Surgery Type 0/1", "Tumor Size Missing")

colnames(out1) <- c("Mean Below", "Mean Above", "Diff. in Means", "p-value")
rownames(out1) <- var_names

## Table 7
out1
temp$sumstats[2,1] + temp$sumstats[2,2]

##### Balance for Window 24 -- 27
out2 <- matrix(NA, n_runs,4)

for(i in 1:n_runs){
  temp <- rdrandinf(Z[,i], data$onc_score, cutoff = 25.9, seed= 5023, wl=24, wr=27)
  out2[i,1] <- temp$sumstats[3,1]
  out2[i,2] <- temp$sumstats[3,2]
  out2[i,3] <- temp$obs.stat
  out2[i,4] <- temp$p.value
}	

colnames(out2) <- c("Mean Below", "Mean Above", "Diff. in Means", "p-value")
rownames(out2) <- var_names

# Table 8
out2
temp$sumstats[2,1] + temp$sumstats[2,2]

##### Balance for Window 23 -- 28
out3 <- matrix(NA, n_runs,4)

for(i in 1:n_runs){
  temp <- rdrandinf(Z[,i], data$onc_score, cutoff = 24.9, seed= 5023, wl=23, wr=28)
  out3[i,1] <- temp$sumstats[3,1]
  out3[i,2] <- temp$sumstats[3,2]
  out3[i,3] <- temp$obs.stat
  out3[i,4] <- temp$p.value
}	

colnames(out3) <- c("Mean Below", "Mean Above", "Diff. in Means", "p-value")
rownames(out3) <- var_names

out3
temp$sumstats[2,1] + temp$sumstats[2,2]

##### Balance for Window 22 -- 29
out4 <- matrix(NA, n_runs,4)

#memory.limit(30300000)

for(i in 1:n_runs){
  temp <- rdrandinf(Z[,i], data$onc_score, cutoff = 24.9, seed= 5023, wl=22, wr=29)
  out4[i,1] <- temp$sumstats[3,1]
  out4[i,2] <- temp$sumstats[3,2]
  out4[i,3] <- temp$obs.stat
  out4[i,4] <- temp$p.value
}	

colnames(out4) <- c("Mean Below", "Mean Above", "Diff. in Means", "p-value")
rownames(out4) <- var_names

temp$sumstats[2,1] + temp$sumstats[2,2]

##### Balance for Window 21 -- 30
out5 <- matrix(NA, n_runs,4)

#memory.limit(30300000)

for(i in 1:n_runs){
  temp <- rdrandinf(Z[,i], data$onc_score, cutoff = 24.9, seed= 5023, wl=21, wr=30)
  out5[i,1] <- temp$sumstats[3,1]
  out5[i,2] <- temp$sumstats[3,2]
  out5[i,3] <- temp$obs.stat
  out5[i,4] <- temp$p.value
}	

pvals <- c(min(out1[,4]), min(out2[,4]), min(out3[,4]), min(out4[,4]), min(out5[,4]))
win1 <- seq(1:5)

df <- as.data.frame(cbind(pvals, win1))
df$window <- ordered(df$win1 , levels = c(1,2,3,4,5),
				labels = c("Window: 25-26", "Window 24-27", "Window 23-28", "Window 22-29", "Window 21-30"))

p <- ggplot(df, aes(x=window, y=pvals, group = 1)) + 
   geom_point() +
   ylab(" Minimum p-value") + 
   xlab("Window Length") +
   theme_bw() +
   theme(legend.title=element_blank())

## Figure 7
p

## Table 9
####### Local Neighborhood 1
ci_vec = c(0.05, seq(from = -1, to = 1, by = 0.01))

itt.chemo.w1 <- rdrandinf(data$chemo, data$onc_score, 
                        cutoff = 25.9, seed= 5023, wl=25, wr=26, ci=ci_vec, 
                        wmasspoints=TRUE)

####### Local Neighborhood 2
itt.chemo.w2 <- rdrandinf(data$chemo, data$onc_score, 
                        cutoff = 25.9, seed= 5023, wl=24, wr=27, ci=ci_vec, 
                        wmasspoints=TRUE)

### Large Sample Results
data.win1 <- data %>%  filter(win1==1)
n.w1 <- nrow(data.win1)
data.win2 <- data %>%  filter(win2==1)
n.w2 <- nrow(data.win2)
table(data.win2$onc_score)
fs.w1 <- lm(chemo ~ above, data=data.win1)
fs.w2 <- lm(chemo ~ above, data=data.win2)


col0 <- c(itt.chemo.w1$obs.stat, itt.chemo.w1$p.value, summary(fs.w1)$r.squared, summary(fs.w1)$fstatistic[1], 
          summary(fs.w1)$coef[2,4], length(summary(fs.w1)$resid))
          
col1 <- c(itt.chemo.w2$obs.stat, itt.chemo.w2$p.value, summary(fs.w2)$r.squared, summary(fs.w2)$fstatistic[1], 
          summary(fs.w2)$coef[2,4], length(summary(fs.w2)$resid))


tab3 <- cbind(col0, col1)
rownames(tab3) <- c("1st Stage Estimate", "Fisherian p-value", "R-squared", "F-statistic", "F-test p-value", "Sample Size")
tab3

## Table 10

####### Local Neighborhood 1
ci_vec = c(0.05, seq(from = -1, to = 1, by = .01))

itt.canc.w1 <- rdrandinf(data$cancer2, data$onc_score, 
                        cutoff = 25.9, seed= 5023, wl=25, wr=26, ci=ci_vec)

row.1.w1 <- c(itt.canc.w1$obs.stat, itt.canc.w1$ci[1], itt.canc.w1$ci[2], itt.canc.w1$p.value)

####### Local Neighborhood 2
itt.canc.w2 <- rdrandinf(data$cancer2, data$onc_score, 
                        cutoff = 25.9, seed= 5023, wl=24, wr=27, ci=ci_vec, 
                        wmasspoints=TRUE)

row.1.w2 <- c(itt.canc.w2$obs.stat, itt.canc.w2$ci[1], itt.canc.w2$ci[2], itt.canc.w2$p.value)

round(rbind(row.1.w1, row.1.w2), 3)


