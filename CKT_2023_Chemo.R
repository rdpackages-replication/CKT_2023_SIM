###################--------------------------##############
# Replication files for "A Guide to Regression Discontinuity Designs in Medical Applications"
# Chemotherapy application
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
library(dplyr)
library(ggplot2)
library(rdrobust)

data <- read.dta("CKT_2023_Chemo.dta", warn.missing.labels = FALSE)
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

# --------------------------------#
# Figure 6(a): RD Plot ITT
# --------------------------------#
ggplot(plot_data, aes(x=onc_score, y=comp)) +
      #geom_line(aes(y=comp, color="orangered2"), show.legend=FALSE) +
      geom_point(aes(y=comp), show.legend=FALSE) +
      xlab("Oncotype Score") +
      ylab("Proprotion Receiving Chemotherapy") +
      theme(legend.position="none") +
      geom_vline(aes(xintercept=26), linetype =2) +
      geom_errorbar(aes(ymax = ch.ci.1, ymin = ch.ci.2), 
      position=position_dodge(0.05), width = .1) +
      theme_bw() + scale_x_continuous(breaks=c(15, 20, 25, 26, 30)) 

ggsave("output/Chemo_Fig6a.pdf", width = 7, height = (13/16)*7, units = "in")


# --------------------------------#
# Figure 6(a): RD Plot ITT
# --------------------------------#
ggplot(plot_data, aes(x=onc_score, y=canc)) +
  #geom_line(aes(y=canc, color="orangered2"), show.legend=FALSE) +
  geom_point(aes(y=canc), show.legend=FALSE) +
  xlab("Oncotype Score") +
  ylab("Proportion with Cancer Recurrence") +
  theme(legend.position="none") +
  geom_vline(aes(xintercept=26), linetype =2) +
  geom_errorbar(aes(ymax = ca.ci.1, ymin = ca.ci.2), 
  position=position_dodge(0.05), width = .1) +
  theme_bw() + scale_x_continuous(breaks=c(15, 20, 25, 26, 30))  

ggsave("output/Chemo_Fig6b.pdf", width = 7, height = (13/16)*7, units = "in")


# --------------------------------#
# Density Test, reported in text
# Check via Binomial Test
# --------------------------------#
library(rddensity)
mtest <- rddensity(X = data$onc_score, vce="jackknife", c=25.5, binoWStep=1, q=2, h=c(5,5))
summary(mtest)
## Table of Output
cbind(mtest$bino$LeftN, mtest$bino$RightN, round(mtest$bino$pval, 3))

# --------------------------------#
# Pre-intervention Covariate Diagnostics
# --------------------------------#
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
  temp <- rdrandinf(Z[,i], data$onc_score, cutoff = 25.5, seed= 5023, wl=25, wr=26)
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

out1
temp$sumstats[2,1] + temp$sumstats[2,2]

##### Balance for Window 24 -- 27
out2 <- matrix(NA, n_runs,4)

for(i in 1:n_runs){
  temp <- rdrandinf(Z[,i], data$onc_score, cutoff = 25.5, seed= 5023, wl=24, wr=27)
  out2[i,1] <- temp$sumstats[3,1]
  out2[i,2] <- temp$sumstats[3,2]
  out2[i,3] <- temp$obs.stat
  out2[i,4] <- temp$p.value
}	

colnames(out2) <- c("Mean Below", "Mean Above", "Diff. in Means", "p-value")
rownames(out2) <- var_names

out2
temp$sumstats[2,1] + temp$sumstats[2,2]

##### Balance for Window 23 -- 28
out3 <- matrix(NA, n_runs,4)

for(i in 1:n_runs){
  temp <- rdrandinf(Z[,i], data$onc_score, cutoff = 25.5, seed= 5023, wl=23, wr=28)
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
  temp <- rdrandinf(Z[,i], data$onc_score, cutoff = 25.5, seed= 5023, wl=22, wr=29)
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

for(i in 1:n_runs){
  temp <- rdrandinf(Z[,i], data$onc_score, cutoff = 25.5, seed= 5023, wl=21, wr=30)
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

ggplot(df, aes(x=window, y=pvals, group = 1)) + 
   geom_point() +
   ylab(" Minimum p-value") + 
   xlab("Window Length") +
   theme_bw() +
   theme(legend.title=element_blank())


