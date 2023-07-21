library(caret)
library(quantreg)
library(mnormt)
library(ks)
library(evd)
library(rdist)
library(evmix)
library(evgam)
library(EXRQ)
library(gbex)
library(treeClust)
library(POT)
library(erf)

################
# evaluation
################
source("methods.R")
source("evaluation.R")

data <- read.csv("../data/Amaurot.csv", stringsAsFactors = FALSE)
data[which(data$Season=="S1"),"Season"] <- 0
data[which(data$Season=="S2"),"Season"] <- 1
data$Season <- as.numeric(data$Season)

data2 <- data[complete.cases(data),] # remove NA's 
data3 <- data2[-which.max(data2$Y),] # remove max



# Daouia et al. (2011)
calculate_cond_SCV1(EQ.kernel.eva, tau0=1-1/10000, y=data2$Y, X=data2[,-1], seed=1)
calculate_cond_SCV2(EQ.kernel.eva, tau0=1-1/10000, y=data2$Y, X=data2[,-1], seed=1)

calculate_cond_SCV1(EQ.kernel.eva, tau0=1-1/10000, y=data3$Y, X=data3[,-1], seed=1)
calculate_cond_SCV2(EQ.kernel.eva, tau0=1-1/10000, y=data3$Y, X=data3[,-1], seed=1)

# TIR method of Wang and Tsai (2009)
calculate_cond_SCV1(TIR.EX.func.eva, tau0=1-1/10000, y=data2$Y, X=data2[,-1], seed=1) # cannot allocate vector of size
calculate_cond_SCV2(TIR.EX.func.eva, tau0=1-1/10000, y=data2$Y, X=data2[,-1], seed=1)

calculate_cond_SCV1(TIR.EX.func.eva, tau0=1-1/10000, y=data3$Y, X=data2[,-1], seed=1) # cannot allocate vector of size
calculate_cond_SCV2(TIR.EX.func.eva, tau0=1-1/10000, y=data3$Y, X=data2[,-1], seed=1)

# Chernouchokov and Du (2006)
calculate_cond_SCV1(GPD.EX.func1.eva, tau0=1-1/10000, y=data2$Y, X=as.matrix(data2[,-1]), seed=1) 
calculate_cond_SCV2(GPD.EX.func1.eva, tau0=1-1/10000, y=data2$Y, X=as.matrix(data2[,-1]), seed=1)

calculate_cond_SCV1(GPD.EX.func1.eva, tau0=1-1/10000, y=data3$Y, X=as.matrix(data3[,-1]), seed=1) 
calculate_cond_SCV2(GPD.EX.func1.eva, tau0=1-1/10000, y=data3$Y, X=as.matrix(data3[,-1]), seed=1) 

# quantile regression
calculate_cond_SCV1(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=1) 
calculate_cond_SCV2(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=1) 

calculate_cond_SCV1(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=1)
calculate_cond_SCV2(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=1) 

# QR
calculate_cond_SCV1(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=1)
calculate_cond_SCV2(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=1)

calculate_cond_SCV1(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=1)
calculate_cond_SCV2(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=1)

# evgam
calculate_cond_SCV1(evgam.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=1) 
calculate_cond_SCV2(evgam.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=1) 

calculate_cond_SCV1(evgam.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=1) 
calculate_cond_SCV2(evgam.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=1) 

# EQR / Wang et al. (2012)
calculate_cond_SCV1(wrap.TwoStage, tau0=1-1/10000, data2$Y, data2[,-1], seed=1)
calculate_cond_SCV2(wrap.TwoStage, tau0=1-1/10000, data2$Y, data2[,-1], seed=1)

calculate_cond_SCV1(wrap.TwoStage, tau0=1-1/10000, data3$Y, data3[,-1], seed=1)
calculate_cond_SCV2(wrap.TwoStage, tau0=1-1/10000, data3$Y, data3[,-1], seed=1)

# ERF
calculate_cond_SCV1(wrap.erf, tau0=1-1/10000, data2$Y, data2[,-1], seed=1)
calculate_cond_SCV2(wrap.erf, tau0=1-1/10000, data2$Y, data2[,-1], seed=1)

calculate_cond_SCV1(wrap.erf, tau0=1-1/10000, data3$Y, data3[,-1], seed=1)
calculate_cond_SCV2(wrap.erf, tau0=1-1/10000, data3$Y, data3[,-1], seed=1)

# GB
calculate_cond_SCV1(wrap.gb, tau0=1-1/10000, data2$Y, data2[,-1], seed=1)
calculate_cond_SCV2(wrap.gb, tau0=1-1/10000, data2$Y, data2[,-1], seed=1)

calculate_cond_SCV1(wrap.gb, tau0=1-1/10000, data3$Y, data3[,-1], seed=1)
calculate_cond_SCV2(wrap.gb, tau0=1-1/10000, data3$Y, data3[,-1], seed=1)
