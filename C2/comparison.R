library(Expectrem)
library(extremefit)
library(expint)
library(evmix)
library(extRemes)
library(texmex)
library(mev)
library(ggplot2)
library(dplyr)
library(tea)
library(ReIns)

################
# evaluation
################
source("methods.R")
source("evaluation.R")
source("evt0.R")

data <- read.csv("../data/Amaurot.csv")
data3 <- data[-which.max(data$Y),] # remove max


# Weissman-type estimator
calculate_loss1(predfunc=extQuant.eva, tau0=1-1/60000, y=data$Y, seed=1)
calculate_SCV1(extQuant.eva, tau0=1-1/60000, data$Y, seed=1) 
calculate_SCV2(extQuant.eva, tau0=1-1/60000, data$Y, seed=1) 

calculate_loss1(predfunc=extQuant.eva, tau0=1-1/60000, y=data3$Y, seed=1)
calculate_SCV1(extQuant.eva, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(extQuant.eva, tau0=1-1/60000, data3$Y, seed=1) 

# Weissman with optimal k / Caeiro, J. and Gomes, M.I. (2016)
calculate_loss1(predfunc=wrap.hill, tau0=1-1/60000, y=data$Y, seed=1)
calculate_SCV1(wrap.hill, tau0=1-1/60000, data$Y, seed=1) 
calculate_SCV2(wrap.hill, tau0=1-1/60000, data$Y, seed=1) 

calculate_loss1(predfunc=wrap.hill, tau0=1-1/60000, y=data3$Y, seed=1)
calculate_SCV1(wrap.hill, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(wrap.hill, tau0=1-1/60000, data3$Y, seed=1) 

# Composite extreme quantile estimation
calculate_loss1(predfunc=extQuantlp.eva, tau0=1-1/60000, y=data$Y, seed=1)
calculate_SCV1(extQuantlp.eva, tau0=1-1/60000, data$Y, seed=1)
calculate_SCV2(extQuantlp.eva, tau0=1-1/60000, data$Y, seed=1) 

calculate_loss1(predfunc=extQuantlp.eva, tau0=1-1/60000, y=data3$Y, seed=1)
calculate_SCV1(extQuantlp.eva, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(extQuantlp.eva, tau0=1-1/60000, data3$Y, seed=1)

# expectile
calculate_loss1(predfunc=extExpect.eva, tau0=1-1/60000, y=data$Y, seed=1)
calculate_SCV1(extExpect.eva, tau0=1-1/60000, data$Y, seed=1)
calculate_SCV2(extExpect.eva, tau0=1-1/60000, data$Y, seed=1)

calculate_loss1(predfunc=extExpect.eva, tau0=1-1/60000, y=data3$Y, seed=1)
calculate_SCV1(extExpect.eva, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(extExpect.eva, tau0=1-1/60000, data3$Y, seed=1) 

# extremefit
calculate_loss1(predfunc=extremefit.eva, tau0 = 1-1/60000, data$Y, seed=1)
calculate_SCV1(extremefit.eva, tau0=1-1/60000, data$Y, seed=1) 
calculate_SCV2(extremefit.eva, tau0=1-1/60000, data$Y, seed=1) 

calculate_loss1(predfunc=extremefit.eva, tau0 = 1-1/60000, data3$Y, seed=1)
calculate_SCV1(extremefit.eva, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(extremefit.eva, tau0=1-1/60000, data3$Y, seed=1) 

# evmix
calculate_loss1(predfunc=wrap.evmix, tau0 = 1-1/60000, data$Y, seed=1)
calculate_SCV1(wrap.evmix, tau0=1-1/60000, data$Y, seed=1) 
calculate_SCV2(wrap.evmix, tau0=1-1/60000, data$Y, seed=1) 

calculate_loss1(predfunc=wrap.evmix, tau0 = 1-1/60000, data3$Y, seed=1)
calculate_SCV1(wrap.evmix, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(wrap.evmix, tau0=1-1/60000, data3$Y, seed=1) 

# GPD (extRemes)
calculate_loss1(predfunc=wrap.fevd, tau0=1-1/60000, y=data$Y, seed=1)
calculate_SCV1(wrap.fevd, tau0=1-1/60000, data$Y, seed=1) 
calculate_SCV2(wrap.fevd, tau0=1-1/60000, data$Y, seed=1) 

calculate_loss1(predfunc=wrap.fevd, tau0=1-1/60000, y=data3$Y, seed=1)
calculate_SCV1(wrap.fevd, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(wrap.fevd, tau0=1-1/60000, data3$Y, seed=1) 


# GPD (texmex)
calculate_loss1(predfunc=wrap.evm, tau0=1-1/60000, y=data$Y, seed=1)
calculate_SCV1(wrap.evm, tau0=1-1/60000, data$Y, seed=1) 
calculate_SCV2(wrap.evm, tau0=1-1/60000, data$Y, seed=1) 

calculate_loss1(predfunc=wrap.evm, tau0=1-1/60000, y=data3$Y, seed=1)
calculate_SCV1(wrap.evm, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(wrap.evm, tau0=1-1/60000, data3$Y, seed=1) 

# refined Weissman
calculate_loss1(predfunc=estimate.rw, tau0=1-1/60000, y=data$Y, seed=1) ### loss for RW
calculate_SCV1(estimate.rw, tau0=1-1/60000, data$Y, seed=1) 
calculate_SCV2(estimate.rw, tau0=1-1/60000, data$Y, seed=1) 

calculate_loss1(predfunc=estimate.rw, tau0=1-1/60000, y=data3$Y, seed=1)
calculate_SCV1(estimate.rw, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(estimate.rw, tau0=1-1/60000, data3$Y, seed=1) 

# refined Weibull
calculate_loss1(predfunc=estimate.rweibull, tau0=1-1/60000, y=data$Y, seed=1)  ### loss for weibull 
calculate_SCV1(estimate.rweibull, tau0=1-1/60000, data$Y, seed=1) 
calculate_SCV2(estimate.rweibull, tau0=1-1/60000, data$Y, seed=1) 

calculate_loss1(predfunc=estimate.rweibull, tau0=1-1/60000, y=data3$Y, seed=1)
calculate_SCV1(estimate.rweibull, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(estimate.rweibull, tau0=1-1/60000, data3$Y, seed=1) 



#############
# test
#############

# Weissman-type estimator
extQuant(X=data$Y, tau=1-1/60000) 
extQuant(X=data3$Y, tau=1-1/60000) 

# Weissman with optimal k / Caeiro, J. and Gomes, M.I. (2016)
wrap.hill(tau=1-1/60000, data=data$Y) 
wrap.hill(tau=1-1/60000, data=data3$Y)

# Composite extreme quantile estimation (selected, )
extQuantlp(X=data$Y,tau=1-1/60000,k=50,p=1.5,estim="lpindex",br=TRUE) 
extQuantlp(X=data3$Y,tau=1-1/60000,k=50,p=1.5,estim="lpindex",br=TRUE) 

# expectile
extExpect.eva(tau=1-1/60000, train=data$Y)
extExpect.eva(tau=1-1/60000, train=data3$Y) 

# extremefit
extremefit.eva(tau=1-1/60000, train=data$Y) 
extremefit.eva(tau=1-1/60000, train=data3$Y) 

# evmix
wrap.evmix(tau=1-1/60000, data=data$Y)
wrap.evmix(tau=1-1/60000, data=data3$Y)

# GPD (extRemes)
wrap.fevd(tau=1-1/60000, data=data$Y)
wrap.fevd(tau=1-1/60000, data=data3$Y)

# GPD (texmex)
wrap.evm(tau=1-1/60000, data=data$Y)
wrap.evm(tau=1-1/60000, data=data3$Y)

# refined Weissman
estimate.rw(x=data$Y, tau=1-1/60000) ### estimate for RW
estimate.rw(x=data3$Y, tau=1-1/60000)

# refined Weibull
estimate.rweibull(x=data$Y, tau=1-1/60000) ### estimate for weibull
estimate.rweibull(x=data3$Y, tau=1-1/60000)
