library(caret)
library(quantreg)
library(mnormt)
library(ks)
library(evd)
library(rdist)
library(evmix)
library(evgam)

load(".Rdata")
source("evaluation.R")


#################################################
#  Kernel estimator of extreme cond. quantile   #
#		for heavy tailed distributions 	#
#			Daouia et al. (2011)	#
#          main function: EQ.kernel             #
#################################################
#beginning of the subroutine functions and function EQ.kernel
bi.kernel = function(x)
{
  15/16*(1-x^2)^2*(abs(x)<=1)
}

mypkde = function(Y, X, newy, newx, h)
{
  #kernel estimator of cumulative distribution function F(ynew|X=newx)=P(Y<=ynew |X=newx)
  #kernel function: product normal kernel
  
  #Inputs
  #Y: observed response vector, n*1
  #X: observed covariate matrix, n*p (no intercept)
  #newy: a n*1 vector
  #newx: a p-dimensional vector
  #h: bandwidth
  
  X = as.matrix(X)
  Y = as.vector(Y)
  newx = as.vector(newx)
  n = nrow(X)
  p = ncol(X)
  if(p != length(newx)) warning("the dimensions of X and newx do not match")
  Dx = t(t(X) - newx)
  
  K = apply(Dx/h, 2, bi.kernel)
  K = apply(K, 1, prod)/h^p
  K = K/sum(K)
  D = outer(Y, newy, "<=")
  Fhat =apply(D*K,2,sum)
  out = list(Fhat=Fhat, eval.points=newy)
  return(out)
}

mypkde.cv = function(Y, X, nh = 50, max.h=0.5)
{
  #select h by cross validation (very slow)
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  
  hs = seq(1/(5*log(n)), max.h, length=nh)
  obj = NULL
  for(h in hs)
  {
    # #leave the i-th observation out
    tmp.obj=0
    for(i in 1:n)
    {
      K = apply(t(t(X[-i,]) - X[i,])/h, 2, bi.kernel)
      K = apply(K, 1, prod)/h^p
      K = K/sum(K)
      D = outer(Y[-i], Y, "<=")
      Fhat.loo= apply(D*K,2,sum)
      tmp.obj=tmp.obj + sum((1*(Y[i]<=Y) - Fhat.loo)^2)
    }	
    obj = c(obj, tmp.obj)
  }
  idx = which.min(obj)
  h = hs[idx]
  return(h)
}

myqkde <- function (tau, Fhat)
{
  #tau-th conditional quantile based on kernel CDF estimation Fhat: output from mypkde function
  #tau can be a vector
  #tau = 0.75
  if (any(tau > 1) | any(tau < 0))
    stop("p must be <= 1 and >= 0")
  CDF = Fhat$Fhat
  ind <- findInterval(x = tau, vec = CDF)
  
  quant <- rep(0, length(ind))
  for (j in 1:length(ind)) {
    i <- ind[j]
    if (i == 0)
      quant[j] <- Fhat$eval.points[1]
    else if (i >= length(Fhat$eval.points))
      quant[j] <- Fhat$eval.points[length(Fhat$eval.points)]
    else {
      quant1 <- Fhat$eval.points[i]
      quant2 <- Fhat$eval.points[i + 1]
      prob1 <- CDF[i]
      prob2 <- CDF[i + 1]
      alpha <- (tau[j] - prob2)/(prob1 - prob2)
      quant[j] <- quant1 * alpha + quant2 * (1 - alpha)
    }
  }
  names(quant) = paste(round(tau,5))
  return(quant)
}

EQ.kernel = function(Y, X, newx, tau.e, k, h=NULL, cv=TRUE, nh, max.h)
{
  #Inputs
  #Y: observed response vector, n*1
  #X: observed covariate matrix, n*p (no intercept)
  #newy: a n*1 vector
  #newx: a m*p matrix
  #h: bandwidth
  
  #estimate the cond quantiles for not so extreme quantiles (level=alpha_n, the same as 1-k/n in our notation)
  n = length(Y)
  alpha.n = 1-k/n
  newy = binning(Y)$eval.points
  # select h
  if(is.null(h)==TRUE | cv==TRUE)
  {
    h= mypkde.cv(Y, X, nh, max.h)	
  }
  
  nx = nrow(newx)
  tauj = 1/(1:9)
  Q.e = NULL
  Gamma.x = rep(0, nx)
  for(i in 1:nx)
  {
    Fhat = mypkde(Y, X, newy, newx=as.vector(newx[i,]), h)
    Q1 = myqkde(1-tauj*(1-alpha.n), Fhat)
    gamma.x = sum(log(Q1)-log(Q1[1]))/sum(log(1/tauj))
    #obtained extrapolated extreme quantile estimator Q.e
    Q.e = rbind(Q.e, Q1[1]*((1-alpha.n)/(1-tau.e))^gamma.x)
    Gamma.x[i] = gamma.x
  }
  out = list(Q.e=Q.e, gamma.x = Gamma.x, h=h)
  return(out)
}

EQ.kernel.eva <- function(tau, train_y, train_x, valid_x, nh = 50, max.h=0.5){
  return(EQ.kernel(Y=train_y, X=train_x, newx=valid_x, tau.e=tau, k=1, h=NULL, cv=TRUE, nh=nh, max.h=max.h))
}

###################################################
#        Pareto distrbution quantile              #
###################################################

qpareto = function(tau, gamma)
{
  (1-tau)^(-gamma)
}
rpareto = function(n, gamma)
{
  u = runif(n, 0, 1)
  (1-u)^(-gamma)
}

########################################
#  TIR method of Wang and Tsai (2009)  #
########################################
TIR.coef.func=function(y, X, wn, theta0)
{
  #estimate theta for a given threshold wn
  #first element of X is 1 corresponding to the intercept
  #theta0: initial value of the parameter
  like.func = function(theta, y, X, wn)
  {
    xtheta = X%*%theta
    idx = which(y>wn)
    tmp = exp(xtheta[idx])*(log(pmax(y[idx],10e-5))-log(wn))-xtheta[idx]
    sum(tmp)
  }
  theta.hat=optim(theta0, like.func, y=y, X=X, wn=wn)$par
  return(theta.hat)
}

###select the optimal wn, return the estimated coefficient theta based on the optimal wn
TIP.func = function(y, X, theta0, nw=50, p0=1)
{
  p=ncol(X)
  n = length(y)
  wns = quantile(y, 1-seq(4*p/n, p0,length=nw))
  wns=pmax(wns,0.0001)
  gof = rep(0, length(wns))
  thetas = NULL
  for(i in 1:length(wns))
  {
    wn=wns[i]
    theta=TIR.coef.func(y, X, wn, theta0)
    xtheta = X%*%theta
    Ui = exp(-exp(xtheta)*(log(pmax(y,10e-10))-log(wn)))
    n0 = sum(y>wn)
    idx = which(y>wn)
    Ui = Ui[idx]
    Fn = rank(Ui)/n0
    # plot(Fn[idx]~Ui[idx],main=paste("wn=",wn)); abline(0,1,col=2)
    gof[i] = mean((Ui-Fn)^2)
    thetas = rbind(thetas, theta)
  }
  idx = order(gof)[1]
  wn = wns[idx]
  #plot(gof~wns); abline(v=wn, col=2)
  theta = thetas[idx,]
  return(theta)
}


#### estimate the extreme conditional quantiles based on TIR estimation
TIR.EX.func = function(xstar, new.tau, y, X, theta0, nw=50, p0=1)
{
  # estimate the new.tau-th conditional quantile of Y given xstar
  # xstar: first element is 1
  # k: is chosen to be as.integer(n^(1/3)*4.5)
  # X: first element is 1
  # theta0: initial parameter value, same dimension as X
  xstar = as.matrix(xstar)
  n = length(y)
  nx = nrow(xstar)
  theta = TIP.func(y, X, theta0, nw, p0)
  hat.gamma = 1/exp(xstar%*%theta)
  #extrapolation
  k = as.integer(n^(1/3)*4.5)
  tau0 = 1-k/n
  q.tau = xstar%*% rq(y~X-1, tau=tau0)$coef
  Qhat.TIR = matrix(0, nrow=nrow(xstar), ncol=length(new.tau))
  for(j in 1:nx)
  {
    Qhat.TIR[j,]= ((1-tau0)/(1-new.tau))^(hat.gamma[j]) * q.tau[j]
  }
  return(list(Qhat=Qhat.TIR, theta.hat=theta))
}

TIR.EX.func.eva <- function(tau, train_y, train_x, valid_x, 
                            theta0=matrix(nrow=dim(train_x)[1], ncol=dim(train_x)[2]), nw=50, p0=1){
  return(TIR.EX.func(xstar = valid_x, new.tau=tau, y=train_y, X=train_x,
                     theta0=theta0, nw=nw, p0=p0))
}

####### Method in Chernouchokov and Du (2006) based on Hill's estimator
CDu.EX.func = function(xstar, new.tau, y, X)
{
  # xstar: first element is 1
  # k: is chosen to be as.integer(n^(1/3)*4.5)
  # X: first element is 1
  
  n = length(y)
  k = as.integer(n^(1/3)*4.5)
  tau0 = 1-k/n
  rqfit = rq(y~X-1,tau=tau0)
  bhat = rqfit$coef
  yhat = rqfit$fitted
  idx = which(yhat>0 & y>yhat)
  EVI = sum(log(y[idx])-log(yhat[idx]))/(n*(1-tau0))
  Q.tau0 = xstar%*%bhat
  Qhat = as.vector(Q.tau0) %*% t(as.vector(((1-tau0)/(1-new.tau))^EVI))
  return(list(Qhat=Qhat, EVI=EVI))
}

GPD.loglike = function(par, z, X)
{
  p = ncol(X)
  beta = par[1:p]
  theta = par[-(1:p)]
  sigma = exp(X%*%beta)
  xi = exp(X%*%theta)
  f = 1/sigma*(1+z*xi/sigma)^(-1/xi-1)
  logf = log(f)
  return(-sum(logf))
}


GPD.loglike2 = function(par, z, X)
{
  p = ncol(X)
  beta = par[1:p]
  theta = par[-(1:p)]
  sigma = exp(X%*%beta)
  xi = exp(X%*%theta)
  obj = 0 #negative log likelihood
  m = length(z)
  for(j in 1:m)
  {
    obj = obj - dgpd(z[j], loc=0, scale=sigma[j], shape=xi[j], log=TRUE)
  }
  return(obj)
}



GPD.EX.func1 = function(xstar, new.tau, y, X, theta0)
{
  #assume log(scale)=X*beta; log(EVI)=X*theta
  n = length(y)
  k = as.integer(n^(1/3)*4.5)
  tau0 = 1-k/n
  rqfit = rq(y~X-1,tau=tau0)
  u = rqfit$fitted  #threshold
  #exceedances
  idx = which(y>u)
  z = y[idx]-u[idx]
  X2 = X[idx,]
  p = ncol(X)
  beta0 = rep(0,p)
  ###obtain MLE	
  tt = optim(c(beta0, theta0), fn=GPD.loglike, z=z, X=X2)
  beta = tt$par[1:p]
  theta = tt$par[-(1:p)]
  ### estimate quantile
  sigma = exp(xstar%*%beta)
  xi = exp(xstar%*%theta)
  p.u = 1-tau0
  ntau = length(new.tau)
  Qhat = matrix(0, nrow=nrow(xstar), ncol=ntau)
  for(j in 1:ntau)
    # Qhat[,j] = u + sigma*((xi==0)*(log(1-new.tau[j])-log(p.u)) + (xi!=0)/xi*(((1-new.tau[j])/p.u)^(-xi)-1))
    Qhat[,j] = xstar%*%rqfit$coefficients + sigma*((xi==0)*(log(1-new.tau[j])-log(p.u)) + (xi!=0)/xi*(((1-new.tau[j])/p.u)^(-xi)-1))
  return(list(Qhat=Qhat, par=tt$par))
}

GPD.EX.func1.eva <-function(tau, train_y, train_x, valid_x){
  return(GPD.EX.func1(xstar=valid_x, new.tau=tau, y=train_y, X=valid_x, theta0 = rep(0, ncol(valid_x))))
}


########################
#  Quantile regression #
########################
qrfunc.eva <- function(tau, train_y, train_x, valid_x){
  data <- cbind(train_y, train_x)
  colnames(data) <- c("Y", colnames(train_x))
  rqfit <- rq(Y ~ ., tau = tau, data = data)
  return(predict(rqfit, newdata = valid_x))
}

qrfunc.eva.ci <- function(tau, train_y, train_x, valid_x){
  data <- cbind(train_y, train_x)
  colnames(data) <- c("Y", colnames(train_x))
  rqfit <- rq(Y ~ ., tau = tau, data = data)
  return(predict(rqfit, newdata = valid_x, interval = "confidence", level=0.5, se="boot"))
}


############
#  EVA GAM #
############
evgam.eva <-function(tau,train_y, train_x, valid_x){
  train<-as.data.frame(cbind(train_y,train_x))
  form<-list(train_y~ V1+V2+V3+V4+Season+WindDirection+WindSpeed+Atmosphere , ~ 1 )
  newdata<-as.data.frame(valid_x)
  ald<-evgam(form, data=train, family="ald", ald.args=list(tau=0.7))
  train$excess<-train_y-predict(ald)$location
  is.na(train$excess[train$excess<=0])<-TRUE
  form<-list(excess~ V1+V2+V3+V4+Season+WindDirection+WindSpeed+Atmosphere , ~ 1 )
  mod<-evgam(form, data=train, family="gpd")
  return(unlist(predict(mod,newdata=newdata,type="quantile",prob=.9999)[,1]+predict(ald, newdata = newdata)$location)) 
}


###############
data <- read.csv("../data/Amaurot.csv", stringsAsFactors = FALSE)
data[which(data$Season=="S1"),"Season"] <- 0
data[which(data$Season=="S2"),"Season"] <- 1
data$Season <- as.numeric(data$Season)

data2 <- data[complete.cases(data),] # remove NA's 
data3 <- data2[-which.max(data2$Y),] # remove max



# Daouia et al. (2011)
calculate_cond_SCV1(EQ.kernel.eva, tau0=1-1/10000, y=data2$Y, X=data2[,-1], seed=1)


# TIR method of Wang and Tsai (2009)
calculate_cond_SCV1(TIR.EX.func.eva, tau0=1-1/10000, y=data2$Y, X=data2[,-1], seed=1)
# cannot allocate vector of size

# Chernouchokov and Du (2006)
calculate_cond_SCV1(GPD.EX.func1.eva, tau0=1-1/10000, y=data2$Y, X=as.matrix(data2[,-1]), seed=1) 
calculate_cond_SCV1(GPD.EX.func1.eva, tau0=1-1/10000, y=data2$Y, X=as.matrix(data2[,-1]), seed=100) 
calculate_cond_SCV1(GPD.EX.func1.eva, tau0=1-1/10000, y=data2$Y, X=as.matrix(data2[,-1]), seed=200) 
calculate_cond_SCV1(GPD.EX.func1.eva, tau0=1-1/10000, y=data2$Y, X=as.matrix(data2[,-1]), seed=300) 

calculate_cond_SCV2(GPD.EX.func1.eva, tau0=1-1/10000, y=data2$Y, X=as.matrix(data2[,-1]), seed=1)
calculate_cond_SCV2(GPD.EX.func1.eva, tau0=1-1/10000, y=data2$Y, X=as.matrix(data2[,-1]), seed=100) 
calculate_cond_SCV2(GPD.EX.func1.eva, tau0=1-1/10000, y=data2$Y, X=as.matrix(data2[,-1]), seed=200) 
calculate_cond_SCV2(GPD.EX.func1.eva, tau0=1-1/10000, y=data2$Y, X=as.matrix(data2[,-1]), seed=300)

calculate_cond_SCV1(GPD.EX.func1.eva, tau0=1-1/10000, y=data3$Y, X=as.matrix(data3[,-1]), seed=1) 
calculate_cond_SCV1(GPD.EX.func1.eva, tau0=1-1/10000, y=data3$Y, X=as.matrix(data3[,-1]), seed=100) 
calculate_cond_SCV1(GPD.EX.func1.eva, tau0=1-1/10000, y=data3$Y, X=as.matrix(data3[,-1]), seed=200) 
calculate_cond_SCV1(GPD.EX.func1.eva, tau0=1-1/10000, y=data3$Y, X=as.matrix(data3[,-1]), seed=300) 

calculate_cond_SCV2(GPD.EX.func1.eva, tau0=1-1/10000, y=data3$Y, X=as.matrix(data3[,-1]), seed=1) 
calculate_cond_SCV2(GPD.EX.func1.eva, tau0=1-1/10000, y=data3$Y, X=as.matrix(data3[,-1]), seed=100) 
calculate_cond_SCV2(GPD.EX.func1.eva, tau0=1-1/10000, y=data3$Y, X=as.matrix(data3[,-1]), seed=200) 
calculate_cond_SCV2(GPD.EX.func1.eva, tau0=1-1/10000, y=data3$Y, X=as.matrix(data3[,-1]), seed=300) 

# quantile regression
calculate_cond_SCV1(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=1) 
calculate_cond_SCV1(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=100) 
calculate_cond_SCV1(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=200)
calculate_cond_SCV1(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=300) 

calculate_cond_SCV2(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=1) 
calculate_cond_SCV2(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=100)
calculate_cond_SCV2(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=200) 
calculate_cond_SCV2(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=300) 

calculate_cond_SCV1(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=1)
calculate_cond_SCV1(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=100) 
calculate_cond_SCV1(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=200)
calculate_cond_SCV1(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=300) 

calculate_cond_SCV2(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=1) 
calculate_cond_SCV2(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=100) 
calculate_cond_SCV2(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=200)
calculate_cond_SCV2(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=300) 


# QR
calculate_cond_SCV1_kernel(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=1)
calculate_cond_SCV1_kernel(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=100)
calculate_cond_SCV2_kernel(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=1) 
calculate_cond_SCV2_kernel(qrfunc.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=100)

calculate_cond_SCV1_kernel(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=1)
calculate_cond_SCV1_kernel(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=100) 
calculate_cond_SCV2_kernel(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=1) 
calculate_cond_SCV2_kernel(qrfunc.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=100) 

# evgam
calculate_cond_SCV1_kernel(evgam.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=1) 
calculate_cond_SCV1_kernel(evgam.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=100) 
calculate_cond_SCV2_kernel(evgam.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=1) 
calculate_cond_SCV2_kernel(evgam.eva, tau0=1-1/10000, data2$Y, data2[,-1], seed=100) 

calculate_cond_SCV1_kernel(evgam.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=1) 
calculate_cond_SCV1_kernel(evgam.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=100) 
calculate_cond_SCV2_kernel(evgam.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=1) 
calculate_cond_SCV2_kernel(evgam.eva, tau0=1-1/10000, data3$Y, data3[,-1], seed=100) 
