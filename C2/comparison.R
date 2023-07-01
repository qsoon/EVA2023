library(Expectrem)
library(extremefit)
################
CIextExpect.modify <- function (X, k = trunc(length(X)/10), tau, estim="Hill",method = "direct", 
                                ci.level = 0.95) 
{
  n = length(X)
  mopest = mop(X, 1:(n - 1), 0, method = "RBMOP")
  if (k[1] == "kopt" && estim == "Hill") {
    k = trunc(((1 - mopest$rho)^2/(-2 * mopest$rho * mopest$beta^2))^(1/(1 - 
                                                                           2 * mopest$rho)) * n^(-2 * mopest$rho/(1 - 2 * mopest$rho)))
  }
  if (k[1] == "kopt" && estim == "tindexp") {
    k = trunc(((1 - mopest$rho)^2/(-2 * mopest$rho * mopest$beta^2))^(1/(1 - 
                                                                           2 * mopest$rho)) * n^(-2 * mopest$rho/(1 - 2 * mopest$rho)))
    k = min(trunc(((1/mopest$EVI[k] - 1)^(2 * mopest$rho - 
                                            1) * (1 - mopest$EVI[k] - mopest$rho)^2/(-2 * mopest$rho * 
                                                                                       mopest$beta^2 * abs(1 - 2 * mopest$EVI[k])))^(1/(1 - 
                                                                                                                                          2 * mopest$rho)) * n^(-2 * mopest$rho/(1 - 2 * mopest$rho))), 
            trunc(n/2) - 1)
  }
  if (length(k) > 1) {
    stop("k must be of length 1.")
  }
  if (length(tau) > 1) {
    stop("tau must be of length 1.")
  }
  if (k > n - 1 || k < 1) {
    stop("k must be between 1 and n-1.")
  }
  if (method != "direct" && method != "indirect") {
    stop("method may be either direct or indirect.")
  }
  if (method == "direct") {
    qtp = expect(X, 1 - k/n)
    gammahat = tindexp(X, k, br = T)
    if (gammahat > 0.5) {
      print("WARNING : Tail index above 1/2 ! Use the indirect method rather ?")
    }
    r = (1 - mean(X)/qtp) * (n/(n - 2 * k)) * (1 + mopest$beta * 
                                                 gammahat * Fbar(X, qtp)^(-mopest$rho)/(gammahat * 
                                                                                          (1 - mopest$rho - gammahat)))^(-1)
    rbet = (1 - mean(X)/(qtp * (k/(n * (1 - tau)))^(gammahat))) * 
      (1/(2 * tau - 1)) * (1 + mopest$beta * gammahat * 
                             (1/gammahat - 1)^(-mopest$rho) * (1 - tau)^(-mopest$rho)/(gammahat * 
                                                                                         (1 - mopest$rho - gammahat)))^(-1)
    estimpoint = qtp * (k/(n * (1 - tau)))^gammahat * (1 + 
                                                         ((k/(n * (1 - tau)))^mopest$rho - 1)/mopest$rho * 
                                                         mopest$beta * gammahat * (n/k)^mopest$rho) * 
      (r/rbet)^gammahat * (1 + ((1/gammahat - 1)^(-mopest$rho) * 
                                  rbet^(-mopest$rho) - 1)/mopest$rho * mopest$beta * 
                             gammahat * (1 - tau)^(-mopest$rho))/(1 + ((1/gammahat - 
                                                                          1)^(-mopest$rho) * r^(-mopest$rho) - 1)/mopest$rho * 
                                                                    mopest$beta * gammahat * (k/n)^(-mopest$rho))
    psihat = mean((X - qtp) * (X > qtp))
    Fbarhat = mean(X > qtp)
    quantint = quantile(X, 1 - Fbarhat)
    psihatquant = mean((X - quantint) * (X > quantint))
    psihatbis = mean((X - max(qtp, quantint)) * (X > max(qtp, 
                                                         quantint)))
    Fbarhatbis = mean(X > max(qtp, quantint))
    psi2hat = min(max(2 * Fbarhat * qtp^2 * gammahat^2 * 
                        (1/abs((1 - gammahat) * (1 - 2 * gammahat)) + Fbarhat^(-mopest$rho) * 
                           mopest$beta/mopest$rho * (1/abs((1 - gammahat - 
                                                              mopest$rho) * (1 - 2 * gammahat - mopest$rho)) - 
                                                       1/abs((1 - gammahat) * (1 - 2 * gammahat)))), 
                      psihat^2), sqrt(mean((X - qtp)^4 * (X > qtp))))
    mhat = mean(X)
    alphan = 1 - Fbarhat
    M11phi = k/n * (psi2hat/psihat^2 - 1)
    M12phi = k/n * ((psihatbis + (apply(cbind(qtp, quantint), 
                                        1, max) - qtp) * Fbarhatbis)/(psihat * (1 - alphan)) - 
                      1)
    M22phi = k/n * alphan/(1 - alphan)
    M11xi = (psihat/qtp)^2 * ((qtp - mhat)^2 * M11phi)/(psihat + 
                                                          Fbarhat * (qtp - mhat))^2
    M12xi = gammahat * psihat/qtp * ((qtp - mhat) * M12phi)/(psihat + 
                                                               Fbarhat * (qtp - mhat))
    M22xi = gammahat^2 * M22phi
    epsilon = (1/gammahat - 1)^(-gammahat) * (1 - mhat/qtp)^(-gammahat) * 
      (1 - 2 * k/n)^gammahat * (Fbarhat/(k/n))^gammahat
    M11 = 1/gammahat^2 * (epsilon^2 * (Fbarhat/(k/n))^2 * 
                            M11xi - 2 * epsilon * (Fbarhat/(k/n)) * M12xi + 
                            M22xi)
    M12 = 1/gammahat * (M12xi - epsilon * (Fbarhat/(k/n)) * 
                          M11xi)
    M22 = M11xi
    S11 = M11 * 1/(1 + Fbarhat/(k/n))^4 * (1 + 8 * M11 * 
                                             (1 + Fbarhat/(k/n))^(-2)/(k))
    S12 = -1/(1 + Fbarhat/(k/n))^2 * M12 * (1 + 3 * M11 * 
                                              (1 + Fbarhat/(k/n))^(-2)/k)
    S22 = M22
    nablah1 = qtp * (1 - 2 * k/n) * (qtp - mhat)/((mhat - 
                                                     2 * qtp * k/n) * 1/(1 + Fbarhat/(k/n)) + 2 * qtp * 
                                                    k/n - qtp)^2
    nablah2 = qtp * (1 - 2 * k/n) * (1 - 1/(1 + Fbarhat/(k/n))) * 
      1/(1 + Fbarhat/(k/n)) * mhat/(2 * qtp * k/n * 1/(1 + 
                                                         Fbarhat/(k/n)) + qtp - 2 * qtp * k/n - 1/(1 + Fbarhat/(k/n)) * 
                                      mhat)^2
    S11prime = nablah1^2 * S11 + nablah2^2 * S22 + 2 * nablah1 * 
      nablah2 * S12
    S12prime = nablah1 * S12 + nablah2 * S22
    S22prime = S22
    h = gammahat
    nabla1 = log((2 * tau - 1)/(1 - 2 * k/n)) - log((1 - 
                                                       tau)/(k/n)) + log(1 - mhat/qtp) - log(1 - mhat/(qtp * 
                                                                                                         ((1 - tau)/(k/n))^(-h))) - h * log((1 - tau)/(k/n)) * 
      mhat/(mhat - qtp * ((1 - tau)/(k/n))^(-h))
    nabla2 = 1 - h * mhat/(qtp * ((1 - tau)/(k/n))^(-h) - 
                             mhat) + h * mhat/(qtp - mhat)
    varcorr = nabla1^2 * S11prime/log(k/(n * (1 - tau)))^2 + 
      nabla2^2 * S22prime/log(k/(n * (1 - tau)))^2 + 2 * 
      nabla1 * nabla2 * S12prime/log(k/(n * (1 - tau)))^2
    estimup = estimpoint * exp(-(sqrt(varcorr)/sqrt(k) * 
                                   log((k/n)/(1 - tau)) * qnorm((1 - ci.level)/2)))
    estimdown = estimpoint * exp(-(sqrt(varcorr)/sqrt(k) * 
                                     log((k/n)/(1 - tau)) * qnorm((1 + ci.level)/2)))
    return(list(Lower_bound = estimdown, Point_estimate = estimpoint, 
                Upper_bound = estimup, gammahat=gammahat))
  }
  if (method == "indirect") {
    gammahat = mopest$EVI[k]
    if (gammahat > 1) {
      stop("Tail index greater than 1 ! Expectile does probably not exist !")
    }
    qtp = quantile(X, 1 - k/n) * (1/gammahat - 1)^(-gammahat)
    r = (1 - mean(X)/qtp) * (n/(n - 2 * k)) * (1 + mopest$beta * 
                                                 gammahat * Fbar(X, qtp)^(-mopest$rho)/(gammahat * 
                                                                                          (1 - mopest$rho - gammahat)))^(-1)
    rbet = (1 - mean(X)/(qtp * (k/(n * (1 - tau)))^(gammahat))) * 
      (1/(2 * tau - 1)) * (1 + mopest$beta * gammahat * 
                             (1/gammahat - 1)^(-mopest$rho) * (1 - tau)^(-mopest$rho)/(gammahat * 
                                                                                         (1 - mopest$rho - gammahat)))^(-1)
    estimpoint = (1/gammahat - 1)^(-gammahat) * quantile(X, 
                                                         1 - k/n) * (k/(n * (1 - tau)))^gammahat * (1 + ((k/(n * 
                                                                                                               (1 - tau)))^mopest$rho - 1)/mopest$rho * mopest$beta * 
                                                                                                      gammahat * (n/k)^mopest$rho) * (1 + ((1/gammahat - 
                                                                                                                                              1)^(-mopest$rho) * rbet^(-mopest$rho) - 1)/mopest$rho * 
                                                                                                                                        mopest$beta * gammahat * (1 - tau)^(-mopest$rho))/(rbet^gammahat)
    quanthat = quantile(X, 1 - k/n)
    expectint = expect(X, 1 - k/n)
    mhat = mean(X)
    V11 = 1
    V12 = (1/(1 - gammahat) - log(1/gammahat - 1))
    V12 = V12 + (3 * gammahat - 1)/(2 * (1 - gammahat)^3 * 
                                      k)
    V12 = V12 + 3 * (10 * gammahat^3 - 10 * gammahat^2 + 
                       5 * gammahat - 1)/(4 * (1 - gammahat)^5 * k^2)
    V22 = (1 + (1/(1 - gammahat) - log(1/gammahat - 1))^2)
    V22 = V22 + (1/(1 - gammahat) - log(1/gammahat - 1)) * 
      (3 * gammahat - 1)/((1 - gammahat)^3 * k) + 0.5/((1 - 
                                                          gammahat)^4 * k)
    V22 = V22 + 36 * (10 * gammahat^3 - 10 * gammahat^2 + 
                        5 * gammahat - 1) * (1 - (1 - gammahat) * log(1/gammahat - 
                                                                        1))/(24 * (1 - gammahat)^6 * k^2) + (6 * gammahat^2 - 
                                                                                                               4 * gammahat + 1)/((1 - gammahat)^6 * k^2) + 20 * 
      (3 * gammahat - 1)^2/(48 * (1 - gammahat)^6 * k^2)
    nabla1 = 1 - gammahat * mhat/((n * (1 - tau)/k)^(-gammahat) * 
                                    (1/gammahat - 1)^(-gammahat) * quanthat - mhat) + 
      (log(2 * tau - 1) - log(1 - mhat/((n * (1 - tau)/k)^(-gammahat) * 
                                          (1/gammahat - 1)^(-gammahat) * quanthat)))/log(k/(n * 
                                                                                              (1 - tau)))
    nabla2 = (1 - mhat * gammahat/(quanthat * (1/gammahat - 
                                                 1)^(-gammahat) * (k/(n * (1 - tau)))^gammahat - 
                                     mhat))/log(k/(n * (1 - tau)))
    varcorr = nabla1^2 * V11 * gammahat^2 + nabla2^2 * V22 * 
      gammahat^2 + 2 * nabla1 * nabla2 * V12 * gammahat^2
    estimup = estimpoint * exp(-(sqrt(varcorr)/sqrt(k) * 
                                   log((k/n)/(1 - tau)) * qnorm((1 - ci.level)/2)))
    estimdown = estimpoint * exp(-(sqrt(varcorr)/sqrt(k) * 
                                     log((k/n)/(1 - tau)) * qnorm((1 + ci.level)/2)))
    return(list(Lower_bound = estimdown, Point_estimate = estimpoint, 
                Upper_bound = estimup, gammahat=gammahat))
  }
}


CIextExpect.eva <- function(tau, train){
  res <- CIextExpect.modify(X=train, tau=tau, k="kopt")
  expect <- res$Point_estimate
  tail.ind <- res$gammahat
  return(expect * (1/tail.ind -1)^tail.ind)
}

extExpect.modify <- function (X, k = "kopt", tau, estim = "Hill", method = "direct", 
                              br = FALSE) 
{
  n = length(X)
  mopest = mop(X, 1:(n - 1), 0, method = "RBMOP")
  if (k[1] == "kopt" && estim == "Hill") {
    k = trunc(((1 - mopest$rho)^2/(-2 * mopest$rho * mopest$beta^2))^(1/(1 - 
                                                                           2 * mopest$rho)) * n^(-2 * mopest$rho/(1 - 2 * mopest$rho)))
  }
  if (k[1] == "kopt" && estim == "tindexp") {
    k = trunc(((1 - mopest$rho)^2/(-2 * mopest$rho * mopest$beta^2))^(1/(1 - 
                                                                           2 * mopest$rho)) * n^(-2 * mopest$rho/(1 - 2 * mopest$rho)))
    k = min(trunc(((1/mopest$EVI[k] - 1)^(2 * mopest$rho - 
                                            1) * (1 - mopest$EVI[k] - mopest$rho)^2/(-2 * mopest$rho * 
                                                                                       mopest$beta^2 * abs(1 - 2 * mopest$EVI[k])))^(1/(1 - 
                                                                                                                                          2 * mopest$rho)) * n^(-2 * mopest$rho/(1 - 2 * mopest$rho))), 
            trunc(n/2) - 1)
  }
  if (any(k > n - 1) || any(k < 1)) {
    stop("k must be between 1 and n-1.")
  }
  if (estim != "Hill" && estim != "tindexp") {
    stop("estim may be either Hill or tindexp.")
  }
  if (method != "direct" && method != "indirect") {
    stop("method may be either direct or indirect.")
  }
  if (!is.logical(br)) {
    stop("br must be boolean.")
  }
  if (estim == "Hill" && br == TRUE) {
    gammahat = mopest$EVI[k]
  }
  if (estim == "Hill" && br == FALSE) {
    gammahat = mop(X, 1:(n - 1), 0, method = "MOP")$EVI[k]
  }
  if (estim == "tindexp") {
    gammahat = tindexp(X, k, br)
  }
  qtp = expect(X, 1 - k/n)
  r = (1 - mean(X)/qtp) * (n/(n - 2 * k)) * (1 + mopest$beta * 
                                               gammahat * Fbar(X, qtp)^(-mopest$rho)/(gammahat * (1 - 
                                                                                                    mopest$rho - gammahat)))^(-1)
  rbet = (1 - mean(X)/(qtp * (k/(n * (1 - tau)))^(gammahat))) * 
    (1/(2 * tau - 1)) * (1 + mopest$beta * gammahat * (1/gammahat - 
                                                         1)^(-mopest$rho) * (1 - tau)^(-mopest$rho)/(gammahat * 
                                                                                                       (1 - mopest$rho - gammahat)))^(-1)
  if (br == FALSE && method == "direct") {
    return(c(gammahat, qtp * (k/(n * (1 - tau)))^gammahat))
  }
  if (br == FALSE && method == "indirect") {
    return(c(gammahat, (1/gammahat - 1)^(-gammahat) * quantile(X, 1 - 
                                                                 k/n) * (k/(n * (1 - tau)))^gammahat))
  }
  if (br == TRUE && method == "direct") {
    return(c(gammahat, qtp * (k/(n * (1 - tau)))^gammahat * (1 + ((k/(n * 
                                                                        (1 - tau)))^mopest$rho - 1)/mopest$rho * mopest$beta * 
                                                               gammahat * (n/k)^mopest$rho) * (r/rbet)^gammahat * 
               (1 + ((1/gammahat - 1)^(-mopest$rho) * rbet^(-mopest$rho) - 
                       1)/mopest$rho * mopest$beta * gammahat * (1 - 
                                                                   tau)^(-mopest$rho))/(1 + ((1/gammahat - 1)^(-mopest$rho) * 
                                                                                               r^(-mopest$rho) - 1)/mopest$rho * mopest$beta * 
                                                                                          gammahat * (k/n)^(-mopest$rho))))
  }
  if (br == TRUE && method == "indirect") {
    return(c(gammahat,(1/gammahat - 1)^(-gammahat) * quantile(X, 1 - 
                                                                k/n) * (k/(n * (1 - tau)))^gammahat * (1 + ((k/(n * 
                                                                                                                  (1 - tau)))^mopest$rho - 1)/mopest$rho * mopest$beta * 
                                                                                                         gammahat * (n/k)^mopest$rho) * (1 + ((1/gammahat - 
                                                                                                                                                 1)^(-mopest$rho) * rbet^(-mopest$rho) - 1)/mopest$rho * 
                                                                                                                                           mopest$beta * gammahat * (1 - tau)^(-mopest$rho))/(rbet^gammahat)))
  }
}


extExpect.eva <- function(tau,train){
  # tail.ind <- lpindex(X=train,k=50,p=1.5,br=TRUE)
  # tail.ind <- tindexp(X=train, k=50,br=TRUE)
  expect <- extExpect.modify(X=train, k="kopt", tau=tau, estim = "Hill", method="direct", br = TRUE)
  tail.ind <- expect[1]
  res <- expect[2] * (1/tail.ind -1)^tail.ind
  return(res)
}


extQuant.eva <- function(tau, train){
  return(extQuant(X=train, tau=tau, k=10, estim="tindexp", br=TRUE))
}


extQuantlp.eva <- function(tau, train){
  return(extQuantlp(X=train,tau=tau,k=10,p=1.2,estim="lpindex",br=TRUE))
}

extremefit.eva <- function(tau, train){
  HH <- hill.adapt(train, weights=rep(1, length(train)), initprop = 0.1,
                   gridlen = 100 , r1 = 0.25, r2 = 0.05, CritVal=10)
  
  return(predict(HH, tau, type = "quantile")$y)
}
################
source("evaluation.R")

data <- read.csv("../data/Amaurot.csv")
data[which(data$Season=="S1"),"Season"] <- 0
data[which(data$Season=="S2"),"Season"] <- 1
data$Season <- as.numeric(data$Season)

data2 <- data[complete.cases(data),] # remove NA's 
data3 <- data2[-which.max(data2$Y),] # remove max


# Weissman-type estimator
calculate_SCV1(extQuant.eva, tau0=1-1/60000, data2$Y, seed=1) 
calculate_SCV1(extQuant.eva, tau0=1-1/60000, data2$Y, seed=100) 
calculate_SCV1(extQuant.eva, tau0=1-1/60000, data2$Y, seed=200)
calculate_SCV1(extQuant.eva, tau0=1-1/60000, data2$Y, seed=300) 

calculate_SCV2(extQuant.eva, tau0=1-1/60000, data2$Y, seed=1) 
calculate_SCV2(extQuant.eva, tau0=1-1/60000, data2$Y, seed=100)
calculate_SCV2(extQuant.eva, tau0=1-1/60000, data2$Y, seed=200) 
calculate_SCV2(extQuant.eva, tau0=1-1/60000, data2$Y, seed=300) 

calculate_SCV1(extQuant.eva, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV1(extQuant.eva, tau0=1-1/60000, data3$Y, seed=100) 
calculate_SCV1(extQuant.eva, tau0=1-1/60000, data3$Y, seed=200) 
calculate_SCV1(extQuant.eva, tau0=1-1/60000, data3$Y, seed=300)

calculate_SCV2(extQuant.eva, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(extQuant.eva, tau0=1-1/60000, data3$Y, seed=100)
calculate_SCV2(extQuant.eva, tau0=1-1/60000, data3$Y, seed=200) 
calculate_SCV2(extQuant.eva, tau0=1-1/60000, data3$Y, seed=300)

# Composite extreme quantile estimation
calculate_SCV1(extQuantlp.eva, tau0=1-1/60000, data2$Y, seed=1)
calculate_SCV1(extQuantlp.eva, tau0=1-1/60000, data2$Y, seed=100)
calculate_SCV1(extQuantlp.eva, tau0=1-1/60000, data2$Y, seed=200) 
calculate_SCV1(extQuantlp.eva, tau0=1-1/60000, data2$Y, seed=300) 

calculate_SCV2(extQuantlp.eva, tau0=1-1/60000, data2$Y, seed=1) 
calculate_SCV2(extQuantlp.eva, tau0=1-1/60000, data2$Y, seed=100) 
calculate_SCV2(extQuantlp.eva, tau0=1-1/60000, data2$Y, seed=200) 
calculate_SCV2(extQuantlp.eva, tau0=1-1/60000, data2$Y, seed=300) 

calculate_SCV1(extQuantlp.eva, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV1(extQuantlp.eva, tau0=1-1/60000, data3$Y, seed=100) 
calculate_SCV1(extQuantlp.eva, tau0=1-1/60000, data3$Y, seed=200) 
calculate_SCV1(extQuantlp.eva, tau0=1-1/60000, data3$Y, seed=300)

calculate_SCV2(extQuantlp.eva, tau0=1-1/60000, data3$Y, seed=1)
calculate_SCV2(extQuantlp.eva, tau0=1-1/60000, data3$Y, seed=100) 
calculate_SCV2(extQuantlp.eva, tau0=1-1/60000, data3$Y, seed=200) 
calculate_SCV2(extQuantlp.eva, tau0=1-1/60000, data3$Y, seed=300) 

# expectile
calculate_SCV1(extExpect.eva, tau0=1-1/60000, data2$Y, seed=1)
calculate_SCV1(extExpect.eva, tau0=1-1/60000, data2$Y, seed=100)
calculate_SCV1(extExpect.eva, tau0=1-1/60000, data2$Y, seed=200) 
calculate_SCV1(extExpect.eva, tau0=1-1/60000, data2$Y, seed=300) 

calculate_SCV2(extExpect.eva, tau0=1-1/60000, data2$Y, seed=1)
calculate_SCV2(extExpect.eva, tau0=1-1/60000, data2$Y, seed=100) 
calculate_SCV2(extExpect.eva, tau0=1-1/60000, data2$Y, seed=200) 
calculate_SCV2(extExpect.eva, tau0=1-1/60000, data2$Y, seed=300) 

calculate_SCV1(extExpect.eva, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV1(extExpect.eva, tau0=1-1/60000, data3$Y, seed=100) 
calculate_SCV1(extExpect.eva, tau0=1-1/60000, data3$Y, seed=200) 
calculate_SCV1(extExpect.eva, tau0=1-1/60000, data3$Y, seed=300) 

calculate_SCV2(extExpect.eva, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(extExpect.eva, tau0=1-1/60000, data3$Y, seed=100) 
calculate_SCV2(extExpect.eva, tau0=1-1/60000, data3$Y, seed=200) 
calculate_SCV2(extExpect.eva, tau0=1-1/60000, data3$Y, seed=300) 

# extremefit
calculate_SCV1(extremefit.eva, tau0=1-1/60000, data2$Y, seed=1) 
calculate_SCV1(extremefit.eva, tau0=1-1/60000, data2$Y, seed=100) 
calculate_SCV1(extremefit.eva, tau0=1-1/60000, data2$Y, seed=200) 
calculate_SCV1(extremefit.eva, tau0=1-1/60000, data2$Y, seed=300) 

calculate_SCV2(extremefit.eva, tau0=1-1/60000, data2$Y, seed=1) 
calculate_SCV2(extremefit.eva, tau0=1-1/60000, data2$Y, seed=100)
calculate_SCV2(extremefit.eva, tau0=1-1/60000, data2$Y, seed=200) 
calculate_SCV2(extremefit.eva, tau0=1-1/60000, data2$Y, seed=300)

calculate_SCV1(extremefit.eva, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV1(extremefit.eva, tau0=1-1/60000, data3$Y, seed=100) 
calculate_SCV1(extremefit.eva, tau0=1-1/60000, data3$Y, seed=200)
calculate_SCV1(extremefit.eva, tau0=1-1/60000, data3$Y, seed=300) 

calculate_SCV2(extremefit.eva, tau0=1-1/60000, data3$Y, seed=1) 
calculate_SCV2(extremefit.eva, tau0=1-1/60000, data3$Y, seed=100) 
calculate_SCV2(extremefit.eva, tau0=1-1/60000, data3$Y, seed=200) 
calculate_SCV2(extremefit.eva, tau0=1-1/60000, data3$Y, seed=300) 

# test
# Weissman-type estimator
extQuant(X=data2$Y, tau=1-1/60000) 
extQuant(X=data3$Y, tau=1-1/60000) 

# Composite extreme quantile estimation (selected, )
extQuantlp(X=data2$Y,tau=1-1/60000,k=50,p=1.5,estim="lpindex",br=TRUE) 
extQuantlp(X=data3$Y,tau=1-1/60000,k=50,p=1.5,estim="lpindex",br=TRUE) 

# expectile
extExpect.eva(tau=1-1/60000, train=data2$Y)
extExpect.eva(tau=1-1/60000, train=data3$Y) 

# extremefit
extremefit.eva(tau=1-1/60000, train=data2$Y) 
extremefit.eva(tau=1-1/60000, train=data3$Y) 


