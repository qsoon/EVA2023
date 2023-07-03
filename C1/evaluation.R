library(caret)
library(evmix)

## marginal quantile estimation (C2) ##

qtloss <- function(predictor, tau, data){
  return((data-predictor)*(tau-as.numeric(data-predictor<0)))
}

SCV <- function(predfunc, tau, train, valid){
  predictor <- predfunc(tau, train) # predfunc input = (tau, data)
  return(mean(qtloss(predictor, tau, valid)))
}

calculate_QS <- function(predfunc, tau0, y){
  return(SCV(predfunc, tau0, y, y))
}

# method 1
calculate_SCV1 <- function(predfunc, tau0=1-1/60000, y, seed=1){
  n <- length(y)
  res <- 0
  alpha_list <- 2**c(0:floor(log(n^0.25,2)))
  l <- length(alpha_list)
  
  for(alpha in alpha_list){
    tauc <- tau0 - alpha / n
    nc <- n / (1 + alpha/n/(1-tau0))
    k <- floor(n/nc)
    if(k==1){
      l <- l-1
      next
    }
    scvsum <- 0
    set.seed(seed)
    flds <- createFolds(1:n, k = k, list = TRUE, returnTrain = FALSE)
    K <- k
    for(j in 1:k){
      train_ind <- flds[[j]]
      train <- y[train_ind]
      valid <- y[-train_ind]
      tmp <- SCV(predfunc, tauc, train, valid)
      if(is.na(tmp)){
        K <- K-1
        next
      }
      scvsum <- scvsum + tmp
    }
    if(K==0){
      l <- l-1
      next
    }
    res <- res + scvsum/K
  }
  return(res/l)
}

c2loss <- function(predfunc, tau, train, valid){
  predictor <- predfunc(tau, train) # predfunc input = (tau, data)
  if(0.99*quantile(valid, probs=tau)>predictor){
    return(0.9*(0.99*quantile(valid, probs=tau) - predictor))
  } 
  else if(1.01*quantile(valid, probs=tau)<predictor){
    return(0.1*(predictor-1.01*quantile(valid, probs=tau)))
  }
  else{
    return(0)
  }
}


calculate_loss1 <- function(predfunc, tau0=1-1/60000, y, seed=1){
  n <- length(y)
  res <- 0
  alpha_list <- 2**c(0:floor(log(n^0.25,2)))
  l <- length(alpha_list)
  
  scvlist <- list()
  k_list <- c()
  ind=1
  for(alpha in alpha_list){
    tauc <- tau0 - alpha / n
    nc <- n / (1 + alpha/n/(1-tau0))
    k <- floor(n/nc)
    if(k==1){
      l <- l-1
      next
    }
    scvsum <- 0
    set.seed(seed)
    flds <- createFolds(1:n, k = k, list = TRUE, returnTrain = FALSE)
    K <- k
    if(k >= 4){
      k_list <- c(k_list, k)
      for(j in 1:k){
        train_ind <- flds[[j]]
        train <- y[train_ind]
        valid <- y[-train_ind]
        tmp <- c2loss(predfunc, tauc, train, valid)
        if(is.na(tmp)){
          K <- K-1
          next
        }
        scvsum <- scvsum + tmp
      }
      if(K==0){
        l <- l-1
        next
      }
      scvlist[[ind]]<- scvsum/K
      ind <- ind + 1
    }
  }
  names(scvlist) <- k_list
  return(scvlist)
}

# method 2
calculate_SCV2 <- function(predfunc, tau0=1-1/60000, y, seed=1){
  n <- length(y)
  res <- 0
  alpha_list <- 2**c(0:floor(log(n^0.25,2)))
  l <- length(alpha_list)
  
  for(alpha in alpha_list){
    tauc1 <- tau0 - alpha / n
    nc <- n / (1 + alpha/n/(1-tau0))
    k <- floor(n/nc)
    
    if(k==1){
      l <- l-1
      next
    }
    alpha2 <- n*(1-tau0)/(k-1)
    tauc <- tau0 - alpha2/n
    
    scvsum <- 0
    set.seed(seed)
    flds <- createFolds(1:n, k = k, list = TRUE, returnTrain = FALSE)
    K <- k
    for(j in 1:k){
      valid_ind <- flds[[j]]
      train <- y[-valid_ind]
      valid <- y[valid_ind]
      tmp <- SCV(predfunc, tauc, train, valid)
      if(is.na(tmp)){
        K <- K-1
        next
      }
      scvsum <- scvsum + tmp
    }
    if(K==0){
      l <- l-1
      next
    }
    res <- res + scvsum/K
  }
  return(res/l)
}


## condition quantile estimation (C1) ##


cond_SCV <- function(predfunc_cond, tau, train_y, train_x, valid_y, valid_x){
  predictor <- predfunc_cond(tau, train_y, train_x, valid_x) 
  # predfunc input = (tau, train_y, train_x, valid_x), output = prediction (at valid_x)vector of size = length(valid set)
  
  return(mean(qtloss(predictor, tau, valid_y)))
}

cond_SCV_kernel <- function(predfunc_cond, tau, train_y, train_x, valid_y, valid_x, train_ind, valid_ind, distmat){
  predictor <- predfunc_cond(tau, train_y, train_x, train_x) 
  # predfunc input = (tau, train_y, train_x, valid_x), output = prediction (at valid_x)vector of size = length(valid set)
  
  res <- 0
  
  for(i in 1:nrow(train_x)){
    res <- res + mean(kdgaussian(distmat[train_ind[i], valid_ind], bw=5)*qtloss(predictor[i], tau, valid_y))
  }
  return(res / nrow(train_x))
}

# method 1
calculate_cond_SCV1 <- function(predfunc_cond, tau0=1-1/10000, y, X, seed=1){
  n <- length(y)
  res <- 0
  alpha_list <- 2**c(0:floor(log(n^0.25,2)))
  l <- length(alpha_list)
  
  for(alpha in alpha_list){
    tauc <- tau0 - alpha / n
    nc <- n / (1 + alpha/n/(1-tau0))
    k <- floor(n/nc)
    if(k==1){
      l <- l-1
      next
    }
    scvsum <- 0
    
    set.seed(seed)
    flds <- createFolds(1:n, k = k, list = TRUE, returnTrain = FALSE)
    K <- k
    for(j in 1:k){
      train_ind <- flds[[j]]
      train_y <- y[train_ind]
      train_x <- X[train_ind,]
      valid_y <- y[-train_ind]
      valid_x <- X[-train_ind,]
      tmp <- cond_SCV(predfunc_cond, tauc, train_y, train_x, valid_y, valid_x)
      if(is.na(tmp)){
        K <- K-1
        next
      }
      scvsum <- scvsum + tmp
    }
    if(K==0){
      l <- l-1
      next
    }
    res <- res + scvsum/K
  }
  return(res/l)
}

calculate_cond_SCV1_kernel <- function(predfunc_cond, tau0=1-1/10000, y, X, distmat, seed=1){
  n <- length(y)
  res <- 0
  alpha_list <- 2**c(0:floor(log(n^0.25,2)))
  l <- length(alpha_list)
  
  for(alpha in alpha_list){
    tauc <- tau0 - alpha / n
    nc <- n / (1 + alpha/n/(1-tau0))
    k <- floor(n/nc)
    if(k==1){
      l <- l-1
      next
    }
    scvsum <- 0
    
    set.seed(seed)
    flds <- createFolds(1:n, k = k, list = TRUE, returnTrain = FALSE)
    K <- k
    for(j in 1:k){
      train_ind <- flds[[j]]
      valid_ind <- c(1:n)[-train_ind]
      train_y <- y[train_ind]
      train_x <- X[train_ind,]
      valid_y <- y[-train_ind]
      valid_x <- X[-train_ind,]
      tmp <- cond_SCV_kernel(predfunc_cond, tauc, train_y, train_x, valid_y, valid_x, train_ind, valid_ind, distmat)
      if(is.na(tmp)){
        K <- K-1
        next
      }
      scvsum <- scvsum + tmp
    }
    if(K==0){
      l <- l-1
      next
    }
    res <- res + scvsum/K
  }
  return(res/l)
}


# method 2
calculate_cond_SCV2 <- function(predfunc_cond, tau0=1-1/10000, y, X, seed=1){
  n <- length(y)
  res <- 0
  alpha_list <- 2**c(0:floor(log(n^0.25,2)))
  l <- length(alpha_list)
  
  for(alpha in alpha_list){
    tauc1 <- tau0 - alpha / n
    nc <- n / (1 + alpha/n/(1-tau0))
    k <- floor(n/nc)
    
    if(k==1){
      l <- l-1
      next
    }
    alpha2 <- n*(1-tau0)/(k-1)
    tauc <- tau0 - alpha2/n
    
    scvsum <- 0
    set.seed(seed)
    flds <- createFolds(1:n, k = k, list = TRUE, returnTrain = FALSE)
    K <- k
    for(j in 1:k){
      valid_ind <- flds[[j]]
      train_y <- y[-valid_ind]
      train_x <- X[-valid_ind,]
      valid_y <- y[valid_ind]
      valid_x <- X[valid_ind,]
      tmp <- cond_SCV(predfunc_cond, tauc, train_y, train_x, valid_y, valid_x)
      if(is.na(tmp)){
        K <- K-1
        next
      }
      scvsum <- scvsum + tmp
    }
    res <- res + scvsum/K
  }
  return(res/l)
}

calculate_cond_SCV2_kernel <- function(predfunc_cond, tau0=1-1/10000, y, X, distmat, seed=1){
  n <- length(y)
  res <- 0
  alpha_list <- 2**c(0:floor(log(n^0.25,2)))
  l <- length(alpha_list)
  
  for(alpha in alpha_list){
    tauc1 <- tau0 - alpha / n
    nc <- n / (1 + alpha/n/(1-tau0))
    k <- floor(n/nc)
    
    if(k==1){
      l <- l-1
      next
    }
    alpha2 <- n*(1-tau0)/(k-1)
    tauc <- tau0 - alpha2/n
    
    scvsum <- 0
    set.seed(seed)
    flds <- createFolds(1:n, k = k, list = TRUE, returnTrain = FALSE)
    K <- k
    for(j in 1:k){
      valid_ind <- flds[[j]]
      train_ind <- c(1:n)[-valid_ind]
      train_y <- y[-valid_ind]
      train_x <- X[-valid_ind,]
      valid_y <- y[valid_ind]
      valid_x <- X[valid_ind,]
      tmp <- cond_SCV_kernel(predfunc_cond, tauc, train_y, train_x, valid_y, valid_x, train_ind, valid_ind, distmat)
      if(is.na(tmp)){
        K <- K-1
        next
      }
      scvsum <- scvsum + tmp
    }
    res <- res + scvsum/K
  }
  return(res/l)
}

