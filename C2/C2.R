# For C2, we selected Lp composite extreme quantile estimation method

data <- read.csv("../data/Amaurot.csv")

wrap.extQuantlp <- function(tau, y, nk, p, method="lpindex", is.br=FALSE) {
  n = length(y)
  k = round(n/nk)
  res = extQuantlp(X=y, tau=tau, k=k, p=p, estim=method, br=is.br)
  res
}


AnswerC2 = wrap.extQuantlp(tau=1-1/60000, y=data$Y, nk=200, p=1.1, method="lpindex", is.br=TRUE)
save(AnswerC2, file = "AnswerC2.RData")