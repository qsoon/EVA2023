library(texmex)
library(matrixStats)
library(igraph)

source("../texmex-funs.R")
data1<-read.csv("../../data/UtopulaU1.csv")
data2<-read.csv("../../data/UtopulaU2.csv")

for (i in 1:25) {
  colnames(data1)[i] <- paste("1", colnames(data1)[i], sep="")
}
for (i in 1:25) {
  colnames(data2)[i] <- paste("2", colnames(data2)[i], sep="")
}

data<-cbind(data1,data2)

pearson <- cor(data, method = 'pearson')

hist(abs(as.vector(pearson)[!abs(as.vector(pearson))==1]), main = "Pearson's rho", xlab = "rho")

pearsonedges <- (pearson>0.2)

pearsonedges <- c()

for (i in 1:50) {
  for (j in 1:50) {
    if (pearson[i,j]>0.2 & i!=j) {
      pearsonedges <- c(pearsonedges,i,j)
    }
  }
}

pearsongraph <- graph(pearsonedges, directed=FALSE)

plot(pearsongraph,vertex.size=8)



k<-5

cl<-vector(mode='list',length=k)

cl[[1]]<-c(29,15,47,3,18,22,45,10)
cl[[2]]<-c(21,42,50,40,32,41,8,25,34,48,33,26,49)
cl[[3]]<-c(17,46,7,1,9,20,2,5,31)
cl[[4]]<-c(11,23,35,39,36,16,27,13,24,6,12,37)
cl[[5]]<-c(30,14,19,28,38,43,4,44)

phi1<-1/300
phi2<-12*phi1
#phi2<-phi1

s1<--log(-log(1-phi1))
s2<--log(-log(1-phi2))

datas<-vector(mode='list',length=k)
one<-vector(mode='list',length=k)
nsim<-10000000
p<-vector(mode='list', length=k)
comb<-vector(mode='list', length=k)
each<-vector(mode='list', length=k)
params<-c(0.8)

for (i in 1:k) {
  datas[[i]]<-data[,cl[[i]]]
  one<-matrix(rep(cl[[i]]<=25,10000),nrow=10000,byrow=TRUE)
  #print(sum(rowProds(as.matrix(datas[[i]]>(s1*one+s2*(1-one))))/10000))
  marg<-migpd(datas[[i]],mqu=0.7,penalty="none")
  pqu<-(1-phi1)*(cl[[i]]<=25)+(1-phi2)*(cl[[i]]>25)
  comb[[i]]<-vector(length=length(params))
  each[[i]]<-vector(mode='list',length=ncol(datas[[i]]))
  for (j in 1:ncol(datas[[i]])) {
    each[[i]][[j]]<-vector(length=length(params))
  }
  for (u in 1:length(params)) {
    models<-vector(mode='list',length=ncol(datas[[i]]))
    p[[i]]<-vector(length=ncol(datas[[i]]))
    for (j in 1:ncol(datas[[i]])) {
      models<-mexDependence(marg,which=colnames(datas[[i]])[j],dqu=params[u])
      sims<-predict(models,pqu=pqu[j],nsim=nsim,smoothZdistribution=T)
      temp<-sims$data$simulated
      one<-matrix(rep(substr(colnames(temp),1,1)=="1",nsim),nrow=nsim,byrow=TRUE)
      p[[i]][j]<-mean(rowProds(temp>(s1*one+s2*(1-one)))*sims$data$CondLargest)
      each[[i]][[j]][u] <- mean(rowProds(temp>(s1*one+s2*(1-one))))*(1-pqu[j]) # self-consistency
    }
    comb[[i]][u]<-sum(p[[i]]*(phi1*(cl[[i]]<=25)+phi2*(cl[[i]]>25)))
  }
  min <- min(comb[[i]])
  max <- max(comb[[i]])
  for (j in 1:ncol(datas[[i]])) {
    min <- min(c(min,each[[i]][[j]]))
    max <- max(c(max,each[[i]][[j]]))
  }
}

AnswerC4 = comb[[1]]*comb[[2]]*comb[[3]]*comb[[4]]*comb[[5]]
save(AnswerC4, file = "AnswerC4.RData")