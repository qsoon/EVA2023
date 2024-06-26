library(matrixStats)
library(igraph)
library(texmex)
library(tidyverse)

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
spearman <- cor(data, method = 'spearman')

#hist(abs(as.vector(pearson)[!abs(as.vector(pearson))==1]), main = "Pearson's rho", xlab = "rho")

# get only elements above the main diagonal to avoid duplications
hist(abs(as.vector(pearson))[which((1:2500)%%50>(1:2500)%/%50+1|((1:2500)%%50==0)&(1:2500)!=2500)], main = "Pearson's rho", xlab = "rho")
hist(abs(as.vector(spearman))[which((1:2500)%%50>(1:2500)%/%50+1|((1:2500)%%50==0)&(1:2500)!=2500)], main = "Spearman's rho", xlab = "rho")

utopula <- cbind(data1, data2)

chi_pair<- matrix(0, nrow = 50, ncol = 50)
for(i in 1:49){
  for(j in (i+1):50){
    chi_mat <- texmex::chi(cbind(utopula[,i], utopula[,j]))$chi
    chi_pair[i,j] <- chi_mat[100,2]
  }
}

chi_pair[upper.tri(chi_pair)] %>%
  hist(main = "Extremal dependence")

################################################################################

pearsonedges <- c()
spearmanedges <- c()

for (i in 1:50) {
  for (j in 1:50) {
    if (pearson[i,j]>0.2 & i!=j) {
      pearsonedges <- c(pearsonedges,i,j)
    }
    if (spearman[i,j]>0.2 & i!=j) {
      spearmanedges <- c(spearmanedges,i,j)
    }
  }
}

pearsongraph <- graph(pearsonedges, directed=FALSE)
spearmangraph <- graph(spearmanedges, directed=FALSE)

plot(pearsongraph,vertex.size=8)
plot(spearmangraph,vertex.size=8)




