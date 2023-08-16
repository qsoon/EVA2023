### EDA of the sub-challenge C3

library(texmex)
library(gridExtra)

c3 <- read.csv("../../data/Coputopia.csv")


### C3 in uniform margin space

par(mfrow=c(1,3))
plot(exp(-exp(-c3$Y1)),exp(-exp(-c3$Y2)),pch=".",main="(a) F(Y1) vs. F(Y2)",asp=1, xlab="F(Y1)",ylab="F(Y2)")
plot(exp(-exp(-c3$Y1)),exp(-exp(-c3$Y3)),pch=".",main="(b) F(Y1) vs. F(Y3)",asp=1, xlab="F(Y1)",ylab="F(Y3)")
plot(exp(-exp(-c3$Y2)),exp(-exp(-c3$Y3)),pch=".",main="(c) F(Y2) vs. F(Y3)",asp=1, xlab="F(Y2)",ylab="F(Y3)")


### Empirical and theoretic probability (under independence assumption)

par(mfrow=c(1,2))
plot(seq(2,10,by=0.1),
     sapply(seq(2,10,by=0.1),\(r)mean(rowSums(c3[,c(3,4,5)]>r)==3)),type="l",
     xlab="r",ylab="Probability",main="(a) P(Y1>r,Y2>r,Y3>r)",
     xlim=c(2,10),ylim=c(0.0,0.01))
lines(seq(2,10,by=0.1),
      sapply(seq(2,10,by=0.1),\(r)(1-exp(-exp(-r)))^3),
      col="red")

m <- -log(log(2))
plot(seq(2,10,by=0.1),
     sapply(seq(2,10,by=0.1),\(r)mean(rowSums(c3[,c(3,4)]>r)==2 & c3[,5]<m)),type="l",
     xlab="r",ylab="Probability",main="(b) P(Y1>r,Y2>r,Y3<m)",
     xlim=c(2,10),ylim=c(0.0,0.01))
lines(seq(2,10,by=0.1),
      sapply(seq(2,10,by=0.1),\(r)(1-exp(-exp(-r)))^2/2),
      col="red")
legend("topright",c("Empirical probability", "Independence model"),col=c("black","red"),lty=1)

### EDA on the covariates
################ Atmosphere ####################################################

colors <- c("black","red","green","magenta","cyan","blue","orange")

c3slice6 <- list(c3[c3$Atmosphere <= quantile(c3$Atmosphere,0.2),],
                 c3[c3$Atmosphere > quantile(c3$Atmosphere,0.2) & c3$Atmosphere <= quantile(c3$Atmosphere,0.4),],
                 c3[c3$Atmosphere > quantile(c3$Atmosphere,0.4) & c3$Atmosphere <= quantile(c3$Atmosphere,0.6),],
                 c3[c3$Atmosphere > quantile(c3$Atmosphere,0.6) & c3$Atmosphere <= quantile(c3$Atmosphere,0.8),],
                 c3[c3$Atmosphere > quantile(c3$Atmosphere,0.8) & c3$Atmosphere <= quantile(c3$Atmosphere,0.95),],
                 c3[c3$Atmosphere > quantile(c3$Atmosphere,0.95) & c3$Atmosphere <= quantile(c3$Atmosphere,1.0),])

par(mfrow=c(1,2))

#### P(Y1,Y2,Y3 > r): 6 slice
plot(seq(-2,10,by=0.1),
     sapply(seq(-2,10,by=0.1),\(r)mean(rowSums(c3slice6[[1]][,c(3,4,5)]>r)==3)),type="l",
     xlab="r",ylab="Empirical probability",main="(a) P(Y1>r,Y2>r,Y3>r)",
     xlim=c(1,10),ylim=c(0.0,0.1))
for(i in 2:length(c3slice6)){
  lines(seq(-2,10,by=0.1),
        sapply(seq(-2,10,by=0.1),\(r)mean(rowSums(c3slice6[[i]][,c(3,4,5)]>r)==3)),
        col=rep(colors,each=1)[i])
  Sys.sleep(0.1)
}
legend("topright",legend=c("Atmosphere","0~20%","20~40%","40~60%","60~80%","80~95%","95~100%"),
       col=c("",colors[1:6]),lty=c(0,rep(1,6)),lwd=c(0,rep(2,6)))

################ Season ########################################################

#### P(Y1>r,Y2>r,Y3>r)
plot(seq(-2,10,by=0.1),
     sapply(seq(-2,10,by=0.1),\(r)mean(rowSums(c3[c3$Season=="S1",c(3,4,5)]>r)==3)),type="l",
     xlab="r",ylab="Empirical probability",main="(b) P(Y1>r,Y2>r,Y3>r)",
     xlim=c(1,10),ylim=c(0.0,0.1))
lines(seq(-2,10,by=0.1),
      sapply(seq(-2,10,by=0.1),\(r)mean(rowSums(c3[c3$Season=="S2",c(3,4,5)]>r)==3)), col="red")
legend("topright",legend=c("Season","S1","S2"),
       col=c("","black","red"),lty=c(0,rep(1,2)),lwd=c(0,rep(2,2)))