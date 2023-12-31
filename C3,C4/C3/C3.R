library(texmex)
library(gridExtra)
library(latex2exp)
source("../texmex-funs.R") # replace the marginal estimating part, as we know the exact marg. distn.

c3 <- read.csv("../../data/Coputopia.csv")
c3slice6 <- list(c3[c3$Atmosphere <= quantile(c3$Atmosphere,0.2),],
                 c3[c3$Atmosphere > quantile(c3$Atmosphere,0.2) & c3$Atmosphere <= quantile(c3$Atmosphere,0.4),],
                 c3[c3$Atmosphere > quantile(c3$Atmosphere,0.4) & c3$Atmosphere <= quantile(c3$Atmosphere,0.6),],
                 c3[c3$Atmosphere > quantile(c3$Atmosphere,0.6) & c3$Atmosphere <= quantile(c3$Atmosphere,0.8),],
                 c3[c3$Atmosphere > quantile(c3$Atmosphere,0.8) & c3$Atmosphere <= quantile(c3$Atmosphere,0.95),],
                 c3[c3$Atmosphere > quantile(c3$Atmosphere,0.95) & c3$Atmosphere <= quantile(c3$Atmosphere,1.0),])


### c3 - p1
g <- 1  # vary g = 1,...,6
data <- c3slice6[[g]]
marg <- migpd(data[,3:5], mqu=.7, penalty="none")
sapply(c3slice6,\(ll)dim(ll)[1])

for(dqu.candidate in seq(0.7,0.99,by=0.01)){
  mex.Y1 <- mexDependence(marg, which="Y1", dqu=dqu.candidate)
  mex.Y2 <- mexDependence(marg, which="Y2", dqu=dqu.candidate)
  mex.Y3 <- mexDependence(marg, which="Y3", dqu=dqu.candidate)
  
  set.seed(0)
  nsim <- 1e6
  pqu.common <- exp(-exp(-6))
  pY1 <- predict(mex.Y1, pqu=pqu.common, nsim=nsim, smoothZdistribution=T)
  pY2 <- predict(mex.Y2, pqu=pqu.common, nsim=nsim, smoothZdistribution=T)
  pY3 <- predict(mex.Y3, pqu=pqu.common, nsim=nsim, smoothZdistribution=T)
  
  marg.pr1 <- 1-exp(-exp(-6))
  res <- mean(pY1$data$CondLargest & apply(pY1$data$simulated[,2:3],1,min) > 6) * marg.pr1 +
    mean(pY2$data$CondLargest & apply(pY2$data$simulated[,2:3],1,min) > 6) * marg.pr1 +
    mean(pY3$data$CondLargest & apply(pY3$data$simulated[,2:3],1,min) > 6) * marg.pr1
  print(paste0(dqu.candidate," / ",res))
} 


### c3 - p2
g <- 1 # vary g = 1,...,6
data <- c3slice6[[g]] 
marg <- migpd(data[,3:5], mqu=.7, penalty="none")

for(dqu.candidate in seq(0.7,0.99,by=0.01)){
  mex.Y1 <- mexDependence(marg, which="Y1", dqu=dqu.candidate)
  mex.Y2 <- mexDependence(marg, which="Y2", dqu=dqu.candidate)
  
  set.seed(0)
  nsim <- 1e6
  pqu.common <- exp(-exp(-7))
  pY1b <- predict(mex.Y1, pqu=pqu.common, nsim=nsim, smoothZdistribution=T)
  pY2b <- predict(mex.Y2, pqu=pqu.common, nsim=nsim, smoothZdistribution=T)
  
  marg.pr2 <- 1-exp(-exp(-7))
  m <- -log(log(2))
  res <- mean(pY1b$data$simulated[,1]>pY1b$data$simulated[,2] & pY1b$data$simulated[,2]>7 & pY1b$data$simulated[,3] < m) * marg.pr2 +
    mean(pY2b$data$simulated[,1]>pY2b$data$simulated[,2] & pY2b$data$simulated[,2]>7 & pY2b$data$simulated[,3] < m) * marg.pr2
  print(paste0(dqu.candidate," / ",res))
} 


################ plot ##########################################################
# results of running the code above are recorded in `c3-results.xlsx` file
library(readxl)
res.i <- read_excel("c3-results.xlsx",
                    sheet = "1e6", range = "B1:I31") # P(Y1,Y2,Y3>6)
res.ii <- read_excel("c3-results.xlsx",
                     sheet = "1e6", range = "B34:I64") # P(Y1,Y2>7,Y3<m)

colors <- c("black","red","green","magenta","cyan","blue","orange")
us <- seq(0.7,0.99,by=0.01)

#### plot --- c3(i)
par(mfrow=c(1,2),oma=c(0,0,0,9.5),mar=c(5,4,4,1)+0.1) 

plot(us,res.i$`0~20%`,
     type="l",xlab=TeX("$italic(u)$"),
     ylab="Probability",main=TeX("(a) $\\widehat{Pr}(italic(Y_1>6,Y_2>6,Y_3>6) | Atmosphere\\in O_g)$"),
     ylim=c(0.0,0.0003),col=colors[1]) 
lines(us,res.i$`20~40%`,col=colors[2])
lines(us,res.i$`40~60%`,col=colors[3])
lines(us,res.i$`60~80%`,col=colors[4])
lines(us,res.i$`80~95%`,col=colors[5])
lines(us,res.i$`95~100%`,col=colors[6])
lines(us,res.i$newcombine,lwd=3)  # combine 6 groups
# legend("topleft",legend=c("Atmosphere","0~20%","20~40%","40~60%","60~80%","80~95%","95~100%","weighted average"),
#        col=c("",colors[1:6],"black"),
#        lty=c(0,rep(1,7)),lwd=c(0,rep(1,6),3))

### plot --- c3(ii)
plot(us,res.ii$`0~20%`,
     type="l",xlab=TeX("$italic(u)$"),ylab="Probability",
     main=TeX("(a) $\\widehat{Pr}(italic(Y_1>7,Y_2>7,Y_3<m) | Atmosphere\\in O_g)$"),
     ylim=c(0.0,4e-5),col=colors[1])
lines(us,res.ii$`20~40%`,col=colors[2])
lines(us,res.ii$`40~60%`,col=colors[3])
lines(us,res.ii$`60~80%`,col=colors[4])
lines(us,res.ii$`80~95%`,col=colors[5])
lines(us,res.ii$`95~100%`,col=colors[6])
lines(us,res.ii$newcombine,lwd=3,col="black")
legend("topright",inset=c(-0.55,0),xpd=NA,bty="n",
       legend=c("Atmosphere","0~20% (g=1)","20~40% (g=2)","40~60% (g=3)","60~80% (g=4)","80~95% (g=5)","95~100% (g=6)","weighted average"),
       col=c("",colors[1:6],"black"),
       lty=c(0,rep(1,7)),lwd=c(0,rep(1,6),3))
