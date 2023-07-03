library(evgam)

# selected function forms
f1 <- formula(train_y~V2+V3*Season+WindSpeed+Atmosphere)
f2 <- formula(~ Season)
f3 <- formula(excess~V2+V3*Season+WindSpeed+Atmosphere)
f4 <- formula(~ Season)

tau.ald <- 0.9

get.est<-function(tau,train_y, train_x, valid_x){
  train<-as.data.frame(cbind(train_y,train_x))
  newdata<-as.data.frame(valid_x)
  
  #Step 1. ALD distribution fitting for the estimation of the threshold u(x)
  #list: GAM form for the estimation of 1st parameter: location parameter u(x), 2nd parameter: shape parameter \sigma(x)
  form1<-list(f1, f2 )
  ald<-evgam(form1, data=train, family="ald", ald.args=list(tau=tau.ald))
  
  # save y(x)-\hat{u(x)} as excess & make non-positive elements NA
  train$excess<-train_y-predict(ald)$location
  is.na(train$excess[train$excess<=0])<-TRUE
  
  #Step 2. Estimate the other two parameters in GPD(u(x),\psi(x),\xi(x)) 
  #list: GAM form for the estimation of 1st: location parameter \psi(x), 2nd: shape parameter \xi(x) 
  form2<-list(f3 , f4 )
  mod<-evgam(form2, data=train, family="gpd")
  
  #Step 3. Finally compute F^{-1}(0.9999)+\hat{u(x)} for new covariate x
  return(unlist(predict(mod,newdata=newdata,type="quantile",prob=.9999)[,1]+predict(ald, newdata = newdata)$location)) 
}


test <- read.csv("../data/AmaurotTestSet.csv", stringsAsFactors = TRUE)
data <- read.csv("../data/Amaurot.csv", stringsAsFactors = TRUE)
data2 <- data[complete.cases(data),] # remove NA's 

#### CI with bootstrap
n<-nrow(data2)
B<-2000
res<-list()
set.seed(1)
for(i in 1:B){
  boot.id = sample(1:n, n, replace=TRUE)
  boot.data = data2[boot.id,]
  res[[i]]<-get.est(tau=0.9999, train_y=boot.data[,1], train_x=boot.data[,-1], valid_x=test)
  if(i%%100==0){print(i)}
}
mat<-do.call(rbind, res)
lbound<-apply(mat,2,function(x){quantile(x, 0.25)})
ubound<-apply(mat,2,function(x){quantile(x, 0.75)})
est<-get.est(tau=0.9999, train_y=data2[,1], train_x=data2[,-1], valid_x=test)
result<-as.data.frame(cbind(lbound, est, ubound))

result<-round(result,3)

AnswerC1 = matrix(0,nrow=100, ncol=3)
AnswerC1[,1] = result[,2]
AnswerC1[,2] = result[,1]
AnswerC1[,3] = result[,3]

save(AnswerC1, file = "AnswerC1.RData")