rm(list=ls())
library(quantreg)
library(plm)
library(tidyverse)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
library(moments)

data<-read.csv("/Users/data_invest.csv", header=T)
attach(data)

cl <- makeCluster(8)
registerDoParallel(cl)


# Number quantiles

taus<-seq(0.03,0.97,0.001)
nquantiles<-length(taus)


#FE Canay
#Regression of Tobin's Q on Investment (Q1 on KexpK1)

Zaux <- model.matrix(~as.factor(id)-1)
Z<-Zaux[,-1]
ncolZ<-ncol(Z)

fitols<-summary(lm(Q1 ~ KexpK1+Z))
Q1hat<-Q1-(Z%*%fitols$coef[3:(ncolZ+2)])

score1.q1=matrix(0,nquantiles)
score2.q1=matrix(0,nquantiles)
mu.q1=matrix(0,nquantiles)
mu.se.q1=matrix(0,nquantiles)
alpha.q1=matrix(0,nquantiles)
alpha.se.q1=matrix(0,nquantiles)


for(i in 1:nquantiles){
  
  fit<-summary(rq(Q1hat ~ KexpK1, tau=taus[i]),se="ker")
  
  mu.q1[i]=fit$coef[1]
  alpha.q1[i]=fit$coef[2]
  mu.se.q1[i]=fit$coef[3]
  alpha.se.q1[i]=fit$coef[4]
  t=taus[i] 
  sign=(fit$residuals<0)
  rho=(fit$residuals)*(t-sign)
  score1.q1[i]= (mean(fit$residuals))/mean(rho) 
  score2.q1[i]=(1-2*t)/(t*(1-t))
}

score.q1=abs(score2.q1-score1.q1)
ptau.q1=which.min(score.q1)

tau.hat.q1=taus[ptau.q1]
mu.hat.q1=mu.q1[ptau.q1]
mu.hat.se.q1=mu.se.q1[ptau.q1]
alpha.hat.q1=alpha.q1[ptau.q1]
alpha.hat.se.q1=alpha.se.q1[ptau.q1]

plot(abs(score2.q1-score1.q1))

#Results Q1 on KexpK1
tau.hat.q1
cbind(alpha.hat.q1,alpha.hat.se.q1)

#Result Skewness of hat{u}
fit.final<-summary(rq(Q1hat ~ KexpK1, tau=taus[ptau.q1]),se="ker")
skewness(fit.final$resid)

## Inference tau
#############
## Subsampling

#the number of Subsampling

num.ss <- 199
tau.ss.dist <- array(0,dim=c(num.ss,1))

taus<-seq(0.03,0.97,0.01)
nquantiles<-length(taus)

n.id<-length(unique(id))
n.ss<-as.integer(n.id^(0.8))

for(j in 1:num.ss){

sub.ss <- sample.int(n.id, n.ss, replace = F)
data.ss <- data %>%filter(id %in% sub.ss)
data.ss<- transform(data.ss,id = as.numeric(factor(data.ss$id)))

Zaux.ss <- model.matrix(~as.factor(data.ss$id)-1)
Z.ss<-Zaux.ss[,-1]
ncolZ.ss<-ncol(Z.ss)

fitols.ss<-summary(lm(data.ss$Q1 ~ data.ss$KexpK1+Z.ss))
Q1hat.ss<-data.ss$Q1-(Z.ss%*%fitols.ss$coef[3:(ncolZ.ss+2)])

X.ss<-as.matrix.csr(cbind(1,data.ss$KexpK1,Z.ss))
Y.ss<-Q1hat.ss

tcoefs.ss.result <- foreach(k = 1:nquantiles, .combine='cbind', .packages = c("quantreg", "MASS")) %dopar% {

reg.ss <- rq.fit.sfn(X.ss, Y.ss, tau=taus[k])

  t=taus[k]
  sign=(reg.ss$residuals<0)
  rho=(reg.ss$residuals)*(t-sign)
  score1.q1= (mean(reg.ss$residuals))/mean(rho) 
  score2.q1=(1-2*t)/(t*(1-t))
  tcoefs.b<- abs(score2.q1-score1.q1)

return(tcoefs.b)

}

ptau.q1=which.min(tcoefs.ss.result)
tau.ss.dist[j]=taus[ptau.q1]

}

sd(tau.ss.dist)



#############
# CONSTRAINED
# Regression of Tobin's Q on Investment
#############

############
# Q1 on KexpK1 CONSTRAINED
# Constrained - Small ==1

taus<-seq(0.03,0.97,0.001)
nquantiles<-length(taus)

score1.q1=matrix(0,nquantiles)
score2.q1=matrix(0,nquantiles)
mu.q1=matrix(0,nquantiles)
mu.se.q1=matrix(0,nquantiles)
alpha.q1=matrix(0,nquantiles)
alpha.se.q1=matrix(0,nquantiles)

data.Small<- data %>%filter(data$Small==1)
data.Small<- transform(data.Small,id = as.numeric(factor(data.Small$id)))

Zaux.Small <- model.matrix(~as.factor(data.Small$id)-1)
Z.Small<-Zaux.Small[,-1]
ncolZ.Small<-ncol(Z.Small)

KexpK1.Small<-data.Small$KexpK1
Q1.Small<-data.Small$Q1

fitols.Small<-summary(lm(Q1.Small ~ KexpK1.Small+Z.Small))
Q1hat.Small<-Q1.Small-(Z.Small%*%fitols.Small$coef[3:(ncolZ.Small+2)])


for(i in 1:nquantiles){
  
  fit<-summary(rq(Q1hat.Small ~ KexpK1.Small, tau=taus[i]),se="ker")
  
  mu.q1[i]=fit$coef[1]
  alpha.q1[i]=fit$coef[2]
  mu.se.q1[i]=fit$coef[3]
  alpha.se.q1[i]=fit$coef[4]
  t=taus[i] 
  sign=(fit$residuals<0)
  rho=(fit$residuals)*(t-sign)
  score1.q1[i]= (mean(fit$residuals))/mean(rho) 
  score2.q1[i]=(1-2*t)/(t*(1-t))
}

score.q1=abs(score2.q1-score1.q1)
ptau.q1=which.min(score.q1)

tau.hat.q1=taus[ptau.q1]
mu.hat.q1=mu.q1[ptau.q1]
mu.hat.se.q1=mu.se.q1[ptau.q1]
alpha.hat.q1=alpha.q1[ptau.q1]
alpha.hat.se.q1=alpha.se.q1[ptau.q1]

plot(abs(score2.q1-score1.q1))

#Results Q1 on KexpK1 Small==1
tau.hat.q1
cbind(alpha.hat.q1,alpha.hat.se.q1)

#Result Skewness of hat{u}
fit.final<-summary(rq(Q1hat.Small ~ KexpK1.Small, tau=taus[ptau.q1]),se="ker")
skewness(fit.final$resid)

#############
##Subsampling

#the number of Subsampling

num.ss <- 199
tau.ss.dist <- array(0,dim=c(num.ss,1))

taus<-seq(0.03,0.97,0.01)
nquantiles<-length(taus)

n.id<-max(data.Small$id)
n.ss<-as.integer(n.id^(0.8))

for(j in 1:num.ss){

sub.ss <- sample.int(n.id, n.ss, replace = F)
data.Small.ss <- data.Small %>%filter(data.Small$id %in% sub.ss)
data.Small.ss<- transform(data.Small.ss,id = as.numeric(factor(data.Small.ss$id)))

Zaux.Small.ss <- model.matrix(~as.factor(data.Small.ss$id)-1)
Z.Small.ss<-Zaux.Small.ss[,-1]
ncolZ.Small.ss<-ncol(Z.Small.ss)

KexpK1.Small.ss<-data.Small.ss$KexpK1
Q1.Small.ss<-data.Small.ss$Q1

fitols.Small.ss<-summary(lm(Q1.Small.ss ~ KexpK1.Small.ss+Z.Small.ss))
Q1hat.Small.ss<-Q1.Small.ss-(Z.Small.ss%*%fitols.Small.ss$coef[3:(ncolZ.Small.ss+2)])

X.ss<-as.matrix.csr(cbind(1,KexpK1.Small.ss))
Y.ss<-Q1hat.Small.ss

tcoefs.ss.result <- foreach(k = 1:nquantiles, .combine='cbind', .packages = c("quantreg", "MASS")) %dopar% {

reg.ss <- rq.fit.sfn(X.ss, Y.ss, tau=taus[k])

  t=taus[k]
  sign=(reg.ss$residuals<0)
  rho=(reg.ss$residuals)*(t-sign)
  score1.q1= (mean(reg.ss$residuals))/mean(rho) 
  score2.q1=(1-2*t)/(t*(1-t))
  tcoefs.aux.b<- abs(score2.q1-score1.q1)

return(tcoefs.aux.b)

}

ptau.q1=which.min(tcoefs.ss.result)
tau.ss.dist[j]=taus[ptau.q1]

}

sd(tau.ss.dist)


############
# Q1 on KexpK1 CONSTRAINED
# Constrained - Large==1

taus<-seq(0.03,0.97,0.001)
nquantiles<-length(taus)

score1.q1=matrix(0,nquantiles)
score2.q1=matrix(0,nquantiles)
mu.q1=matrix(0,nquantiles)
mu.se.q1=matrix(0,nquantiles)
alpha.q1=matrix(0,nquantiles)
alpha.se.q1=matrix(0,nquantiles)

data.Large<- data %>%filter(data$Large==1)
data.Large<- transform(data.Large,id = as.numeric(factor(data.Large$id)))

Zaux.Large <- model.matrix(~as.factor(data.Large$id)-1)
Z.Large<-Zaux.Large[,-1]
ncolZ.Large<-ncol(Z.Large)

KexpK1.Large<-data.Large$KexpK1
Q1.Large<-data.Large$Q1

fitols.Large<-summary(lm(Q1.Large ~ KexpK1.Large+Z.Large))
Q1hat.Large<-Q1.Large-(Z.Large%*%fitols.Large$coef[3:(ncolZ.Large+2)])

for(i in 1:nquantiles){
  
  fit<-summary(rq(Q1hat.Large ~ KexpK1.Large, tau=taus[i]),se="ker")
  
  mu.q1[i]=fit$coef[1]
  alpha.q1[i]=fit$coef[2]
  mu.se.q1[i]=fit$coef[3]
  alpha.se.q1[i]=fit$coef[4]
  t=taus[i] 
  sign=(fit$residuals<0)
  rho=(fit$residuals)*(t-sign)
  score1.q1[i]= (mean(fit$residuals))/mean(rho) 
  score2.q1[i]=(1-2*t)/(t*(1-t))
}

score.q1=abs(score2.q1-score1.q1)
ptau.q1=which.min(score.q1)

tau.hat.q1=taus[ptau.q1]
mu.hat.q1=mu.q1[ptau.q1]
mu.hat.se.q1=mu.se.q1[ptau.q1]
alpha.hat.q1=alpha.q1[ptau.q1]
alpha.hat.se.q1=alpha.se.q1[ptau.q1]

plot(abs(score2.q1-score1.q1))


#Results Q1 on KexpK1 Large==1
tau.hat.q1
cbind(alpha.hat.q1,alpha.hat.se.q1)

#Result Skewness of hat{u}
fit.final<-summary(rq(Q1hat.Large ~ KexpK1.Large, tau=taus[ptau.q1]),se="ker")
skewness(fit.final$resid)

#############
##Subsampling

#the number of Subsampling

num.ss <- 199
tau.ss.dist <- array(0,dim=c(num.ss,1))

taus<-seq(0.03,0.97,0.01)
nquantiles<-length(taus)

n.id<-max(data.Large$id)
n.ss<-as.integer(n.id^(0.8))

for(j in 1:num.ss){

sub.ss <- sample.int(n.id, n.ss, replace = F)
data.Large.ss <- data.Large %>%filter(data.Large$id %in% sub.ss)
data.Large.ss<- transform(data.Large.ss,id = as.numeric(factor(data.Large.ss$id)))

Zaux.Large.ss <- model.matrix(~as.factor(data.Large.ss$id)-1)
Z.Large.ss<-Zaux.Large.ss[,-1]
ncolZ.Large.ss<-ncol(Z.Large.ss)

KexpK1.Large.ss<-data.Large.ss$KexpK1
Q1.Large.ss<-data.Large.ss$Q1

fitols.Large.ss<-summary(lm(Q1.Large.ss ~ KexpK1.Large.ss+Z.Large.ss))
Q1hat.Large.ss<-Q1.Large.ss-(Z.Large.ss%*%fitols.Large.ss$coef[3:(ncolZ.Large.ss+2)])

X.ss<-as.matrix.csr(cbind(1,KexpK1.Large.ss))
Y.ss<-Q1hat.Large.ss

tcoefs.ss.result <- foreach(k = 1:nquantiles, .combine='cbind', .packages = c("quantreg", "MASS")) %dopar% {

reg.ss <- rq.fit.sfn(X.ss, Y.ss, tau=taus[k])

  t=taus[k]
  sign=(reg.ss$residuals<0)
  rho=(reg.ss$residuals)*(t-sign)
  score1.q1= (mean(reg.ss$residuals))/mean(rho) 
  score2.q1=(1-2*t)/(t*(1-t))
  tcoefs.aux.b<- abs(score2.q1-score1.q1)

return(tcoefs.aux.b)

}

ptau.q1=which.min(tcoefs.ss.result)
tau.ss.dist[j]=taus[ptau.q1]

}

sd(tau.ss.dist)



############
# Q1 on KexpK1 CONSTRAINED
# Constrained - LowPayout == 1


taus<-seq(0.03,0.97,0.001)
nquantiles<-length(taus)

score1.q1=matrix(0,nquantiles)
score2.q1=matrix(0,nquantiles)
mu.q1=matrix(0,nquantiles)
mu.se.q1=matrix(0,nquantiles)
alpha.q1=matrix(0,nquantiles)
alpha.se.q1=matrix(0,nquantiles)

data.LowPayout<- data %>%filter(data$LowPayout==1)
data.LowPayout<- transform(data.LowPayout,id = as.numeric(factor(data.LowPayout$id)))

Zaux.LowPayout <- model.matrix(~as.factor(data.LowPayout$id)-1)
Z.LowPayout<-Zaux.LowPayout[,-1]
ncolZ.LowPayout<-ncol(Z.LowPayout)

KexpK1.LowPayout<-data.LowPayout$KexpK1
Q1.LowPayout<-data.LowPayout$Q1

fitols.LowPayout<-summary(lm(Q1.LowPayout ~ KexpK1.LowPayout+Z.LowPayout))
Q1hat.LowPayout<-Q1.LowPayout-(Z.LowPayout%*%fitols.LowPayout$coef[3:(ncolZ.LowPayout+2)])


for(i in 1:nquantiles){
  
  fit<-summary(rq(Q1hat.LowPayout ~ KexpK1.LowPayout, tau=taus[i]),se="ker")
  
  mu.q1[i]=fit$coef[1]
  alpha.q1[i]=fit$coef[2]
  mu.se.q1[i]=fit$coef[3]
  alpha.se.q1[i]=fit$coef[4]
  t=taus[i] 
  sign=(fit$residuals<0)
  rho=(fit$residuals)*(t-sign)
  score1.q1[i]= (mean(fit$residuals))/mean(rho) 
  score2.q1[i]=(1-2*t)/(t*(1-t))
}

score.q1=abs(score2.q1-score1.q1)
ptau.q1=which.min(score.q1)

tau.hat.q1=taus[ptau.q1]
mu.hat.q1=mu.q1[ptau.q1]
mu.hat.se.q1=mu.se.q1[ptau.q1]
alpha.hat.q1=alpha.q1[ptau.q1]
alpha.hat.se.q1=alpha.se.q1[ptau.q1]

plot(abs(score2.q1-score1.q1))

#Results Q1 on KexpK1 LowPayout==1
tau.hat.q1
cbind(alpha.hat.q1,alpha.hat.se.q1)

#Result Skewness of hat{u}
fit.final<-summary(rq(Q1hat.LowPayout ~ KexpK1.LowPayout, tau=taus[ptau.q1]),se="ker")
skewness(fit.final$resid)

#############
##Subsampling

#the number of Subsampling

num.ss <- 199
tau.ss.dist <- array(0,dim=c(num.ss,1))

taus<-seq(0.03,0.97,0.01)
nquantiles<-length(taus)

n.id<-max(data.LowPayout$id)
n.ss<-as.integer(n.id^(0.8))

for(j in 1:num.ss){
  
  sub.ss <- sample.int(n.id, n.ss, replace = F)
  data.LowPayout.ss <- data.LowPayout %>%filter(data.LowPayout$id %in% sub.ss)
  data.LowPayout.ss<- transform(data.LowPayout.ss,id = as.numeric(factor(data.LowPayout.ss$id)))
  
  Zaux.LowPayout.ss <- model.matrix(~as.factor(data.LowPayout.ss$id)-1)
  Z.LowPayout.ss<-Zaux.LowPayout.ss[,-1]
  ncolZ.LowPayout.ss<-ncol(Z.LowPayout.ss)
  
  KexpK1.LowPayout.ss<-data.LowPayout.ss$KexpK1
  Q1.LowPayout.ss<-data.LowPayout.ss$Q1
  
  fitols.LowPayout.ss<-summary(lm(Q1.LowPayout.ss ~ KexpK1.LowPayout.ss+Z.LowPayout.ss))
  Q1hat.LowPayout.ss<-Q1.LowPayout.ss-(Z.LowPayout.ss%*%fitols.LowPayout.ss$coef[3:(ncolZ.LowPayout.ss+2)])
  
  X.ss<-as.matrix.csr(cbind(1,KexpK1.LowPayout.ss))
  Y.ss<-Q1hat.LowPayout.ss
  
  tcoefs.ss.result <- foreach(k = 1:nquantiles, .combine='cbind', .packages = c("quantreg", "MASS")) %dopar% {
    
    reg.ss <- rq.fit.sfn(X.ss, Y.ss, tau=taus[k])
    
    t=taus[k]
    sign=(reg.ss$residuals<0)
    rho=(reg.ss$residuals)*(t-sign)
    score1.q1= (mean(reg.ss$residuals))/mean(rho) 
    score2.q1=(1-2*t)/(t*(1-t))
    tcoefs.aux.b<- abs(score2.q1-score1.q1)
    
    return(tcoefs.aux.b)
    
  }
  
  ptau.q1=which.min(tcoefs.ss.result)
  tau.ss.dist[j]=taus[ptau.q1]
  
}

sd(tau.ss.dist)



############
# Q1 on KexpK1 CONSTRAINED
# Constrained - HighPayout==1

taus<-seq(0.03,0.97,0.001)
nquantiles<-length(taus)

score1.q1=matrix(0,nquantiles)
score2.q1=matrix(0,nquantiles)
mu.q1=matrix(0,nquantiles)
mu.se.q1=matrix(0,nquantiles)
alpha.q1=matrix(0,nquantiles)
alpha.se.q1=matrix(0,nquantiles)

data.HighPayout<- data %>%filter(data$HighPayout==1)
data.HighPayout<- transform(data.HighPayout,id = as.numeric(factor(data.HighPayout$id)))

Zaux.HighPayout <- model.matrix(~as.factor(data.HighPayout$id)-1)
Z.HighPayout<-Zaux.HighPayout[,-1]
ncolZ.HighPayout<-ncol(Z.HighPayout)

KexpK1.HighPayout<-data.HighPayout$KexpK1
Q1.HighPayout<-data.HighPayout$Q1

fitols.HighPayout<-summary(lm(Q1.HighPayout ~ KexpK1.HighPayout+Z.HighPayout))
Q1hat.HighPayout<-Q1.HighPayout-(Z.HighPayout%*%fitols.HighPayout$coef[3:(ncolZ.HighPayout+2)])


for(i in 1:nquantiles){
  
  fit<-summary(rq(Q1hat.HighPayout ~ KexpK1.HighPayout, tau=taus[i]),se="ker")
  
  mu.q1[i]=fit$coef[1]
  alpha.q1[i]=fit$coef[2]
  mu.se.q1[i]=fit$coef[3]
  alpha.se.q1[i]=fit$coef[4]
  t=taus[i] 
  sign=(fit$residuals<0)
  rho=(fit$residuals)*(t-sign)
  score1.q1[i]= (mean(fit$residuals))/mean(rho) 
  score2.q1[i]=(1-2*t)/(t*(1-t))
}

score.q1=abs(score2.q1-score1.q1)
ptau.q1=which.min(score.q1)

tau.hat.q1=taus[ptau.q1]
mu.hat.q1=mu.q1[ptau.q1]
mu.hat.se.q1=mu.se.q1[ptau.q1]
alpha.hat.q1=alpha.q1[ptau.q1]
alpha.hat.se.q1=alpha.se.q1[ptau.q1]

plot(abs(score2.q1-score1.q1))


#Results Q1 on KexpK1 HighPayout==1
tau.hat.q1
cbind(alpha.hat.q1,alpha.hat.se.q1)

#Result Skewness of hat{u}
fit.final<-summary(rq(Q1hat.HighPayout ~ KexpK1.HighPayout, tau=taus[ptau.q1]),se="ker")
skewness(fit.final$resid)

#############
##Subsampling

#the number of Subsampling

num.ss <- 199
tau.ss.dist <- array(0,dim=c(num.ss,1))

taus<-seq(0.03,0.97,0.01)
nquantiles<-length(taus)

n.id<-max(data.HighPayout$id)
n.ss<-as.integer(n.id^(0.8))

for(j in 1:num.ss){

sub.ss <- sample.int(n.id, n.ss, replace = F)
data.HighPayout.ss <- data.HighPayout %>%filter(data.HighPayout$id %in% sub.ss)
data.HighPayout.ss<- transform(data.HighPayout.ss,id = as.numeric(factor(data.HighPayout.ss$id)))

Zaux.HighPayout.ss <- model.matrix(~as.factor(data.HighPayout.ss$id)-1)
Z.HighPayout.ss<-Zaux.HighPayout.ss[,-1]
ncolZ.HighPayout.ss<-ncol(Z.HighPayout.ss)

KexpK1.HighPayout.ss<-data.HighPayout.ss$KexpK1
Q1.HighPayout.ss<-data.HighPayout.ss$Q1

fitols.HighPayout.ss<-summary(lm(Q1.HighPayout.ss ~ KexpK1.HighPayout.ss+Z.HighPayout.ss))
Q1hat.HighPayout.ss<-Q1.HighPayout.ss-(Z.HighPayout.ss%*%fitols.HighPayout.ss$coef[3:(ncolZ.HighPayout.ss+2)])

X.ss<-as.matrix.csr(cbind(1,KexpK1.HighPayout.ss))
Y.ss<-Q1hat.HighPayout.ss

tcoefs.ss.result <- foreach(k = 1:nquantiles, .combine='cbind', .packages = c("quantreg", "MASS")) %dopar% {

reg.ss <- rq.fit.sfn(X.ss, Y.ss, tau=taus[k])

  t=taus[k]
  sign=(reg.ss$residuals<0)
  rho=(reg.ss$residuals)*(t-sign)
  score1.q1= (mean(reg.ss$residuals))/mean(rho) 
  score2.q1=(1-2*t)/(t*(1-t))
  tcoefs.aux.b<- abs(score2.q1-score1.q1)

return(tcoefs.aux.b)

}

ptau.q1=which.min(tcoefs.ss.result)
tau.ss.dist[j]=taus[ptau.q1]

}

sd(tau.ss.dist)



############
# Q1 on KexpK1 CONSTRAINED
# Constrained - LowDebt==1


taus<-seq(0.03,0.97,0.001)
nquantiles<-length(taus)

score1.q1=matrix(0,nquantiles)
score2.q1=matrix(0,nquantiles)
mu.q1=matrix(0,nquantiles)
mu.se.q1=matrix(0,nquantiles)
alpha.q1=matrix(0,nquantiles)
alpha.se.q1=matrix(0,nquantiles)

data.LowDebt<- data %>%filter(data$LowDebt==1)
data.LowDebt<- transform(data.LowDebt,id = as.numeric(factor(data.LowDebt$id)))

Zaux.LowDebt <- model.matrix(~as.factor(data.LowDebt$id)-1)
Z.LowDebt<-Zaux.LowDebt[,-1]
ncolZ.LowDebt<-ncol(Z.LowDebt)

KexpK1.LowDebt<-data.LowDebt$KexpK1
Q1.LowDebt<-data.LowDebt$Q1

fitols.LowDebt<-summary(lm(Q1.LowDebt ~ KexpK1.LowDebt+Z.LowDebt))
Q1hat.LowDebt<-Q1.LowDebt-(Z.LowDebt%*%fitols.LowDebt$coef[3:(ncolZ.LowDebt+2)])


for(i in 1:nquantiles){
  
  fit<-summary(rq(Q1hat.LowDebt ~ KexpK1.LowDebt, tau=taus[i]),se="ker")
  
  mu.q1[i]=fit$coef[1]
  alpha.q1[i]=fit$coef[2]
  mu.se.q1[i]=fit$coef[3]
  alpha.se.q1[i]=fit$coef[4]
  t=taus[i] 
  sign=(fit$residuals<0)
  rho=(fit$residuals)*(t-sign)
  score1.q1[i]= (mean(fit$residuals))/mean(rho) 
  score2.q1[i]=(1-2*t)/(t*(1-t))
}

score.q1=abs(score2.q1-score1.q1)
ptau.q1=which.min(score.q1)

tau.hat.q1=taus[ptau.q1]
mu.hat.q1=mu.q1[ptau.q1]
mu.hat.se.q1=mu.se.q1[ptau.q1]
alpha.hat.q1=alpha.q1[ptau.q1]
alpha.hat.se.q1=alpha.se.q1[ptau.q1]

plot(abs(score2.q1-score1.q1))

#Results Q1 on KexpK1 LowDebt==1
tau.hat.q1
cbind(alpha.hat.q1,alpha.hat.se.q1)

#Result Skewness of hat{u}
fit.final<-summary(rq(Q1hat.LowDebt ~ KexpK1.LowDebt, tau=taus[ptau.q1]),se="ker")
skewness(fit.final$resid)

#############
##Subsampling

#the number of Subsampling

num.ss <- 199
tau.ss.dist <- array(0,dim=c(num.ss,1))

taus<-seq(0.03,0.97,0.01)
nquantiles<-length(taus)

n.id<-max(data.LowDebt$id)
n.ss<-as.integer(n.id^(0.8))

for(j in 1:num.ss){

sub.ss <- sample.int(n.id, n.ss, replace = F)
data.LowDebt.ss <- data.LowDebt %>%filter(data.LowDebt$id %in% sub.ss)
data.LowDebt.ss<- transform(data.LowDebt.ss,id = as.numeric(factor(data.LowDebt.ss$id)))

Zaux.LowDebt.ss <- model.matrix(~as.factor(data.LowDebt.ss$id)-1)
Z.LowDebt.ss<-Zaux.LowDebt.ss[,-1]
ncolZ.LowDebt.ss<-ncol(Z.LowDebt.ss)

KexpK1.LowDebt.ss<-data.LowDebt.ss$KexpK1
Q1.LowDebt.ss<-data.LowDebt.ss$Q1

fitols.LowDebt.ss<-summary(lm(Q1.LowDebt.ss ~ KexpK1.LowDebt.ss+Z.LowDebt.ss))
Q1hat.LowDebt.ss<-Q1.LowDebt.ss-(Z.LowDebt.ss%*%fitols.LowDebt.ss$coef[3:(ncolZ.LowDebt.ss+2)])

X.ss<-as.matrix.csr(cbind(1,KexpK1.LowDebt.ss))
Y.ss<-Q1hat.LowDebt.ss

tcoefs.ss.result <- foreach(k = 1:nquantiles, .combine='cbind', .packages = c("quantreg", "MASS")) %dopar% {

reg.ss <- rq.fit.sfn(X.ss, Y.ss, tau=taus[k])

  t=taus[k]
  sign=(reg.ss$residuals<0)
  rho=(reg.ss$residuals)*(t-sign)
  score1.q1= (mean(reg.ss$residuals))/mean(rho) 
  score2.q1=(1-2*t)/(t*(1-t))
  tcoefs.aux.b<- abs(score2.q1-score1.q1)

return(tcoefs.aux.b)

}

ptau.q1=which.min(tcoefs.ss.result)
tau.ss.dist[j]=taus[ptau.q1]

}

sd(tau.ss.dist)



############
# Q1 on KexpK1 CONSTRAINED
# Constrained - HighDebt==1

taus<-seq(0.03,0.97,0.001)
nquantiles<-length(taus)

score1.q1=matrix(0,nquantiles)
score2.q1=matrix(0,nquantiles)
mu.q1=matrix(0,nquantiles)
mu.se.q1=matrix(0,nquantiles)
alpha.q1=matrix(0,nquantiles)
alpha.se.q1=matrix(0,nquantiles)

data.HighDebt<- data %>%filter(data$HighDebt==1)
data.HighDebt<- transform(data.HighDebt,id = as.numeric(factor(data.HighDebt$id)))

Zaux.HighDebt <- model.matrix(~as.factor(data.HighDebt$id)-1)
Z.HighDebt<-Zaux.HighDebt[,-1]
ncolZ.HighDebt<-ncol(Z.HighDebt)

KexpK1.HighDebt<-data.HighDebt$KexpK1
Q1.HighDebt<-data.HighDebt$Q1

fitols.HighDebt<-summary(lm(Q1.HighDebt ~ KexpK1.HighDebt+Z.HighDebt))
Q1hat.HighDebt<-Q1.HighDebt-(Z.HighDebt%*%fitols.HighDebt$coef[3:(ncolZ.HighDebt+2)])


for(i in 1:nquantiles){
  
  fit<-summary(rq(Q1hat.HighDebt ~ KexpK1.HighDebt, tau=taus[i]),se="ker")
  
  mu.q1[i]=fit$coef[1]
  alpha.q1[i]=fit$coef[2]
  mu.se.q1[i]=fit$coef[3]
  alpha.se.q1[i]=fit$coef[4]
  t=taus[i] 
  sign=(fit$residuals<0)
  rho=(fit$residuals)*(t-sign)
  score1.q1[i]= (mean(fit$residuals))/mean(rho) 
  score2.q1[i]=(1-2*t)/(t*(1-t))
}

score.q1=abs(score2.q1-score1.q1)
ptau.q1=which.min(score.q1)

tau.hat.q1=taus[ptau.q1]
mu.hat.q1=mu.q1[ptau.q1]
mu.hat.se.q1=mu.se.q1[ptau.q1]
alpha.hat.q1=alpha.q1[ptau.q1]
alpha.hat.se.q1=alpha.se.q1[ptau.q1]

plot(abs(score2.q1-score1.q1))


#Results Q1 on KexpK1 HighDebt==1
tau.hat.q1
cbind(alpha.hat.q1,alpha.hat.se.q1)

#Result Skewness of hat{u}
fit.final<-summary(rq(Q1hat.HighDebt ~ KexpK1.HighDebt, tau=taus[ptau.q1]),se="ker")
skewness(fit.final$resid)

#############
##Subsampling

#the number of Subsampling

num.ss <- 199
tau.ss.dist <- array(0,dim=c(num.ss,1))

taus<-seq(0.03,0.97,0.01)
nquantiles<-length(taus)

n.id<-max(data.HighDebt$id)
n.ss<-as.integer(n.id^(0.8))

for(j in 1:num.ss){

sub.ss <- sample.int(n.id, n.ss, replace = F)
data.HighDebt.ss <- data.HighDebt %>%filter(data.HighDebt$id %in% sub.ss)
data.HighDebt.ss<- transform(data.HighDebt.ss,id = as.numeric(factor(data.HighDebt.ss$id)))

Zaux.HighDebt.ss <- model.matrix(~as.factor(data.HighDebt.ss$id)-1)
Z.HighDebt.ss<-Zaux.HighDebt.ss[,-1]
ncolZ.HighDebt.ss<-ncol(Z.HighDebt.ss)

KexpK1.HighDebt.ss<-data.HighDebt.ss$KexpK1
Q1.HighDebt.ss<-data.HighDebt.ss$Q1

fitols.HighDebt.ss<-summary(lm(Q1.HighDebt.ss ~ KexpK1.HighDebt.ss+Z.HighDebt.ss))
Q1hat.HighDebt.ss<-Q1.HighDebt.ss-(Z.HighDebt.ss%*%fitols.HighDebt.ss$coef[3:(ncolZ.HighDebt.ss+2)])

X.ss<-as.matrix.csr(cbind(1,KexpK1.HighDebt.ss))
Y.ss<-Q1hat.HighDebt.ss

tcoefs.ss.result <- foreach(k = 1:nquantiles, .combine='cbind', .packages = c("quantreg", "MASS")) %dopar% {

reg.ss <- rq.fit.sfn(X.ss, Y.ss, tau=taus[k])

  t=taus[k]
  sign=(reg.ss$residuals<0)
  rho=(reg.ss$residuals)*(t-sign)
  score1.q1= (mean(reg.ss$residuals))/mean(rho) 
  score2.q1=(1-2*t)/(t*(1-t))
  tcoefs.aux.b<- abs(score2.q1-score1.q1)

return(tcoefs.aux.b)

}

ptau.q1=which.min(tcoefs.ss.result)
tau.ss.dist[j]=taus[ptau.q1]

}

sd(tau.ss.dist)



############
# Q1 on KexpK1 CONSTRAINED
# Constrained - LowCash==1

taus<-seq(0.03,0.97,0.001)
nquantiles<-length(taus)

score1.q1=matrix(0,nquantiles)
score2.q1=matrix(0,nquantiles)
mu.q1=matrix(0,nquantiles)
mu.se.q1=matrix(0,nquantiles)
alpha.q1=matrix(0,nquantiles)
alpha.se.q1=matrix(0,nquantiles)

data.LowCash<- data %>%filter(data$LowCash==1)
data.LowCash<- transform(data.LowCash,id = as.numeric(factor(data.LowCash$id)))

Zaux.LowCash <- model.matrix(~as.factor(data.LowCash$id)-1)
Z.LowCash<-Zaux.LowCash[,-1]
ncolZ.LowCash<-ncol(Z.LowCash)

KexpK1.LowCash<-data.LowCash$KexpK1
Q1.LowCash<-data.LowCash$Q1

fitols.LowCash<-summary(lm(Q1.LowCash ~ KexpK1.LowCash+Z.LowCash))
Q1hat.LowCash<-Q1.LowCash-(Z.LowCash%*%fitols.LowCash$coef[3:(ncolZ.LowCash+2)])


for(i in 1:nquantiles){
  
  fit<-summary(rq(Q1hat.LowCash ~ KexpK1.LowCash, tau=taus[i]),se="ker")
  
  mu.q1[i]=fit$coef[1]
  alpha.q1[i]=fit$coef[2]
  mu.se.q1[i]=fit$coef[3]
  alpha.se.q1[i]=fit$coef[4]
  t=taus[i] 
  sign=(fit$residuals<0)
  rho=(fit$residuals)*(t-sign)
  score1.q1[i]= (mean(fit$residuals))/mean(rho) 
  score2.q1[i]=(1-2*t)/(t*(1-t))
}

score.q1=abs(score2.q1-score1.q1)
ptau.q1=which.min(score.q1)

tau.hat.q1=taus[ptau.q1]
mu.hat.q1=mu.q1[ptau.q1]
mu.hat.se.q1=mu.se.q1[ptau.q1]
alpha.hat.q1=alpha.q1[ptau.q1]
alpha.hat.se.q1=alpha.se.q1[ptau.q1]

plot(abs(score2.q1-score1.q1))

#Results Q1 on KexpK1 LowCash==1
tau.hat.q1
cbind(alpha.hat.q1,alpha.hat.se.q1)

#Result Skewness of hat{u}
fit.final<-summary(rq(Q1hat.LowCash ~ KexpK1.LowCash, tau=taus[ptau.q1]),se="ker")
skewness(fit.final$resid)

#############
##Subsampling

#the number of Subsampling

num.ss <- 199
tau.ss.dist <- array(0,dim=c(num.ss,1))

taus<-seq(0.03,0.97,0.01)
nquantiles<-length(taus)

n.id<-max(data.LowCash$id)
n.ss<-as.integer(n.id^(0.8))

for(j in 1:num.ss){

sub.ss <- sample.int(n.id, n.ss, replace = F)
data.LowCash.ss <- data.LowCash %>%filter(data.LowCash$id %in% sub.ss)
data.LowCash.ss<- transform(data.LowCash.ss,id = as.numeric(factor(data.LowCash.ss$id)))

Zaux.LowCash.ss <- model.matrix(~as.factor(data.LowCash.ss$id)-1)
Z.LowCash.ss<-Zaux.LowCash.ss[,-1]
ncolZ.LowCash.ss<-ncol(Z.LowCash.ss)

KexpK1.LowCash.ss<-data.LowCash.ss$KexpK1
Q1.LowCash.ss<-data.LowCash.ss$Q1

fitols.LowCash.ss<-summary(lm(Q1.LowCash.ss ~ KexpK1.LowCash.ss+Z.LowCash.ss))
Q1hat.LowCash.ss<-Q1.LowCash.ss-(Z.LowCash.ss%*%fitols.LowCash.ss$coef[3:(ncolZ.LowCash.ss+2)])

X.ss<-as.matrix.csr(cbind(1,KexpK1.LowCash.ss))
Y.ss<-Q1hat.LowCash.ss

tcoefs.ss.result <- foreach(k = 1:nquantiles, .combine='cbind', .packages = c("quantreg", "MASS")) %dopar% {

reg.ss <- rq.fit.sfn(X.ss, Y.ss, tau=taus[k])

  t=taus[k]
  sign=(reg.ss$residuals<0)
  rho=(reg.ss$residuals)*(t-sign)
  score1.q1= (mean(reg.ss$residuals))/mean(rho) 
  score2.q1=(1-2*t)/(t*(1-t))
  tcoefs.aux.b<- abs(score2.q1-score1.q1)

return(tcoefs.aux.b)

}

ptau.q1=which.min(tcoefs.ss.result)
tau.ss.dist[j]=taus[ptau.q1]

}

sd(tau.ss.dist)



############
# Q1 on KexpK1 CONSTRAINED
# Constrained - HighCash==1

taus<-seq(0.03,0.97,0.001)
nquantiles<-length(taus)

score1.q1=matrix(0,nquantiles)
score2.q1=matrix(0,nquantiles)
mu.q1=matrix(0,nquantiles)
mu.se.q1=matrix(0,nquantiles)
alpha.q1=matrix(0,nquantiles)
alpha.se.q1=matrix(0,nquantiles)

data.HighCash<- data %>%filter(data$HighCash==1)
data.HighCash<- transform(data.HighCash,id = as.numeric(factor(data.HighCash$id)))

Zaux.HighCash <- model.matrix(~as.factor(data.HighCash$id)-1)
Z.HighCash<-Zaux.HighCash[,-1]
ncolZ.HighCash<-ncol(Z.HighCash)

KexpK1.HighCash<-data.HighCash$KexpK1
Q1.HighCash<-data.HighCash$Q1

fitols.HighCash<-summary(lm(Q1.HighCash ~ KexpK1.HighCash+Z.HighCash))
Q1hat.HighCash<-Q1.HighCash-(Z.HighCash%*%fitols.HighCash$coef[3:(ncolZ.HighCash+2)])


for(i in 1:nquantiles){
  
  fit<-summary(rq(Q1hat.HighCash ~ KexpK1.HighCash, tau=taus[i]),se="ker")
  
  mu.q1[i]=fit$coef[1]
  alpha.q1[i]=fit$coef[2]
  mu.se.q1[i]=fit$coef[3]
  alpha.se.q1[i]=fit$coef[4]
  t=taus[i] 
  sign=(fit$residuals<0)
  rho=(fit$residuals)*(t-sign)
  score1.q1[i]= (mean(fit$residuals))/mean(rho) 
  score2.q1[i]=(1-2*t)/(t*(1-t))
}

score.q1=abs(score2.q1-score1.q1)
ptau.q1=which.min(score.q1)

tau.hat.q1=taus[ptau.q1]
mu.hat.q1=mu.q1[ptau.q1]
mu.hat.se.q1=mu.se.q1[ptau.q1]
alpha.hat.q1=alpha.q1[ptau.q1]
alpha.hat.se.q1=alpha.se.q1[ptau.q1]

plot(abs(score2.q1-score1.q1))

#Results Q1 on KexpK1 HighCash==1
tau.hat.q1
cbind(alpha.hat.q1,alpha.hat.se.q1)

#Result Skewness of hat{u}
fit.final<-summary(rq(Q1hat.HighCash ~ KexpK1.HighCash, tau=taus[ptau.q1]),se="ker")
skewness(fit.final$resid)

#############
##Subsampling

#the number of Subsampling

num.ss <- 199
tau.ss.dist <- array(0,dim=c(num.ss,1))

taus<-seq(0.03,0.97,0.01)
nquantiles<-length(taus)

n.id<-max(data.HighCash$id)
n.ss<-as.integer(n.id^(0.8))

for(j in 1:num.ss){

sub.ss <- sample.int(n.id, n.ss, replace = F)
data.HighCash.ss <- data.HighCash %>%filter(data.HighCash$id %in% sub.ss)
data.HighCash.ss<- transform(data.HighCash.ss,id = as.numeric(factor(data.HighCash.ss$id)))

Zaux.HighCash.ss <- model.matrix(~as.factor(data.HighCash.ss$id)-1)
Z.HighCash.ss<-Zaux.HighCash.ss[,-1]
ncolZ.HighCash.ss<-ncol(Z.HighCash.ss)

KexpK1.HighCash.ss<-data.HighCash.ss$KexpK1
Q1.HighCash.ss<-data.HighCash.ss$Q1

fitols.HighCash.ss<-summary(lm(Q1.HighCash.ss ~ KexpK1.HighCash.ss+Z.HighCash.ss))
Q1hat.HighCash.ss<-Q1.HighCash.ss-(Z.HighCash.ss%*%fitols.HighCash.ss$coef[3:(ncolZ.HighCash.ss+2)])

X.ss<-as.matrix.csr(cbind(1,KexpK1.HighCash.ss))
Y.ss<-Q1hat.HighCash.ss

tcoefs.ss.result <- foreach(k = 1:nquantiles, .combine='cbind', .packages = c("quantreg", "MASS")) %dopar% {

reg.ss <- rq.fit.sfn(X.ss, Y.ss, tau=taus[k])

  t=taus[k]
  sign=(reg.ss$residuals<0)
  rho=(reg.ss$residuals)*(t-sign)
  score1.q1= (mean(reg.ss$residuals))/mean(rho) 
  score2.q1=(1-2*t)/(t*(1-t))
  tcoefs.aux.b<- abs(score2.q1-score1.q1)

return(tcoefs.aux.b)

}

ptau.q1=which.min(tcoefs.ss.result)
tau.ss.dist[j]=taus[ptau.q1]

}

sd(tau.ss.dist)


