##Import Rent Data
data<-read.csv("/Users/dshartra/Documents/GT-Spring-2016/Adv-Stat-Models/7401Project/RentData.csv")
attach(data)
library(car)
scatterplotMatrix(data)
##Plots of Potential Colinearity
plot(FamDep,Age)
plot(Supp,IPE)
plot(Supp,FamDep)
plot(IPE,FamDep)
plot(Supp,MI)
rentreg<-lm(Rent~.,data=data)
summary(rentreg)
X<-model.matrix(rentreg)[,-1]
round(cor(X),2)
v=numeric(6)
for(i in 1:6)
  v[i]=1/(1-summary(lm(X[,i] ~ X[,-i]))$r.squared)
v
##Residual Plots, Clear Heteroscedasticity, Non-Normality on Right Side
plot(rentreg)
##Variable Selection##
library(quadprog)
require(leaps)
sub<-regsubsets(Rent~.,data=data)
subsum<-summary(sub)
plot(1:6,subsum$cp,xlab="No. of Parameters",ylab="Cp Statistic")
abline(0,1)
##CP says use 5 predictors, excluding age
rentreg.step=step(rentreg)
summary(rentreg.step)
##AIC says use only IPE, FamDep, and Roommates
rentdata=data.frame(data)
rent.sca=as.data.frame(scale(as.matrix(rentdata)))
a=lm(Rent~.-1,data=rent.sca)
y = rent.sca$Rent
summary(a)
B=diag(a$coef)
Z=model.matrix(a)%*%B
D=t(Z)%*%Z
d=t(Z)%*%y
A=cbind(-1,diag(6))
M=seq(0.01,6,length=100) ##number of rows
gcv=numeric(100)
for(i in 1:100){
  b0=c(-M[i],rep(0,6))
  coef.nng=solve.QP(D,d,A,b0)$sol
  e=y-Z%*%coef.nng
  gcv[i]=sum(e^2)/(59*(1-M[i]/59)^2)
}
plot(M,gcv,type="l",main="Selection from Non-Negative Garrote",xlab="Number of Predictors",ylab="Generalized Cross-Validation")
M=M[which.min(gcv)]
b0=c(-M,rep(0,6))
coef.nng=round(solve.QP(D,d,A,b0)$sol,8)
coef.nng
beta.nng=B%*%coef.nng
round(cbind(a$coef,beta.nng),5)
e=y-Z%*%coef.nng
1-sum(e^2)/sum(y^2)
##NNG suggests the 3 variables as well
##Below is Strong Heredity Function
gstrong = function(p)
{
  mat = mat.or.vec(p,p*(p-1))
  c = 1
  m = p
  while(m > 1)
  {
    for (j in 1:(m-1))
    {
      mat[,c]= c(numeric(p - m),1,numeric(m-1))
      mat[,(c+p*(p-1)/2)] = c(numeric(p - m+j),1,numeric(m-(j+1)))
      c = c+1
    }
    m = m-1
  }
  cbind(-1,diag(p+p*(p-1)/2),rbind(mat,-cbind(diag(p*(p-1)/2),diag(p*(p-1)/2))))
}
A1=gstrong(3)
a1=lm(Rent~(IPE+FamDep+Roommates)^2.-1,data=rent.sca)
y = rent.sca$Rent
summary(a1)
B=diag(a1$coef)
Z=model.matrix(a1)%*%B
D=t(Z)%*%Z
d=t(Z)%*%y
M=seq(0.01,3,length=100) ##number of rows
gcv=numeric(100)
for(i in 1:100){
  b0=c(-M[i],rep(0,12))
  coef.nng=solve.QP(D,d,A1,b0)$sol
  e=y-Z%*%coef.nng
  gcv[i]=sum(e^2)/(59*(1-M[i]/59)^2)
}
plot(M,gcv,type="l")
M=M[which.min(gcv)]
b0=c(-M,rep(0,12))
coef.nng=round(solve.QP(D,d,A1,b0)$sol,8)
coef.nng
beta.nng=B%*%coef.nng
round(cbind(a1$coef,beta.nng),5)
e=y-Z%*%coef.nng
1-sum(e^2)/sum(y^2)
##No signficant interactions
##True residual analysis
finalreg <-lm(Rent~IPE+FamDep+Roommates,data=data1)
summary(finalreg)
plot(finalreg)
Xfinal=model.matrix(finalreg)
H=Xfinal%*%solve(t(Xfinal)%*%Xfinal)%*%t(Xfinal)
std.res=finalreg$res/(summary(a)$sig*sqrt(1-diag(H)))
qqnorm(std.res)
qqline(std.res)
##Still have bad residuals
##Need to transform
require(MASS)
library(faraway)
boxcox(finalreg, lambda=seq(-1,1,by=0.05),plotit=T) ##shows CI, shows evidence that log transformation can be used
X = model.matrix(finalreg)
ha = hat(X)
finalreg$res/(summary(finalreg)$sig*sqrt(1-ha))
r=rstandard(finalreg)
halfnorm(r)
data[11,] ##outlier 1
data[14,] ##someone put down that they had 14 roommates, potential data error, will remove
data[47,] ##outlier 2
data1<-(data[-c(14),])

reg<-lm(Rent~IPE+FamDep+Roommates,data=data1)
plot(reg)
boxcox(reg, lambda=seq(-1,1,by=0.05),plotit=T)
##USE THIS BOXCOX
logreg <-lm(log(Rent)~IPE+FamDep+Roommates,data=data1)
plot(logreg)
res = resid(logreg)
plot(data1$IPE,res)
plot(data1$FamDep,res)
plot(data1$Roommates,res)
plot((data1$FamDep),res)
##outliers possible because Im not considering location of housing
summary(logreg)


##Rerun interaction NNG after data point 14 has been removed.
rentdata1=data.frame(data1)
rent.sca1=as.data.frame(scale(as.matrix(rentdata1)))
A1=gstrong(3)
a1=lm(Rent~(IPE+FamDep+Roommates)^2.-1,data=rent.sca1)
y = rent.sca1$Rent
summary(a1)
B=diag(a1$coef)
Z=model.matrix(a1)%*%B
D=t(Z)%*%Z
d=t(Z)%*%y
M=seq(0.01,3,length=100) ##number of rows
gcv=numeric(100)
for(i in 1:100){
  b0=c(-M[i],rep(0,12))
  coef.nng=solve.QP(D,d,A1,b0)$sol
  e=y-Z%*%coef.nng
  gcv[i]=sum(e^2)/(59*(1-M[i]/59)^2)
}
plot(M,gcv,type="l")
M=M[which.min(gcv)]
b0=c(-M,rep(0,12))
coef.nng=round(solve.QP(D,d,A1,b0)$sol,8)
coef.nng
beta.nng=B%*%coef.nng
round(cbind(a1$coef,beta.nng),5)
e=y-Z%*%coef.nng
1-sum(e^2)/sum(y^2)
count = 0
for (i in 1:58){
  if (data1[i,5] > 0)
  count = count + 1
}
##Lasso finds IPE, FamDep and Roommates as variables
require(lars)
lmod <- lars(as.matrix(data[,-1]),data$Rent)
plot(lmod)
set.seed(123)
cvlmod <- cv.lars(as.matrix(data[,-1]),data$Rent)
cvlmod$index[which.min(cvlmod$cv)]
predict(lmod,s=0.5151515,type="coef",mode="fraction")$coef




