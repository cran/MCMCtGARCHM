#######################################################
# 28/10/2009 #
# Estimate the GARCH(1,1)-in-mean model with t errors #
# using the random-walk Metropolis-Hastings algorithm #
#######################################################
rm(list=ls(all=TRUE))
##mcmctgarchm

mctgarchm=function(y, m)
{

# log posterior of parameters
garch11.post<-function(xp)
{
a0=exp(xp[1])
a1=exp(xp[2])/(1+exp(xp[2]))
b1=exp(xp[3])/(1+exp(xp[3]))
delta=xp[4]
nv=exp(xp[5])

h=vector(mode="numeric",length=n)
h[1]=a0+a1*y0*y0+b1*h0
for(i in 2:n)
{
h[i]<-a0+a1*y.dat[i-1]*y.dat[i-1]+b1*h[i-1]
}
#log-likelihood of tgarch
tem2<-dt((y.dat-delta*sqrt(h))/sqrt(h),df=nv)/sqrt(h)
logpost<-sum(log(tem2))

#Jacobi
logpost=logpost+log(a1-a1*a1)+log(b1-b1*b1)+xp[1]+xp[4]
#prior of a0 is log-normal(0,0.5)
std=0.5
tem2=-0.5*log(a0)^2/std^2-log(a0*std*sqrt(2*pi))
logpost=logpost+tem2
#prior of nv
#tem2<-dnorm(nv,mean=8,sd=9)
logpost<-logpost-log(1+nv*nv)
#prior of delta
logpost<-logpost+dnorm(delta,mean=0,sd=3,log=TRUE)

return(logpost)
}
# Randon-walk Metropolis-Hastings Algorithm
rw.mh<-function(xp,lnpost)
{
fa<-lnpost
rn=rnorm(4,mean=0,sd=1)
rn=rn/sqrt(sum(rn*rn))
tem2=xp[1:4]+tune*rn
tem2=c(tem2,xp[5])
fb=garch11.post(tem2)
r=fb-fa
accept=0
if(r>0) accept=1
else
{ un=runif(1,min=0,max=1)
if(un<exp(r)) accept=1
}
if (accept==1)
{ xp=tem2
fa=fb
}
return(list(parameter=xp,accept=accept,logpost=fa))
}

rw.mh.nv<-function(xp,lnpost)
{
fa<-lnpost
rn=rnorm(1,mean=0,sd=1)
tem2=xp[5]
xp[5]=xp[5]+tune.nv*rnorm(1,mean=0,sd=1)
fb=garch11.post(xp)
r=fb-fa
accept=0
if(r>0) accept=1
else
{ un=runif(1,min=0,max=1)
if(un<exp(r)) accept=1
}
if(accept==1) fa=fb
else xp[5]=tem2

return(list(parameter=xp,accept=accept,logpost=fa))
}

#Main Program
set.seed(500)
n=length(y)-1
y.dat<-y[2:(n+1)]
y.dat<-y.dat-mean(y.dat)
y0<-y[1]
xp=numeric(4)
xp[1]<-log(0.05) #omega=0.05
xp[2]<-log(0.05/(1-0.05)) #alpha=0.05
xp[3]<-log(0.85/(1-0.85)) #beta =0.85
xp[4]<-0 #coefficient for garch-in-mean:gamma
xp[5]<-log(10) #nv=10
h0<-mean(y.dat*y.dat)
lnpost<-garch11.post(xp)
warm<-m/10
tune<-0.12
tune.nv<-0.8
para.matrix<-matrix(0,nr=m,nc=5)
accept.rate<-accept.rate.nv<-0
for(ks in 1:warm)
{
tem<-rw.mh(xp,lnpost)
xp<-tem$parameter
accept.rate<-accept.rate+tem[[2]]
lnpost<-tem[[3]]
tem<-rw.mh.nv(xp,lnpost)
xp<-tem$parameter
accept.rate.nv<-accept.rate.nv+tem[[2]]
lnpost<-tem[[3]]
}

for(ks in 1:m)
{
tem<-rw.mh(xp,lnpost)
xp<-tem$parameter
accept.rate<-accept.rate+tem[[2]]
lnpost<-tem[[3]]
tem<-rw.mh.nv(xp,lnpost)
xp<-tem$parameter
accept.rate.nv<-accept.rate.nv+tem[[2]]
lnpost<-tem[[3]]
para.matrix[ks,1]<-exp(xp[1])
para.matrix[ks,2]<-exp(xp[2])/(1+exp(xp[2]))
para.matrix[ks,3]<-exp(xp[3])/(1+exp(xp[3]))
para.matrix[ks,4]<-xp[4]
para.matrix[ks,5]<-exp(xp[5])
}


#Compute SIF
dm<-length(xp)
num.batch<-50
siz.batch<-m/num.batch
hat.para<-vector(mode="numeric",length=dm)
for(i in 1:dm)
{
hat.para[i]=mean(para.matrix[,i])
}
sif<-vector(mode="numeric",length=dm)
for(i in 1:dm)
{
tem<-matrix(para.matrix[,i],nc=num.batch,byrow=F)
batch.mn<-vector(mode="numeric",length=num.batch)
for(j in 1:num.batch)
{
batch.mn[j]<-mean(tem[,j])
}
var.batch<-siz.batch*sum((batch.mn-hat.para[i])^2)/(num.batch-1)
var.total<-sum((tem-hat.para[i])^2)/(m-1)
sif[i]<-var.batch/var.total
sif[i]<-round(sif[i],digit=2)
}
#Simulation Inefficient Factors
sif
round(hat.para,digit=4)
#Acceptance rate
accept.rate/(m+warm)
accept.rate.nv/(m+warm)


cat("omega:",hat.para[1],"sd of omega:",sd(para.matrix[,1]),"sif of omega:",sif[1],"\n") 
cat("alpha:",hat.para[2],"sd of alpha:",sd(para.matrix[,2]),"sif of alpha:",sif[2],"\n")
cat("beta:",hat.para[3],"sd of beta:",sd(para.matrix[,3]),"sif of beta",sif[3],"\n")
cat("gamma:",hat.para[4],"sd of gamma:",sd(para.matrix[,4]),"sif of beta",sif[4],"\n")
cat("nv",hat.para[5],"sd of nv:",sd(para.matrix[,5]),"sif of nv",sif[5],"\n")

return(para.matrix)
}

###
plotpara<-function(para.matrix, m)
{

#Keep one draw for every 10 draws
ik=1:(m/10)
ik=ik*10
tem=para.matrix[ik,]
par(mfrow=c(5,1),mar=c(4,4,4,4))
plot(tem[,1],typ='l',lty=1,col=4,xlab="Iteration",ylab="para",main="a0")
plot(tem[,2],typ='l',lty=1,col=4,xlab="Iteration",ylab="para",main="a1")
plot(tem[,3],typ='l',lty=1,col=4,xlab="Iteration",ylab="para",main="b1")

plot(tem[,4],typ='l',lty=1,col=4,xlab="Iteration",ylab="para",main="delta")
plot(tem[,5],typ='l',lty=1,col=4,xlab="Iteration",ylab="para",main="nv")
}

###Example:
#xt<-matrix(scan(file="\\MCMCTGARCHM\DATA\AUS2005.TXT"),ncol=1,byrow=T)
#m<-10000
#mc=mctgarchm(xt, m)
#plotpara(mc, m)



