library(MASS)

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

# size of training data
n <- 10000
# number of variables
p <- 30

# quantile of the distrubution of SDs of columns of U
# to determine affected components
size=.999

# number of outliers
nout=1000
# normal data that will be contaminated with outliers
nin <- 4900

# how many time we want to replicate the experiment
nrep <- 100

# generate data
# Sigma is AR1 with coefficient 0.9
SS <- 0.9 ^ outer(1:p, 1:p, function(aa, bb) abs(aa - bb))
X <- mvrnorm(n=n, mu= rep(0,p), Sigma=SS)

# SVD of normal data
ss=svd(X)
U=ss$u
D=diag(ss$d)
V=ss$v

############################################
# bootstrap distribution
BT2 <- replicate(10000,{
  XB <-  X[sample(1:n,nin+nout,replace=TRUE),]
  UPB <- XB%*%V%*%solve(D)
  apply(UPB*sqrt(n-1),2,sd)
})

#boxplot(BT2~row(BT2),outline=F)
#abline(h=quantile(as.vector(BT2),size))
apply(BT2,1,quantile,size)

#####################################
st <- 1
en <- sum(eigen(cov(scale(X)))$values>1)
#############################################
sumo=data.frame(matrix(0,nin+nout,12))
# components along which we want to make shift
q <- sample(1:p,3,replace=F)

for (i in 1:nrep) 
{
  
  Z <- X%*%V
  index_att <- sample(1:n, nout)
  TL <- numeric(n)
  TL[index_att] <- 1
  Z[index_att, q] <- Z[index_att, q] + matrix(1*sqrt((ss$d^2)/(n-1))[q],nout,length(q),byrow = T)
  NB <- c(which(TL==1), sample(which(TL!=1),nin,replace=F))
  Y <- Z %*% t(V)
  Y=Y[NB,]
  TL <- rep(c(1,0),c(nout,nin))
  
  UP <- Y%*%V%*%solve(D)
  apply(UP*sqrt(n-1),2,sd)
  
  sigcomp <- apply(UP*sqrt(n-1),2,sd) > quantile(as.vector(BT2),size)
  
  # UU <- svd(Y)
  # plot(UU$u[,1],UU$u[,2],col=TL+1)
  
  if(sum(sigcomp)>0){
    SCprop <- rowSums((n-1)*UP[,sigcomp,drop=F]^2)
  } else {
    SCprop <- rowSums((n-1)*UP^2)
  }
  ROCprop <- simple_roc(labels=(TL==1), scores=SCprop)
  
  SCwt <- rowSums((n-1)*sweep(UP,2,apply(UP*sqrt(n-1),2,sd),FUN="*")^2)
  ROCwt <- simple_roc(labels=(TL==1), scores=SCwt)
  
  
  YB <- Y%*%V[,st:en,drop=F]%*%t(V[,st:en,drop=F])
  
  SCwang <- rowSums((Y-YB)^2)
  ROCwang <- simple_roc(labels=(TL==1), scores=SCwang)
  
  mhd <- mahalanobis(Y, rep(0,p), cov=SS)
  ROCmhd <- simple_roc(labels=(TL==1), scores=mhd)
  
  a=data.frame(ROCprop,ROCwang,ROCwt,ROCmhd)
  sumo=sumo+a
}
sum0=sumo/nrep

#############################################
sumo=data.frame(matrix(0,nin+nout,12))

for (i in 1:nrep) 
{
  
  Z <- X%*%V
  index_att <- sample(1:n, nout)
  TL <- numeric(n)
  TL[index_att] <- 1
  Z[index_att, q] <- Z[index_att, q] + matrix(2*sqrt((ss$d^2)/(n-1))[q],nout,length(q),byrow = T)
  NB <- c(which(TL==1), sample(which(TL!=1),nin,replace=F))
  Y <- Z %*% t(V)
  Y=Y[NB,]
  TL <- rep(c(1,0),c(nout,nin))
  
  UP <- Y%*%V%*%solve(D)
  apply(UP*sqrt(n-1),2,sd)
  
  sigcomp <- apply(UP*sqrt(n-1),2,sd) > quantile(as.vector(BT2),size)
  
  # UU <- svd(Y)
  # plot(UU$u[,1],UU$u[,2],col=TL+1)
  
  if(sum(sigcomp)>0){
    SCprop <- rowSums((n-1)*UP[,sigcomp,drop=F]^2)
  } else {
    SCprop <- rowSums((n-1)*UP^2)
  }
  ROCprop <- simple_roc(labels=(TL==1), scores=SCprop)
  
  SCwt <- rowSums((n-1)*sweep(UP,2,apply(UP*sqrt(n-1),2,sd),FUN="*")^2)
  ROCwt <- simple_roc(labels=(TL==1), scores=SCwt)
  
  
  YB <- Y%*%V[,st:en,drop=F]%*%t(V[,st:en,drop=F])
  
  SCwang <- rowSums((Y-YB)^2)
  ROCwang <- simple_roc(labels=(TL==1), scores=SCwang)
  
  mhd <- mahalanobis(Y, rep(0,p), cov=SS)
  ROCmhd <- simple_roc(labels=(TL==1), scores=mhd)
  
  a=data.frame(ROCprop,ROCwang,ROCwt,ROCmhd)
  sumo=sumo+a
}
sum1=sumo/nrep

#############################################
sumo=data.frame(matrix(0,nin+nout,12))

for (i in 1:nrep) 
{
  
  Z <- X%*%V
  index_att <- sample(1:n, nout)
  TL <- numeric(n)
  TL[index_att] <- 1
  Z[index_att, q] <- Z[index_att, q] + matrix(3*sqrt((ss$d^2)/(n-1))[q],nout,length(q),byrow = T)
  NB <- c(which(TL==1), sample(which(TL!=1),nin,replace=F))
  Y <- Z %*% t(V)
  Y=Y[NB,]
  TL <- rep(c(1,0),c(nout,nin))
  
  UP <- Y%*%V%*%solve(D)
  apply(UP*sqrt(n-1),2,sd)
  
  sigcomp <- apply(UP*sqrt(n-1),2,sd) > quantile(as.vector(BT2),size)
  
  # UU <- svd(Y)
  # plot(UU$u[,1],UU$u[,2],col=TL+1)
  
  if(sum(sigcomp)>0){
    SCprop <- rowSums((n-1)*UP[,sigcomp,drop=F]^2)
  } else {
    SCprop <- rowSums((n-1)*UP^2)
  }
  ROCprop <- simple_roc(labels=(TL==1), scores=SCprop)
  
  SCwt <- rowSums((n-1)*sweep(UP,2,apply(UP*sqrt(n-1),2,sd),FUN="*")^2)
  ROCwt <- simple_roc(labels=(TL==1), scores=SCwt)
  
  
  YB <- Y%*%V[,st:en,drop=F]%*%t(V[,st:en,drop=F])
  
  SCwang <- rowSums((Y-YB)^2)
  ROCwang <- simple_roc(labels=(TL==1), scores=SCwang)
  
  mhd <- mahalanobis(Y, rep(0,p), cov=SS)
  ROCmhd <- simple_roc(labels=(TL==1), scores=mhd)
  
  a=data.frame(ROCprop,ROCwang,ROCwt,ROCmhd)
  sumo=sumo+a
}
sum2=sumo/nrep

pdf(file="Fig1.pdf",width=10, height=14)
par(mfrow=c(3,1), mar=c(6, 6, 4, 2))

plot(sum0[,11],sum0[,10],type="l",col=1,lwd=2,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum0[,2],sum0[,1],type="l",col=2,lwd=2,lty=2)
lines(sum0[,8],sum0[,7],type="l",col=3,lwd=2,lty=3)
lines(sum0[,5],sum0[,4],type="l",col=4,lwd=2,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(a)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2.5, line=3.5)
legend("bottomright", legend=c("True","AAD","WAAD", "WBPCA"), 
       lty=c(1,2,3,4), col=c(1,2,3,4), 
       cex = 2,lwd=c(2,2,2,2))

plot(sum1[,11],sum1[,10],type="l",col=1,lwd=2,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum1[,2],sum1[,1],type="l",col=2,lwd=2,lty=2)
lines(sum1[,8],sum1[,7],type="l",col=3,lwd=2,lty=3)
lines(sum1[,5],sum1[,4],type="l",col=4,lwd=2,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(b)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2.5, line=3.5)


plot(sum2[,11],sum2[,10],type="l",col=1,lwd=2,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum2[,2],sum2[,1],type="l",col=2,lwd=2,lty=2)
lines(sum2[,8],sum2[,7],type="l",col=3,lwd=2,lty=3)
lines(sum2[,5],sum2[,4],type="l",col=4,lwd=2,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(c)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2.5, line=3.5)

dev.off()

setwd("C:\\Users\\Admin\\Desktop\\IEEEtran")

postscript(file="Fig1a.eps")
par(mar=c(6, 7, 4, 2))
plot(sum0[,11],sum0[,10],type="l",col=1,lwd=3,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum0[,2],sum0[,1],type="l",col=2,lwd=3,lty=2)
lines(sum0[,8],sum0[,7],type="l",col=3,lwd=3,lty=3)
lines(sum0[,5],sum0[,4],type="l",col=4,lwd=3,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(a)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2, line=5)
legend("bottomright", legend=c("True","AAD","WAAD", "WBPCA"), 
       lty=c(1,2,3,4), col=c(1,2,3,4), 
       cex = 2,lwd=c(3,3,3,3))
dev.off()

postscript(file="Fig1b.eps")
par(mar=c(6, 7, 4, 2))
plot(sum1[,11],sum1[,10],type="l",col=1,lwd=3,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum1[,2],sum1[,1],type="l",col=2,lwd=3,lty=2)
lines(sum1[,8],sum1[,7],type="l",col=3,lwd=3,lty=3)
lines(sum1[,5],sum1[,4],type="l",col=4,lwd=3,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(b)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2, line=5)
dev.off()

postscript(file="Fig1c.eps")
par(mar=c(6, 7, 4, 2))
plot(sum2[,11],sum2[,10],type="l",col=1,lwd=3,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum2[,2],sum2[,1],type="l",col=2,lwd=3,lty=2)
lines(sum2[,8],sum2[,7],type="l",col=3,lwd=3,lty=3)
lines(sum2[,5],sum2[,4],type="l",col=4,lwd=3,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(c)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2, line=5)
dev.off()
##########################################################################################
##########################################################################################
##########################################################################################

sumo=data.frame(matrix(0,nin+nout,12))
# components along which we want to make shift
q <- c(p-2,p-1,p)

for (i in 1:nrep) 
{
  
  Z <- X%*%V
  index_att <- sample(1:n, nout)
  TL <- numeric(n)
  TL[index_att] <- 1
  Z[index_att, q] <- Z[index_att, q] + matrix(1*sqrt((ss$d^2)/(n-1))[q],nout,length(q),byrow = T)
  NB <- c(which(TL==1), sample(which(TL!=1),nin,replace=F))
  Y <- Z %*% t(V)
  Y=Y[NB,]
  TL <- rep(c(1,0),c(nout,nin))
  
  UP <- Y%*%V%*%solve(D)
  apply(UP*sqrt(n-1),2,sd)
  
  sigcomp <- apply(UP*sqrt(n-1),2,sd) > quantile(as.vector(BT2),size)
  
  # UU <- svd(Y)
  # plot(UU$u[,1],UU$u[,2],col=TL+1)
  
  if(sum(sigcomp)>0){
    SCprop <- rowSums((n-1)*UP[,sigcomp,drop=F]^2)
  } else {
    SCprop <- rowSums((n-1)*UP^2)
  }
  ROCprop <- simple_roc(labels=(TL==1), scores=SCprop)
  
  SCwt <- rowSums((n-1)*sweep(UP,2,apply(UP*sqrt(n-1),2,sd),FUN="*")^2)
  ROCwt <- simple_roc(labels=(TL==1), scores=SCwt)
  
  
  YB <- Y%*%V[,st:en,drop=F]%*%t(V[,st:en,drop=F])
  
  SCwang <- rowSums((Y-YB)^2)
  ROCwang <- simple_roc(labels=(TL==1), scores=SCwang)
  
  mhd <- mahalanobis(Y, rep(0,p), cov=SS)
  ROCmhd <- simple_roc(labels=(TL==1), scores=mhd)
  
  a=data.frame(ROCprop,ROCwang,ROCwt,ROCmhd)
  sumo=sumo+a
}
sum0=sumo/nrep

#############################################
sumo=data.frame(matrix(0,nin+nout,12))

for (i in 1:nrep) 
{
  
  Z <- X%*%V
  index_att <- sample(1:n, nout)
  TL <- numeric(n)
  TL[index_att] <- 1
  Z[index_att, q] <- Z[index_att, q] + matrix(2*sqrt((ss$d^2)/(n-1))[q],nout,length(q),byrow = T)
  NB <- c(which(TL==1), sample(which(TL!=1),nin,replace=F))
  Y <- Z %*% t(V)
  Y=Y[NB,]
  TL <- rep(c(1,0),c(nout,nin))
  
  UP <- Y%*%V%*%solve(D)
  apply(UP*sqrt(n-1),2,sd)
  
  sigcomp <- apply(UP*sqrt(n-1),2,sd) > quantile(as.vector(BT2),size)
  
  # UU <- svd(Y)
  # plot(UU$u[,1],UU$u[,2],col=TL+1)
  
  if(sum(sigcomp)>0){
    SCprop <- rowSums((n-1)*UP[,sigcomp,drop=F]^2)
  } else {
    SCprop <- rowSums((n-1)*UP^2)
  }
  ROCprop <- simple_roc(labels=(TL==1), scores=SCprop)
  
  SCwt <- rowSums((n-1)*sweep(UP,2,apply(UP*sqrt(n-1),2,sd),FUN="*")^2)
  ROCwt <- simple_roc(labels=(TL==1), scores=SCwt)
  
  
  YB <- Y%*%V[,st:en,drop=F]%*%t(V[,st:en,drop=F])
  
  SCwang <- rowSums((Y-YB)^2)
  ROCwang <- simple_roc(labels=(TL==1), scores=SCwang)
  
  mhd <- mahalanobis(Y, rep(0,p), cov=SS)
  ROCmhd <- simple_roc(labels=(TL==1), scores=mhd)
  
  a=data.frame(ROCprop,ROCwang,ROCwt,ROCmhd)
  sumo=sumo+a
}
sum1=sumo/nrep

#############################################
sumo=data.frame(matrix(0,nin+nout,12))

for (i in 1:nrep) 
{
  
  Z <- X%*%V
  index_att <- sample(1:n, nout)
  TL <- numeric(n)
  TL[index_att] <- 1
  Z[index_att, q] <- Z[index_att, q] + matrix(3*sqrt((ss$d^2)/(n-1))[q],nout,length(q),byrow = T)
  NB <- c(which(TL==1), sample(which(TL!=1),nin,replace=F))
  Y <- Z %*% t(V)
  Y=Y[NB,]
  TL <- rep(c(1,0),c(nout,nin))
  
  UP <- Y%*%V%*%solve(D)
  apply(UP*sqrt(n-1),2,sd)
  
  sigcomp <- apply(UP*sqrt(n-1),2,sd) > quantile(as.vector(BT2),size)
  
  # UU <- svd(Y)
  # plot(UU$u[,1],UU$u[,2],col=TL+1)
  
  if(sum(sigcomp)>0){
    SCprop <- rowSums((n-1)*UP[,sigcomp,drop=F]^2)
  } else {
    SCprop <- rowSums((n-1)*UP^2)
  }
  ROCprop <- simple_roc(labels=(TL==1), scores=SCprop)
  
  SCwt <- rowSums((n-1)*sweep(UP,2,apply(UP*sqrt(n-1),2,sd),FUN="*")^2)
  ROCwt <- simple_roc(labels=(TL==1), scores=SCwt)
  
  
  YB <- Y%*%V[,st:en,drop=F]%*%t(V[,st:en,drop=F])
  
  SCwang <- rowSums((Y-YB)^2)
  ROCwang <- simple_roc(labels=(TL==1), scores=SCwang)
  
  mhd <- mahalanobis(Y, rep(0,p), cov=SS)
  ROCmhd <- simple_roc(labels=(TL==1), scores=mhd)
  
  a=data.frame(ROCprop,ROCwang,ROCwt,ROCmhd)
  sumo=sumo+a
}
sum2=sumo/nrep

pdf(file="Fig2.pdf",width=10, height=14)
par(mfrow=c(3,1), mar=c(6, 6, 4, 2))

plot(sum0[,11],sum0[,10],type="l",col=1,lwd=2,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum0[,2],sum0[,1],type="l",col=2,lwd=2,lty=2)
lines(sum0[,8],sum0[,7],type="l",col=3,lwd=2,lty=3)
lines(sum0[,5],sum0[,4],type="l",col=4,lwd=2,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(a)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2.5, line=3.5)
legend("bottomright", legend=c("True","AAD","WAAD", "WBPCA"), 
       lty=c(1,2,3,4), col=c(1,2,3,4), 
       cex = 2,lwd=c(2,2,2,2))

plot(sum1[,11],sum1[,10],type="l",col=1,lwd=2,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum1[,2],sum1[,1],type="l",col=2,lwd=2,lty=2)
lines(sum1[,8],sum1[,7],type="l",col=3,lwd=2,lty=3)
lines(sum1[,5],sum1[,4],type="l",col=4,lwd=2,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(b)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2.5, line=3.5)


plot(sum2[,11],sum2[,10],type="l",col=1,lwd=2,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum2[,2],sum2[,1],type="l",col=2,lwd=2,lty=2)
lines(sum2[,8],sum2[,7],type="l",col=3,lwd=2,lty=3)
lines(sum2[,5],sum2[,4],type="l",col=4,lwd=2,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(c)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2.5, line=3.5)

dev.off()


postscript(file="Fig2a.eps")
par(mar=c(6, 7, 4, 2))
plot(sum0[,11],sum0[,10],type="l",col=1,lwd=3,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum0[,2],sum0[,1],type="l",col=2,lwd=3,lty=2)
lines(sum0[,8],sum0[,7],type="l",col=3,lwd=3,lty=3)
lines(sum0[,5],sum0[,4],type="l",col=4,lwd=3,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(a)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2, line=5)
legend("bottomright", legend=c("True","AAD","WAAD", "WBPCA"), 
       lty=c(1,2,3,4), col=c(1,2,3,4), 
       cex = 2,lwd=c(3,3,3,3))
dev.off()

postscript(file="Fig2b.eps")
par(mar=c(6, 7, 4, 2))
plot(sum1[,11],sum1[,10],type="l",col=1,lwd=3,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum1[,2],sum1[,1],type="l",col=2,lwd=3,lty=2)
lines(sum1[,8],sum1[,7],type="l",col=3,lwd=3,lty=3)
lines(sum1[,5],sum1[,4],type="l",col=4,lwd=3,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(b)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2, line=5)
dev.off()

postscript(file="Fig2c.eps")
par(mar=c(6, 7, 4, 2))
plot(sum2[,11],sum2[,10],type="l",col=1,lwd=3,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum2[,2],sum2[,1],type="l",col=2,lwd=3,lty=2)
lines(sum2[,8],sum2[,7],type="l",col=3,lwd=3,lty=3)
lines(sum2[,5],sum2[,4],type="l",col=4,lwd=3,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(c)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2, line=5)
dev.off()
##########################################################################################
##########################################################################################
##########################################################################################

sumo=data.frame(matrix(0,nin+nout,12))
# components along which we want to make shift
q <- c(1,2,3)

for (i in 1:nrep) 
{
  
  Z <- X%*%V
  index_att <- sample(1:n, nout)
  TL <- numeric(n)
  TL[index_att] <- 1
  Z[index_att, q] <- Z[index_att, q] + matrix(1*sqrt((ss$d^2)/(n-1))[q],nout,length(q),byrow = T)
  NB <- c(which(TL==1), sample(which(TL!=1),nin,replace=F))
  Y <- Z %*% t(V)
  Y=Y[NB,]
  TL <- rep(c(1,0),c(nout,nin))
  
  UP <- Y%*%V%*%solve(D)
  apply(UP*sqrt(n-1),2,sd)
  
  sigcomp <- apply(UP*sqrt(n-1),2,sd) > quantile(as.vector(BT2),size)
  
  # UU <- svd(Y)
  # plot(UU$u[,1],UU$u[,2],col=TL+1)
  
  if(sum(sigcomp)>0){
    SCprop <- rowSums((n-1)*UP[,sigcomp,drop=F]^2)
  } else {
    SCprop <- rowSums((n-1)*UP^2)
  }
  ROCprop <- simple_roc(labels=(TL==1), scores=SCprop)
  
  SCwt <- rowSums((n-1)*sweep(UP,2,apply(UP*sqrt(n-1),2,sd),FUN="*")^2)
  ROCwt <- simple_roc(labels=(TL==1), scores=SCwt)
  
  
  YB <- Y%*%V[,st:en,drop=F]%*%t(V[,st:en,drop=F])
  
  SCwang <- rowSums((Y-YB)^2)
  ROCwang <- simple_roc(labels=(TL==1), scores=SCwang)
  
  mhd <- mahalanobis(Y, rep(0,p), cov=SS)
  ROCmhd <- simple_roc(labels=(TL==1), scores=mhd)
  
  a=data.frame(ROCprop,ROCwang,ROCwt,ROCmhd)
  sumo=sumo+a
}
sum0=sumo/nrep

#############################################
sumo=data.frame(matrix(0,nin+nout,12))

for (i in 1:nrep) 
{
  
  Z <- X%*%V
  index_att <- sample(1:n, nout)
  TL <- numeric(n)
  TL[index_att] <- 1
  Z[index_att, q] <- Z[index_att, q] + matrix(2*sqrt((ss$d^2)/(n-1))[q],nout,length(q),byrow = T)
  NB <- c(which(TL==1), sample(which(TL!=1),nin,replace=F))
  Y <- Z %*% t(V)
  Y=Y[NB,]
  TL <- rep(c(1,0),c(nout,nin))
  
  UP <- Y%*%V%*%solve(D)
  apply(UP*sqrt(n-1),2,sd)
  
  sigcomp <- apply(UP*sqrt(n-1),2,sd) > quantile(as.vector(BT2),size)
  
  # UU <- svd(Y)
  # plot(UU$u[,1],UU$u[,2],col=TL+1)
  
  if(sum(sigcomp)>0){
    SCprop <- rowSums((n-1)*UP[,sigcomp,drop=F]^2)
  } else {
    SCprop <- rowSums((n-1)*UP^2)
  }
  ROCprop <- simple_roc(labels=(TL==1), scores=SCprop)
  
  SCwt <- rowSums((n-1)*sweep(UP,2,apply(UP*sqrt(n-1),2,sd),FUN="*")^2)
  ROCwt <- simple_roc(labels=(TL==1), scores=SCwt)
  
  
  YB <- Y%*%V[,st:en,drop=F]%*%t(V[,st:en,drop=F])
  
  SCwang <- rowSums((Y-YB)^2)
  ROCwang <- simple_roc(labels=(TL==1), scores=SCwang)
  
  mhd <- mahalanobis(Y, rep(0,p), cov=SS)
  ROCmhd <- simple_roc(labels=(TL==1), scores=mhd)
  
  a=data.frame(ROCprop,ROCwang,ROCwt,ROCmhd)
  sumo=sumo+a
}
sum1=sumo/nrep

#############################################
sumo=data.frame(matrix(0,nin+nout,12))

for (i in 1:nrep) 
{
  
  Z <- X%*%V
  index_att <- sample(1:n, nout)
  TL <- numeric(n)
  TL[index_att] <- 1
  Z[index_att, q] <- Z[index_att, q] + matrix(3*sqrt((ss$d^2)/(n-1))[q],nout,length(q),byrow = T)
  NB <- c(which(TL==1), sample(which(TL!=1),nin,replace=F))
  Y <- Z %*% t(V)
  Y=Y[NB,]
  TL <- rep(c(1,0),c(nout,nin))
  
  UP <- Y%*%V%*%solve(D)
  apply(UP*sqrt(n-1),2,sd)
  
  sigcomp <- apply(UP*sqrt(n-1),2,sd) > quantile(as.vector(BT2),size)
  
  # UU <- svd(Y)
  # plot(UU$u[,1],UU$u[,2],col=TL+1)
  
  if(sum(sigcomp)>0){
    SCprop <- rowSums((n-1)*UP[,sigcomp,drop=F]^2)
  } else {
    SCprop <- rowSums((n-1)*UP^2)
  }
  ROCprop <- simple_roc(labels=(TL==1), scores=SCprop)
  
  SCwt <- rowSums((n-1)*sweep(UP,2,apply(UP*sqrt(n-1),2,sd),FUN="*")^2)
  ROCwt <- simple_roc(labels=(TL==1), scores=SCwt)
  
  
  YB <- Y%*%V[,st:en,drop=F]%*%t(V[,st:en,drop=F])
  
  SCwang <- rowSums((Y-YB)^2)
  ROCwang <- simple_roc(labels=(TL==1), scores=SCwang)
  
  mhd <- mahalanobis(Y, rep(0,p), cov=SS)
  ROCmhd <- simple_roc(labels=(TL==1), scores=mhd)
  
  a=data.frame(ROCprop,ROCwang,ROCwt,ROCmhd)
  sumo=sumo+a
}
sum2=sumo/nrep

pdf(file="Fig3.pdf",width=10, height=14)
par(mfrow=c(3,1), mar=c(6, 6, 4, 2))

plot(sum0[,11],sum0[,10],type="l",col=1,lwd=2,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum0[,2],sum0[,1],type="l",col=2,lwd=2,lty=2)
lines(sum0[,8],sum0[,7],type="l",col=3,lwd=2,lty=3)
lines(sum0[,5],sum0[,4],type="l",col=4,lwd=2,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(a)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2.5, line=3.5)
legend("bottomright", legend=c("True","AAD","WAAD", "WBPCA"), 
       lty=c(1,2,3,4), col=c(1,2,3,4), 
       cex = 2,lwd=c(2,2,2,2))

plot(sum1[,11],sum1[,10],type="l",col=1,lwd=2,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum1[,2],sum1[,1],type="l",col=2,lwd=2,lty=2)
lines(sum1[,8],sum1[,7],type="l",col=3,lwd=2,lty=3)
lines(sum1[,5],sum1[,4],type="l",col=4,lwd=2,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(b)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2.5, line=3.5)


plot(sum2[,11],sum2[,10],type="l",col=1,lwd=2,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum2[,2],sum2[,1],type="l",col=2,lwd=2,lty=2)
lines(sum2[,8],sum2[,7],type="l",col=3,lwd=2,lty=3)
lines(sum2[,5],sum2[,4],type="l",col=4,lwd=2,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(c)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2.5, line=3.5)

dev.off()




postscript(file="Fig3a.eps")
par(mar=c(6, 7, 4, 2))
plot(sum0[,11],sum0[,10],type="l",col=1,lwd=3,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum0[,2],sum0[,1],type="l",col=2,lwd=3,lty=2)
lines(sum0[,8],sum0[,7],type="l",col=3,lwd=3,lty=3)
lines(sum0[,5],sum0[,4],type="l",col=4,lwd=3,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(a)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2, line=5)
legend("bottomright", legend=c("True","AAD","WAAD", "WBPCA"), 
       lty=c(1,2,3,4), col=c(1,2,3,4), 
       cex = 2,lwd=c(3,3,3,3))
dev.off()

postscript(file="Fig3b.eps")
par(mar=c(6, 7, 4, 2))
plot(sum1[,11],sum1[,10],type="l",col=1,lwd=3,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum1[,2],sum1[,1],type="l",col=2,lwd=3,lty=2)
lines(sum1[,8],sum1[,7],type="l",col=3,lwd=3,lty=3)
lines(sum1[,5],sum1[,4],type="l",col=4,lwd=3,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(b)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2, line=5)
dev.off()

postscript(file="Fig3c.eps")
par(mar=c(6, 7, 4, 2))
plot(sum2[,11],sum2[,10],type="l",col=1,lwd=3,xlab="", ylab="", 
     cex.axis=2, las=1,ylim=c(0,1),xlim=c(0,1))

lines(sum2[,2],sum2[,1],type="l",col=2,lwd=3,lty=2)
lines(sum2[,8],sum2[,7],type="l",col=3,lwd=3,lty=3)
lines(sum2[,5],sum2[,4],type="l",col=4,lwd=3,lty=4)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(c)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=2, line = 3.5)
title(ylab="True Positive Rate", cex.lab=2, line=5)
dev.off()