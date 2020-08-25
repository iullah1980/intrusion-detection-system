
# specify a path where you need to save the results
setwd("C:\\KDD")

library(data.table)
library(plyr)
library(expm)


#KDDfull <- read.csv('C:/Users/Admin/Downloads/kddcup.data/kddcup.data.corrected',header=F)
tmp <- tempfile()
download.file("http://kdd.ics.uci.edu/databases/kddcup99/kddcup.data.gz",tmp)
KDDfull <- read.csv(gzfile(tmp),header=F)

KDDfullup <- KDDfull#KDDfull[-c(352754,409970,420131,1401653,1415576,1415577,1417572,3421756,3666531,3927298,3927301,3927302,3927303,3927304),]
cnms <- read.csv('http://kdd.ics.uci.edu/databases/kddcup99/kddcup.names', skip = 1, header=FALSE, sep = ":")

att <- ifelse(KDDfullup$V42=="normal.",0,1)
KDDf <- KDDfullup[,-c(which(cnms$V2==" symbolic."),42)]

KDDnorm <- KDDfullup[KDDfullup$V42=="normal.",-c(which(cnms$V2==" symbolic."),42)]
KDDinf <- KDDfullup[KDDfullup$V42!="normal.",-c(which(cnms$V2==" symbolic."),42)]

colnames(KDDnorm)=cnms$V1[cnms$V2==" continuous."]
colnames(KDDinf)=cnms$V1[cnms$V2==" continuous."]
colnames(KDDf)=cnms$V1[cnms$V2==" continuous."]

# To remove columns than have very small number of non-zero values
which(apply(KDDf!=0,2,sum)<6000)
which(apply(KDDnorm!=0,2,sum)<6000)
which(apply(KDDinf!=0,2,sum)<6000)
# columns 4,5,6,7,8,9,10,11,12,13, and 15 need to be excluded from analysis

# levels(KDDfullup$V42)
# table(KDDfullup$V42)


cleandata <- KDDnorm[,-c(4,5,6,7,8,9,10,11,12,13,15)]
s_cleandata <- scale(cleandata)

ss <- svd(s_cleandata)
#U=ss$u
D=diag(ss$d)
V=ss$v

############################################
BTdelta <- replicate(5000,{
  XB <-  s_cleandata[sample(1:dim(s_cleandata)[1],30000,replace=TRUE),]
  UPB <- XB%*%V%*%diag(1/ss$d)
  apply(UPB*sqrt(dim(s_cleandata)[1]-1),2,sd)
})

delta0=apply(BTdelta,1,quantile,probs = c(0.01,0.99))
hist(BTdelta)
boxplot(BTdelta~row(BTdelta),outline=F)
abline(h=1)


png(file="KDD99outliers.png")
par(mar=c(6, 6, 4, 2))
boxplot(BTdelta~row(BTdelta),outline=F, xlab="", ylab="", cex.axis=1.5, las=1, main="")
abline(h=1,col="gray",lwd=1)
title(xlab="Index", cex.lab=1.5, line = 3.5)
title(ylab=bquote(paste(s[j]^u)), cex.lab=1.5, line=3.5)
dev.off()

#############################################
which(apply(BTdelta,1,quantile,0.75) < 1)
abs(V[,which(apply(BTdelta,1,quantile,0.75) < 1)]) > .3
V[,which(apply(BTdelta,1,quantile,0.75) < 1)]

names(cleandata)

hist(cleandata[,3])


out=c(which(cleandata[,2] > 3000000),
      which(cleandata[,3] > 3000000)
)
cleandata <- cleandata[-out,]
s_cleandata <- scale(cleandata)

ss <- svd(s_cleandata)
#U=ss$u
D=diag(ss$d)
V=ss$v

############################################
BTdelta <- replicate(5000,{
  XB <-  s_cleandata[sample(1:dim(s_cleandata)[1],30000,replace=TRUE),]
  UPB <- XB%*%V%*%diag(1/ss$d)
  apply(UPB*sqrt(dim(s_cleandata)[1]-1),2,sd)
})

delta0=apply(BTdelta,1,quantile,probs = c(0.01,0.99))

hist(BTdelta)
boxplot(BTdelta~row(BTdelta),outline=F)
abline(h=1)

png(file="KDD99wooutliers.png")
par(mar=c(6, 6, 4, 2))
boxplot(BTdelta~row(BTdelta),outline=F, xlab="", ylab="", cex.axis=1.5, las=1, main="")
abline(h=1,col="gray",lwd=1)
title(xlab="Index", cex.lab=1.5, line = 3.5)
title(ylab=bquote(paste(s[j]^u)), cex.lab=1.5, line=3.5)
dev.off()

#############################################
apply(BTdelta,1,quantile,0.75)
which(apply(BTdelta,1,quantile,0.75) < 1)
abs(V[,which(apply(BTdelta,1,quantile,0.75) < 1)]) > .2


# levels(KDDfullup$V42)
# table(KDDfullup$V42)

KDDinfmod <- KDDinf[,-c(4,5,6,7,8,9,10,11,12,13,15)]


nout=100
nin=9900
n=nin+nout

res <- replicate(1000,{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=T),],
                       KDDinfmod[sample(1:dim(KDDinfmod)[1],nout,replace=T),]), 
                 center = attr(s_cleandata, 'scaled:center'),scale=attr(s_cleandata, 'scaled:scale'))
  
  svds2 <- NData%*%V%*%diag(1/ss$d)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  sdsvds2 > delta0[2,]
})

barplot(apply(res,1,sum)/1000)

barplot(table(apply(res,2,sum)))

#######################################################
png(file="barplotms.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(table(apply(res,2,sum))/sum(table(apply(res,2,sum))),
                   ylim=c(0,max(table(apply(res,2,sum))/sum(table(apply(res,2,sum))))+.02),
                   col="gray",las=1,ylab="",cex.lab=1.5,cex.axis=1.5)
axis(1,barCenters, labels = FALSE)
title(ylab = "Proportion", cex.lab=1.5, line=4)
title(xlab="Number of components", cex.lab=1.5)
title(main="(a)",cex.main=2)
box()
#mtext("(a)",at=1, line = 1, cex =2)
dev.off()

png(file="barplotall.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(apply(res,1,sum)/1000*100,ylim=c(0,105),col="gray",
                   las=1,ylab="Percentage",cex.lab=1.5,cex.axis=1.5)
text(x = barCenters[seq(1,32,2)], y = par("usr")[3] - 3, srt = 45,
     adj = 1, labels = paste0("Comp-",seq(1,32,2)), xpd = TRUE,cex.lab=2.5)
axis(1,barCenters, labels = FALSE)
title(main="(b)",cex.main=2)
box()
#mtext("(b)",at=1, line = 1, cex =2)
dev.off()

# boxplot(res~row(res),outline=F)
# abline(h=1)

###################################################
# a function to find ROC curve
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

nout=100
nin=9900
n=nin+nout
sumo=data.frame(matrix(0,n,9))

st <- 1
en <- eigen(cov(s_cleandata))$values>1



for (i in 1:1000) 
{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=T),],
                       KDDinfmod[sample(1:dim(KDDinfmod)[1],nout,replace=T),]), 
                 center = attr(s_cleandata, 'scaled:center'),scale=attr(s_cleandata, 'scaled:scale'))
  
  att=rep(c(0,1),c(nin,nout))
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  SCaffected <- rowSums((dim(s_cleandata)[1]-1)*svds2[,sdsvds2 > delta0[2,],drop=F]^2)#
  ROCaffected <- simple_roc(labels=(att==1), scores=SCaffected)
  
  SCw <- rowSums((dim(s_cleandata)[1]-1)*(sweep(svds2,2,sdsvds2,FUN="*")^2))
  ROCw <- simple_roc(labels=(att==1), scores=SCw)
  
  
  XB <- NData%*%V[,en,drop=F]%*%t(V[,en,drop=F])
  SCWBPCA <- rowSums((NData-XB)^2)
  ROCWBPCA <- simple_roc(labels=(att==1), scores=SCWBPCA)
  
  
  a=data.frame(ROCaffected,ROCw,ROCWBPCA)
  sumo=sumo+a
}
sum0=sumo/1000

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", 
     ylab="", cex.axis=2,xlim=c(0,.05))
lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)

#############################################

png(file="kdd99roc.png")
par(mar=c(6, 7, 4, 2))

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=1.5,xlim=c(0,.05),las=1)

lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(c)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=1.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=1.5, line=3.5)
legend("bottomright", legend=c("AAD", "WBPCA","WAAD"), 
       lty=c(1,2,3), col=c(1,1,1), 
       cex = 1.5,lwd=c(2,2,2))

dev.off()


#####################################################################
#####################################################################
#####################################################################

# BTdelta <- replicate(5000,{
#   XB <-  s_cleandata[sample(1:dim(s_cleandata)[1],30000,replace=TRUE),]
#   UPB <- XB%*%V%*%diag(1/ss$d)
#   apply(UPB*sqrt(dim(s_cleandata)[1]-1),2,sd)
# })
# 
# delta0=apply(BTdelta,1,quantile,probs = c(0.05,0.95))
# 
# hist(BTdelta)
# boxplot(BTdelta~row(BTdelta),outline=F)
# abline(h=1)

levels(KDDfullup$V42)
table(KDDfullup$V42)
KDDinf <- KDDfullup[!KDDfullup$V42%in%c("back.","ipsweep.","neptune.","nmap.",
                                        "normal.","portsweep.","satan.","smurf."),
                    -c(which(cnms$V2==" symbolic."),42)]

# KDDinf <- KDDfullup[KDDfullup$V42%in%c("teardrop.","teardrop.","pod.",
#                                        "guess_passwd."),
#                     -c(which(cnms$V2==" symbolic."),42)]

colnames(KDDinf)=cnms$V1[cnms$V2==" continuous."]
KDDinfmod <- KDDinf[,-c(4,5,6,7,8,9,10,11,12,13,15)]
dim(KDDinfmod)
  
nout=100
nin=10000
n=nin+nout

res <- replicate(1000,{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=T),],
                       KDDinfmod[sample(1:dim(KDDinfmod)[1],nout,replace=T),]),  
                 center = attr(s_cleandata, 'scaled:center'),scale=attr(s_cleandata, 'scaled:scale'))
  
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  sdsvds2 > delta0[2,]
})

barplot(apply(res,1,sum)/1000)

barplot(table(apply(res,2,sum)))

#######################################################
png(file="barplotmsrare.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(table(apply(res,2,sum))/sum(table(apply(res,2,sum))),
                   ylim=c(0,max(table(apply(res,2,sum))/sum(table(apply(res,2,sum))))+.02),
                   col="gray",las=1,ylab="",cex.lab=1.5,cex.axis=1.5)
axis(1,barCenters, labels = FALSE)
title(ylab = "Proportion", cex.lab=1.5, line=4)
title(xlab="Number of components", cex.lab=1.5)
title(main="(a)",cex.main=2)
box()
#mtext("(a)",at=1, line = 1, cex =2)
dev.off()

png(file="barplotrare.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(apply(res,1,sum)/1000*100,ylim=c(0,105),col="gray",
                   las=1,ylab="Percentage",cex.lab=1.5,cex.axis=1.5)
text(x = barCenters[seq(1,32,2)], y = par("usr")[3] - 3, srt = 45,
     adj = 1, labels = paste0("Comp-",seq(1,32,2)), xpd = TRUE,cex.lab=2.5)
axis(1,barCenters, labels = FALSE)
title(main="(b)",cex.main=2)
box()
#mtext("(b)",at=1, line = 1, cex =2)
dev.off()

# boxplot(res~row(res),outline=F)
# abline(h=1)

###################################################
nout=100
nin=9900
n=nin+nout
sumo=data.frame(matrix(0,n,9))

st <- 1
en <- eigen(cov(s_cleandata))$values>1



for (i in 1:1000) 
{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=T),],
                       KDDinfmod[sample(1:dim(KDDinfmod)[1],nout,replace=T),]),  
                 center = attr(s_cleandata, 'scaled:center'),scale=attr(s_cleandata, 'scaled:scale'))
  
  att=rep(c(0,1),c(nin,nout))
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  SCaffected <- rowSums((dim(s_cleandata)[1]-1)*svds2[,sdsvds2 > delta0[2,],drop=F]^2)#
  ROCaffected <- simple_roc(labels=(att==1), scores=SCaffected)
  
  SCw <- rowSums((dim(s_cleandata)[1]-1)*(sweep(svds2,2,sdsvds2,FUN="*")^2))
  ROCw <- simple_roc(labels=(att==1), scores=SCw)
  
  
  XB <- NData%*%V[,en,drop=F]%*%t(V[,en,drop=F])
  SCWBPCA <- rowSums((NData-XB)^2)
  ROCWBPCA <- simple_roc(labels=(att==1), scores=SCWBPCA)
  
  # plot(SCw,col=att+1)
  # plot(SCWBPCA,col=att+1)
  # plot(SCw,SCWBPCA,col=att+1)
  
  a=data.frame(ROCaffected,ROCw,ROCWBPCA)
  sumo=sumo+a
}
sum0=sumo/1000

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=2,xlim=c(0,1))
lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)

#############################################

png(file="kdd99rocrare.png")
par(mar=c(6, 7, 4, 2))

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=1.5,xlim=c(0,.4),las=1)

lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(c)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=1.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=1.5, line=3.5)
legend("bottomright", legend=c("AAD", "WBPCA","WAAD"), 
       lty=c(1,2,3), col=c(1,1,1), 
       cex = 1.5,lwd=c(2,2,2))

dev.off()


#####################################################################
#####################################################################
#####################################################################

levels(KDDfullup$V42)
table(KDDfullup$V42)
KDDinf <- KDDfullup[KDDfullup$V42%in%c("smurf."),-c(which(cnms$V2==" symbolic."),42)]
colnames(KDDinf)=cnms$V1[cnms$V2==" continuous."]
KDDinfmod <- KDDinf[,-c(4,5,6,7,8,9,10,11,12,13,15)]


nout=100
nin=9900
n=nin+nout

res <- replicate(1000,{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=T),],
                       KDDinfmod[sample(1:dim(KDDinfmod)[1],nout,replace=T),]),  
                 center = attr(s_cleandata, 'scaled:center'),scale=attr(s_cleandata, 'scaled:scale'))
  
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  sdsvds2 > delta0[2,]
})

barplot(apply(res,1,sum)/1000)

barplot(table(apply(res,2,sum)))

#######################################################
png(file="barplotmssmurf.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(table(apply(res,2,sum))/sum(table(apply(res,2,sum))),
                   ylim=c(0,max(table(apply(res,2,sum))/sum(table(apply(res,2,sum))))+.02),
                   col="gray",las=1,ylab="",cex.lab=1.5,cex.axis=1.5)
axis(1,barCenters, labels = FALSE)
title(ylab = "Proportion", cex.lab=1.5, line=4)
title(xlab="Number of components", cex.lab=1.5)
title(main="(a)",cex.main=2)
box()
#mtext("(a)",at=1, line = 1, cex =2)
dev.off()

png(file="barplotsmurf.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(apply(res,1,sum)/1000*100,ylim=c(0,105),col="gray",
                   las=1,ylab="Percentage",cex.lab=1.5,cex.axis=1.5)
text(x = barCenters[seq(1,32,2)], y = par("usr")[3] - 3, srt = 45,
     adj = 1, labels = paste0("Comp-",seq(1,32,2)), xpd = TRUE,cex.lab=2.5)
axis(1,barCenters, labels = FALSE)
title(main="(b)",cex.main=2)
box()
#mtext("(b)",at=1, line = 1, cex =2)
dev.off()

# boxplot(res~row(res),outline=F)
# abline(h=1)

###################################################
nout=100
nin=9900
n=nin+nout
sumo=data.frame(matrix(0,n,9))

st <- 1
en <- eigen(cov(s_cleandata))$values>1

for (i in 1:1000) 
{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=T),],
                       KDDinfmod[sample(1:dim(KDDinfmod)[1],nout,replace=T),]),
                 center = attr(s_cleandata, 'scaled:center'),scale=attr(s_cleandata, 'scaled:scale'))
  
  att=rep(c(0,1),c(nin,nout))
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  SCaffected <- rowSums((dim(s_cleandata)[1]-1)*svds2[,sdsvds2 > delta0[2,],drop=F]^2)#
  ROCaffected <- simple_roc(labels=(att==1), scores=SCaffected)
  
  SCw <- rowSums((dim(s_cleandata)[1]-1)*(sweep(svds2,2,sdsvds2,FUN="*")^2))
  ROCw <- simple_roc(labels=(att==1), scores=SCw)
  
  
  XB <- NData%*%V[,en,drop=F]%*%t(V[,en,drop=F])
  SCWBPCA <- rowSums((NData-XB)^2)
  ROCWBPCA <- simple_roc(labels=(att==1), scores=SCWBPCA)
  
  
  a=data.frame(ROCaffected,ROCw,ROCWBPCA)
  sumo=sumo+a
}
sum0=sumo/1000

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=2,xlim=c(0,1))
lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)

#############################################

png(file="kdd99rocsmurf.png")
par(mar=c(6, 7, 4, 2))

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=1.5,xlim=c(0,.05),las=1)

lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(c)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=1.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=1.5, line=3.5)
legend("bottomright", legend=c("AAD", "WBPCA","WAAD"), 
       lty=c(1,2,3), col=c(1,1,1), 
       cex = 1.5,lwd=c(2,2,2))

dev.off()


#####################################################################
#####################################################################
#####################################################################

levels(KDDfullup$V42)
table(KDDfullup$V42)
KDDinf <- KDDfullup[KDDfullup$V42%in%c("neptune."),-c(which(cnms$V2==" symbolic."),42)]
colnames(KDDinf)=cnms$V1[cnms$V2==" continuous."]
KDDinfmod <- KDDinf[,-c(4,5,6,7,8,9,10,11,12,13,15)]


nout=100
nin=9900
n=nin+nout

res <- replicate(1000,{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=T),],
                       KDDinfmod[sample(1:dim(KDDinfmod)[1],nout,replace=T),]),  
                 center = attr(s_cleandata, 'scaled:center'),scale=attr(s_cleandata, 'scaled:scale')) 
  
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  sdsvds2 > delta0[2,]
})

barplot(apply(res,1,sum)/1000)

barplot(table(apply(res,2,sum)))

#######################################################
png(file="barplotmsneptune.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(table(apply(res,2,sum))/sum(table(apply(res,2,sum))),
                   ylim=c(0,max(table(apply(res,2,sum))/sum(table(apply(res,2,sum))))+.02),
                   col="gray",las=1,ylab="",cex.lab=1.5,cex.axis=1.5)
axis(1,barCenters, labels = FALSE)
title(ylab = "Proportion", cex.lab=1.5, line=4)
title(xlab="Number of components", cex.lab=1.5)
title(main="(a)",cex.main=2)
box()
#mtext("(a)",at=1, line = 1, cex =2)
dev.off()

png(file="barplotneptune.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(apply(res,1,sum)/1000*100,ylim=c(0,105),col="gray",
                   las=1,ylab="Percentage",cex.lab=1.5,cex.axis=1.5)
text(x = barCenters[seq(1,32,2)], y = par("usr")[3] - 3, srt = 45,
     adj = 1, labels = paste0("Comp-",seq(1,32,2)), xpd = TRUE,cex.lab=2.5)
axis(1,barCenters, labels = FALSE)
title(main="(b)",cex.main=2)
box()
#mtext("(b)",at=1, line = 1, cex =2)
dev.off()

# boxplot(res~row(res),outline=F)
# abline(h=1)

###################################################
nout=100
nin=9900
n=nin+nout
sumo=data.frame(matrix(0,n,9))

st <- 1
en <- eigen(cov(s_cleandata))$values>1



for (i in 1:1000) 
{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=T),],
                       KDDinfmod[sample(1:dim(KDDinfmod)[1],nout,replace=T),]),  
                 center = attr(s_cleandata, 'scaled:center'),scale=attr(s_cleandata, 'scaled:scale'))
  
  att=rep(c(0,1),c(nin,nout))
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  SCaffected <- rowSums((dim(s_cleandata)[1]-1)*svds2[,sdsvds2 > delta0[2,],drop=F]^2)#
  ROCaffected <- simple_roc(labels=(att==1), scores=SCaffected)
  
  SCw <- rowSums((dim(s_cleandata)[1]-1)*(sweep(svds2,2,sdsvds2,FUN="*")^2))
  ROCw <- simple_roc(labels=(att==1), scores=SCw)
  
  
  XB <- NData%*%V[,en,drop=F]%*%t(V[,en,drop=F])
  SCWBPCA <- rowSums((NData-XB)^2)
  ROCWBPCA <- simple_roc(labels=(att==1), scores=SCWBPCA)
  
  
  a=data.frame(ROCaffected,ROCw,ROCWBPCA)
  sumo=sumo+a
}
sum0=sumo/1000

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=2,xlim=c(0,1))
lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)

#############################################

png(file="kdd99rocneptune.png")
par(mar=c(6, 7, 4, 2))

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=1.5,xlim=c(0,.05),las=1)

lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(c)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=1.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=1.5, line=3.5)
legend("bottomright", legend=c("AAD", "WBPCA","WAAD"), 
       lty=c(1,2,3), col=c(1,1,1), 
       cex = 1.5,lwd=c(2,2,2))

dev.off()


#####################################################################
#####################################################################
#####################################################################

levels(KDDfullup$V42)
table(KDDfullup$V42)
KDDinf <- KDDfullup[KDDfullup$V42%in%c("satan."),-c(which(cnms$V2==" symbolic."),42)]
colnames(KDDinf)=cnms$V1[cnms$V2==" continuous."]
KDDinfmod <- KDDinf[,-c(4,5,6,7,8,9,10,11,12,13,15)]


nout=100
nin=9900
n=nin+nout

res <- replicate(1000,{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=T),],
                       KDDinfmod[sample(1:dim(KDDinfmod)[1],nout,replace=T),]),  
                 center = attr(s_cleandata, 'scaled:center'),scale=attr(s_cleandata, 'scaled:scale'))
  
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  sdsvds2 > delta0[2,]
})

barplot(apply(res,1,sum)/1000)

barplot(table(apply(res,2,sum)))

#######################################################
png(file="barplotmssatan.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(table(apply(res,2,sum))/sum(table(apply(res,2,sum))),
                   ylim=c(0,max(table(apply(res,2,sum))/sum(table(apply(res,2,sum))))+.02),
                   col="gray",las=1,ylab="",cex.lab=1.5,cex.axis=1.5)
axis(1,barCenters, labels = FALSE)
title(ylab = "Proportion", cex.lab=1.5, line=4)
title(xlab="Number of components", cex.lab=1.5)
title(main="(a)",cex.main=2)
box()
#mtext("(a)",at=1, line = 1, cex =2)
dev.off()

png(file="barplotsatan.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(apply(res,1,sum)/1000*100,ylim=c(0,105),col="gray",
                   las=1,ylab="Percentage",cex.lab=1.5,cex.axis=1.5)
text(x = barCenters[seq(1,32,2)], y = par("usr")[3] - 3, srt = 45,
     adj = 1, labels = paste0("Comp-",seq(1,32,2)), xpd = TRUE,cex.lab=2.5)
axis(1,barCenters, labels = FALSE)
title(main="(b)",cex.main=2)
box()
#mtext("(b)",at=1, line = 1, cex =2)
dev.off()

# boxplot(res~row(res),outline=F)
# abline(h=1)

###################################################
nout=100
nin=9900
n=nin+nout
sumo=data.frame(matrix(0,n,9))

st <- 1
en <- eigen(cov(s_cleandata))$values>1



for (i in 1:1000) 
{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=T),],
                       KDDinfmod[sample(1:dim(KDDinfmod)[1],nout,replace=T),]),  
                 center = attr(s_cleandata, 'scaled:center'),scale=attr(s_cleandata, 'scaled:scale'))
  
  att=rep(c(0,1),c(nin,nout))
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  SCaffected <- rowSums((dim(s_cleandata)[1]-1)*svds2[,sdsvds2 > delta0[2,],drop=F]^2)#
  ROCaffected <- simple_roc(labels=(att==1), scores=SCaffected)
  
  SCw <- rowSums((dim(s_cleandata)[1]-1)*(sweep(svds2,2,sdsvds2,FUN="*")^2))
  ROCw <- simple_roc(labels=(att==1), scores=SCw)
  
  
  XB <- NData%*%V[,en,drop=F]%*%t(V[,en,drop=F])
  SCWBPCA <- rowSums((NData-XB)^2)
  ROCWBPCA <- simple_roc(labels=(att==1), scores=SCWBPCA)
  
  
  a=data.frame(ROCaffected,ROCw,ROCWBPCA)
  sumo=sumo+a
}
sum0=sumo/1000

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=2,xlim=c(0,1))
lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)

#############################################

png(file="kdd99rocsatan.png")
par(mar=c(6, 7, 4, 2))

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=1.5,xlim=c(0,.05),las=1)

lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)
#mtext("(a)",at=0, line = 1, cex =2.5)
title(main="(c)",cex.main=2)
title(xlab="False Positive Rate", cex.lab=1.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=1.5, line=3.5)
legend("bottomright", legend=c("AAD", "WBPCA","WAAD"), 
       lty=c(1,2,3), col=c(1,1,1), 
       cex = 1.5,lwd=c(2,2,2))

dev.off()

