setwd("C:\\KDD")

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

library(data.table)
library(plyr)
library(expm)
library(pcaMethods)

ADFA_NB15_dat1 <- read.csv('https://www.unsw.adfa.edu.au/unsw-canberra-cyber/cybersecurity/ADFA-NB15-Datasets/UNSW-NB15_1.csv',header=FALSE)
ADFA_NB15_dat2 <- read.csv('https://www.unsw.adfa.edu.au/unsw-canberra-cyber/cybersecurity/ADFA-NB15-Datasets/UNSW-NB15_2.csv',header=FALSE)
ADFA_NB15_dat3 <- read.csv('https://www.unsw.adfa.edu.au/unsw-canberra-cyber/cybersecurity/ADFA-NB15-Datasets/UNSW-NB15_3.csv',header=FALSE)
ADFA_NB15_dat4 <- read.csv('https://www.unsw.adfa.edu.au/unsw-canberra-cyber/cybersecurity/ADFA-NB15-Datasets/UNSW-NB15_4.csv',header=FALSE)

var_names <- read.csv('https://www.unsw.adfa.edu.au/unsw-canberra-cyber/cybersecurity/ADFA-NB15-Datasets/NUSW-NB15_features.csv',header=TRUE)
  
colnames(ADFA_NB15_dat1)=var_names$Name
colnames(ADFA_NB15_dat2)=var_names$Name
colnames(ADFA_NB15_dat3)=var_names$Name
colnames(ADFA_NB15_dat4)=var_names$Name


vars <- c("dur","Spkts","Dpkts","sbytes","dbytes","sttl","dttl","Sload","Dload","sloss","dloss","Sintpkt",
          "Dintpkt","Sjit","Djit","swin","stcpb","dtcpb","dwin","tcprtt","synack","ackdat","smeansz","dmeansz",
          "trans_depth","res_bdy_len","Label","attack_cat")

ADFA_NB15_dat_q1 <- ADFA_NB15_dat1[,as.character(var_names$Name)[as.character(var_names$Name) %chin%  vars]]
ADFA_NB15_dat_q2 <- ADFA_NB15_dat2[,as.character(var_names$Name)[as.character(var_names$Name) %chin%  vars]]
ADFA_NB15_dat_q3 <- ADFA_NB15_dat3[,as.character(var_names$Name)[as.character(var_names$Name) %chin%  vars]]
ADFA_NB15_dat_q4 <- ADFA_NB15_dat4[,as.character(var_names$Name)[as.character(var_names$Name) %chin%  vars]]

# write.csv(ADFA_NB15_dat_q1, file = "ADFA_NB15_dat_q1.csv",row.names=FALSE, na="")
# write.csv(ADFA_NB15_dat_q2, file = "ADFA_NB15_dat_q2.csv",row.names=FALSE, na="")
# write.csv(ADFA_NB15_dat_q3, file = "ADFA_NB15_dat_q3.csv",row.names=FALSE, na="")
# write.csv(ADFA_NB15_dat_q4, file = "ADFA_NB15_dat_q4.csv",row.names=FALSE, na="")

ADFA_NB15_dat <- rbind(ADFA_NB15_dat_q1,ADFA_NB15_dat_q2,ADFA_NB15_dat_q3,ADFA_NB15_dat_q4)

ADFA_NB15_dat$attack_cat=revalue(ADFA_NB15_dat$attack_cat, c(" Fuzzers "="Fuzzers", 
                                                             " Fuzzers"="Fuzzers", 
                                                             " Reconnaissance "="Reconnaissance",
                                                             " Shellcode "="Shellcode",
                                                             "Backdoor"="Backdoors"))

# levels(ADFA_NB15_dat$attack_cat)
# table(ADFA_NB15_dat$Label)
# table(ADFA_NB15_dat$attack_cat)



#normalD <- which(ADFA_NB15_dat$attack_cat=="")
cleandata <- ADFA_NB15_dat[300001:900000,-c(27,28)]
#cleandata <- ADFA_NB15_dat[300001:900000,-c(8,20,21,27,28)]#200001:1000000
cleandata[is.na(cleandata)]=0# this is because the 0 in some variables 0 is coded as NA
s_cleandata <- scale(cleandata)

ss <- svd(s_cleandata)
U=ss$u
D=diag(ss$d)
V=ss$v


############################################
BTdelta <- replicate(5000,{
  XB <-  s_cleandata[sample(1:dim(s_cleandata)[1],30000,replace=TRUE),]
  UPB <- XB%*%V%*%solve(D)
  #max(diag(cov(UPB)*(dim(s_cleandata)[1]-1)))
  #max(eigen(sqrtm(cov(UPB))*sqrt(dim(s_cleandata)[1]-1))$values)
  apply(UPB*sqrt(dim(s_cleandata)[1]-1),2,sd)
})

delta0=apply(BTdelta,1,quantile,probs = c(0.01,0.99))
dim(BTdelta)
boxplot(BTdelta~row(BTdelta),outline=F)
abline(h=1)

png(file="Fig11.png")
par(mar=c(6, 6, 4, 2))
boxplot(BTdelta~row(BTdelta),outline=F, xlab="", ylab="", cex.axis=1.5, las=1, main="")
abline(h=1,col="gray",lwd=1)
title(xlab="Index", cex.lab=1.5, line = 3.5)
title(ylab=bquote(paste(s[j]^u)), cex.lab=1.5, line=3.5)
dev.off()

which(apply(BTdelta,1,quantile,0.75) < 1)
abs(V[,which(apply(BTdelta,1,quantile,0.75) < 1)]) > .3
V[,which(apply(BTdelta,1,quantile,0.75) < 1)]
names(cleandata)[abs(V[,which(apply(BTdelta,1,quantile,0.75) < 1)]) > .3]

hist(cleandata[,21])
#cleandata <- ADFA_NB15_dat[300001:900000,-c(27,28)]
out=c(which(cleandata[,20]>100000)+300000,
      which(cleandata[,21]>100000)+300000)

ADFA_NB15_dat <- rbind(ADFA_NB15_dat_q1,ADFA_NB15_dat_q2,ADFA_NB15_dat_q3,ADFA_NB15_dat_q4)

ADFA_NB15_dat <- ADFA_NB15_dat[-out,]
ADFA_NB15_dat$attack_cat=revalue(ADFA_NB15_dat$attack_cat, c(" Fuzzers "="Fuzzers", 
                                                             " Fuzzers"="Fuzzers", 
                                                             " Reconnaissance "="Reconnaissance",
                                                             " Shellcode "="Shellcode",
                                                             "Backdoor"="Backdoors"))


# levels(ADFA_NB15_dat$attack_cat)
# table(ADFA_NB15_dat$Label)
# table(ADFA_NB15_dat$attack_cat)



#normalD <- which(ADFA_NB15_dat$attack_cat=="")
cleandata <- ADFA_NB15_dat[300001:900000,-c(27,28)]
#cleandata <- ADFA_NB15_dat[300001:900000,-c(8,20,21,27,28)]#200001:1000000
cleandata[is.na(cleandata)]=0# this is because the 0 in some variables 0 is coded as NA
s_cleandata <- scale(cleandata)

ss <- svd(s_cleandata)
U=ss$u
D=diag(ss$d)
V=ss$v

ERall <- rowSums(U^2)
plot(log(ERall),pch=".")
# rss <- robustSvd(s_cleandata)
############################################
BTdelta <- replicate(5000,{
  XB <-  s_cleandata[sample(1:dim(s_cleandata)[1],30000,replace=TRUE),]
  UPB <- XB%*%V%*%solve(D)
  #max(diag(cov(UPB)*(dim(s_cleandata)[1]-1)))
  #max(eigen(sqrtm(cov(UPB))*sqrt(dim(s_cleandata)[1]-1))$values)
  apply(UPB*sqrt(dim(s_cleandata)[1]-1),2,sd)
})

delta0=apply(BTdelta,1,quantile,probs = c(0.01,0.99))
dim(BTdelta)
boxplot(BTdelta~row(BTdelta),outline=F)
abline(h=1)

png(file="Fig10.png")
par(mar=c(6, 6, 4, 2))
boxplot(BTdelta~row(BTdelta),outline=F, xlab="", ylab="", cex.axis=1.5, las=1, main="")
abline(h=1,col="gray",lwd=1)
title(xlab="Index", cex.lab=1.5, line = 3.5)
title(ylab=bquote(paste(s[j]^u)), cex.lab=1.5, line=3.5)
dev.off()

#############################################


# attackdata <- rbind(ADFA_NB15_dat[ADFA_NB15_dat$attack_cat%in%names(table(ADFA_NB15_dat$attack_cat))[-1],-c(27,28)])
# attlabels <- ADFA_NB15_dat$Label#rep(c(1,0),c(321283,300000))
# 
# attackdata <- rbind(ADFA_NB15_dat[ADFA_NB15_dat$attack_cat=="Analysis",-c(27,28)],
#                     cleandata[sample(1:dim(cleandata)[1],50000),])

attackdata <- ADFA_NB15_dat[,-c(27,28)]
attlabels <- ADFA_NB15_dat$Label
# attackdata <- ADFA_NB15_dat[ADFA_NB15_dat$attack_cat=="Fuzzers",-c(27,28)]
NData <- scale(attackdata, center = attr(s_cleandata, 'scaled:center'),
               scale=attr(s_cleandata, 'scaled:scale'))
dim(NData)

svds <- NData%*%ss$v%*%diag(1/ss$d)
sdsvds <- apply(svds*sqrt(dim(s_cleandata)[1]-1),2,sd)

plot(sdsvds)
abline(h=1)

############################################
############################################
BTw <- replicate(1000,{
  XB <-  s_cleandata[sample(1:dim(s_cleandata)[1],10000,replace=TRUE),]
  UPB <- XB%*%V%*%diag(1/ss$d)
  sc <- log(rowSums((dim(s_cleandata)[1]-1)*(sweep(UPB,2,sdsvds,FUN="*")^2)))
})
thetaw=quantile(as.vector(BTw),0.9999)
############################################

ERw <- rowSums((dim(s_cleandata)[1]-1)*(sweep(svds,2,sdsvds,FUN="*")^2))
# plot(log(ERw),col=(attlabels+1),ylim=c(0,max(log(ERw))),
#      pch=ifelse(attlabels==0,".","."))
# 
# abline(h=thetaw)
table(log(ERw)>thetaw,attlabels)

100-(sum(attlabels==1)-sum(log(ERw)[attlabels==1]>thetaw))/sum(attlabels==1)*100

############################################
############################################
BTall <- replicate(1000,{
  XB <-  s_cleandata[sample(1:dim(s_cleandata)[1],10000,replace=TRUE),]
  UPB <- XB%*%V%*%diag(1/ss$d)
  sc <- log(rowSums((dim(s_cleandata)[1]-1)*UPB^2))
})
thetaall=quantile(as.vector(BTall),0.9999)

ERall <- rowSums((dim(s_cleandata)[1]-1)*svds^2)
# plot(log(ERall),col=attlabels+1,
#      pch=ifelse(attlabels==0,".","."))
# abline(h=thetaall)
table(log(ERall)>thetaall,attlabels)

100-(sum(attlabels==1)-sum(log(ERall)[attlabels==1]>thetaall))/sum(attlabels==1)*100
############################################
############################################
ind <- which(sdsvds > delta0[2,])
BTaffected <- replicate(1000,{
  XB <-  s_cleandata[sample(1:dim(s_cleandata)[1],10000,replace=TRUE),]
  UPB <- XB%*%V%*%diag(1/ss$d)
  sc <- log(rowSums((dim(s_cleandata)[1]-1)*UPB[, ind, drop=F]^2))
})
thetaaffected=quantile(as.vector(BTaffected),0.9999)
############################################
ERaffected <- rowSums((dim(s_cleandata)[1]-1)*svds[, ind, drop=F]^2)
# plot(log(ERaffected),col=attlabels+1,
#      pch=ifelse(attlabels==0,".","."))
# abline(h=thetaaffected)
table(log(ERaffected)>thetaaffected,attlabels)

100-(sum(attlabels==1)-sum(log(ERaffected)[attlabels==1]>thetaaffected))/sum(attlabels==1)*100
############################################
############################################
ind <- which(sdsvds > 30)
BTmaffected <- replicate(1000,{
  XB <-  s_cleandata[sample(1:dim(s_cleandata)[1],10000,replace=TRUE),]
  UPB <- XB%*%V%*%diag(1/ss$d)
  sc <- log(rowSums((dim(s_cleandata)[1]-1)*UPB[, ind, drop=F]^2))
})
thetamaffected=quantile(as.vector(BTmaffected),0.9999)
############################################
ERmaffected <- rowSums((dim(s_cleandata)[1]-1)*svds[, ind, drop=F]^2)
# plot(log(ERmaffected),col=attlabels+1,
#      pch=ifelse(attlabels==0,".","."))
# abline(h=thetamaffected)
table(log(ERmaffected)>thetamaffected,attlabels)

100-(sum(attlabels==1)-sum(log(ERmaffected)[attlabels==1]>thetamaffected))/sum(attlabels==1)*100
############################################
############################################

st <- 1
en <- sum(eigen(cov(s_cleandata))$values>1)
ERT <- rowSums((s_cleandata-s_cleandata%*%V[,st:en,drop=F]%*%t(V[,st:en,drop=F]))^2)

BTWBPACA <- replicate(1000,{
  B <-  quantile(sample(ERT[1:dim(s_cleandata)[1]],10000,replace=TRUE),0.9999)
})
thetaWBPACA <- log(mean(BTWBPACA))

XB <- NData%*%V[,st:en,drop=F]%*%t(V[,st:en,drop=F])
ERWBPACA <- rowSums((NData-XB)^2)

# plot(log(ERWBPACA),col=(attlabels+1),ylim=c(0,max(log(ERWBPACA))),
#      pch=ifelse(attlabels==0,".","."))
# 
# abline(h=thetaWBPACA)

table(log(ERWBPACA)>thetaWBPACA,attlabels)

100-(sum(attlabels==1)-sum(log(ERWBPACA)[attlabels==1]>thetaWBPACA))/sum(attlabels==1)*100

################################################################
################################################################
library(zoo)
png(file="Fig4.png")
par(mar=c(6, 6, 4, 2))
x=index(ERaffected)
plot(x, log(ERaffected),  pch=ifelse(ADFA_NB15_dat$Label==0,".",".") ,col=(ADFA_NB15_dat$Label+1) ,
     cex=2,xlab="", ylab="", cex.axis=1.5, ylim=c(min(log(ERaffected)),max(log(ERaffected))+2))
rect(300001, min(log(ERaffected))-5, 900000, max(log(ERaffected))+5, border = rgb(red = 0, green = 0, blue = 1, alpha = 0.1), 
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.05))
# rect(186789, min(log(ERaffected))-5, 500000, max(log(ERaffected))+5, border = rgb(red = 0, green = 1, blue = 0, alpha = 0.1), 
#      col = rgb(red = 0, green = 1, blue = 0, alpha = 0.05))
#axis(1, at=c(1,5,10,15,20,25)*1e+05, labels=c(1,2,3,4,5,6),cex.axis=2)
abline(h=thetaaffected,col="green",cex=2)
mtext("(a)",at=0, line = 1, cex =2.5)
#title(main="(a)",cex.main=2)
title(xlab="Time index", cex.lab=1.5, line = 3.5)
title(ylab=expression(log(t[i])), cex.lab=1.5, line=3.5)
legend("topleft", legend=c("Normal", "Infected"), 
       pch=c(".","."), col=c(1,2), text.col = c(1,2),
       cex = 1)

dev.off()

################################################################

png(file="Fig5.png")
par(mar=c(6, 6, 4, 2))
x=index(ERmaffected)
plot(x, log(ERmaffected),  pch=ifelse(ADFA_NB15_dat$Label==0,".",".") ,col=(ADFA_NB15_dat$Label+1) ,
     cex=2,xlab="", ylab="", cex.axis=1.5, ylim=c(min(log(ERmaffected)),max(log(ERmaffected))+2))
rect(300001, min(log(ERmaffected))-5, 900000, max(log(ERmaffected))+5, border = rgb(red = 0, green = 0, blue = 1, alpha = 0.1), 
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.05))
# rect(186789, min(log(ERaffected))-5, 500000, max(log(ERaffected))+5, border = rgb(red = 0, green = 1, blue = 0, alpha = 0.1), 
#      col = rgb(red = 0, green = 1, blue = 0, alpha = 0.05))
#axis(1, at=c(1,5,10,15,20,25)*1e+05, labels=c(1,2,3,4,5,6),cex.axis=2)
abline(h=thetamaffected,col="green",cex=2)
mtext("(b)",at=0, line = 1, cex =2.5)
#title(main="(a)",cex.main=2)
title(xlab="Time index", cex.lab=1.5, line = 3.5)
title(ylab=expression(log(t[i])), cex.lab=1.5, line=3.5)
legend("topleft", legend=c("Normal", "Infected"), 
       pch=c(".","."), col=c(1,2), text.col = c(1,2),
       cex = 1)

dev.off()

################################################################

png(file="Fig6.png")
par(mar=c(6, 6, 4, 2))
x=index(ERw)
plot(x, log(ERw),  pch=ifelse(ADFA_NB15_dat$Label==0,".",".") ,col=(ADFA_NB15_dat$Label+1) ,
     cex=2,xlab="", ylab="", cex.axis=1.5, ylim=c(min(log(ERw)),max(log(ERw))+2))
rect(300001, min(log(ERw))-5, 900000, max(log(ERw))+5, border = rgb(red = 0, green = 0, blue = 1, alpha = 0.1), 
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.05))
# rect(186789, min(log(ERw))-5, 500000, max(log(ERw))+5, border = rgb(red = 0, green = 1, blue = 0, alpha = 0.1), 
#      col = rgb(red = 0, green = 1, blue = 0, alpha = 0.05))
#axis(1, at=c(1,5,10,15,20,25)*1e+05, labels=c(1,2,3,4,5,6),cex.axis=2)
abline(h=thetaw,col="green",cex=2)
mtext("(c)",at=0, line = 1, cex =2.5)
#title(main="(a)",cex.main=2)
title(xlab="Time index", cex.lab=1.5, line = 3.5)
title(ylab=expression(log(t[i])), cex.lab=1.5, line=3.5)
legend("topleft", legend=c("Normal", "Infected"), 
       pch=c(".","."), col=c(1,2), text.col = c(1,2),
       cex = 1)

dev.off()

#############################################################
png(file="Fig7.png")
par(mar=c(6, 6, 4, 2))
x=index(ERWBPACA)
plot(x, log(ERWBPACA),  pch=ifelse(ADFA_NB15_dat$Label==0,".",".") ,col=(ADFA_NB15_dat$Label+1) ,
     cex=2,xlab="", ylab="", cex.axis=1.5, ylim=c(min(log(ERWBPACA)),max(log(ERWBPACA))+2))
rect(300001, min(log(ERWBPACA))-5, 900000, max(log(ERWBPACA))+5, border = rgb(red = 0, green = 0, blue = 1, alpha = 0.1), 
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.05))
# rect(186789, min(log(ERWBPACA))-5, 500000, max(log(ERWBPACA))+5, border = rgb(red = 0, green = 1, blue = 0, alpha = 0.1), 
#      col = rgb(red = 0, green = 1, blue = 0, alpha = 0.05))

#axis(1, at=c(1,5,10,15,20,25)*1e+05, labels=c(1,2,3,4,5,6),cex.axis=2)
abline(h=thetaWBPACA,col="green",cex=2)
mtext("(d)",at=0, line = 1, cex =2.5)
#title(main="(a)",cex.main=2)
title(xlab="Time index", cex.lab=1.5, line = 3.5)
title(ylab=expression(log(t[i])), cex.lab=1.5, line=3.5)
legend("topleft", legend=c("Normal", "Infected"), 
       pch=c(".","."), col=c(1,2), text.col = c(1,2),
       cex = 1)

dev.off()



###############################################################################################
ROCw <- simple_roc(labels=(attlabels==1), scores=ERw)
ROCaffected <- simple_roc(labels=(attlabels==1), scores=ERaffected)
ROCWBAPCA <- simple_roc(labels=(attlabels==1), scores=ERWBPACA)

png(file="Fig8.png")
par(mar=c(6, 7, 4, 2))
plot(ROCw$FPR,ROCw$TPR,type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=1.5,xlim=c(0,.1),lty=3, las=1)
lines(ROCaffected$FPR,ROCaffected$TPR,type="l",col=1,lwd=2,lty=1)
lines(ROCWBAPCA$FPR,ROCWBAPCA$TPR,type="l",col=1,lwd=2,lty=2)
title(xlab="False Positive Rate", cex.lab=1.5, line = 3.5)
title(ylab="True Positive Rate", cex.lab=1.5, line=3.5)
legend("bottomright", legend=c("AAD", "WBPCA", "WAAD"), 
       lty=c(1,2,3), col=c(1,1,1), text.col = c(1,1,1),
       cex = 1)

dev.off()

###############################################################################################
png(file="Fig9.png")
par(mar=c(6, 6, 4, 2))
plot(apply(svds*sqrt(dim(s_cleandata)[1]-1),2,sd),pch=16, col=1,lwd=2,xlab="", ylab="", cex.axis=1.5,lty=3, las=1)
#abline(h=1,col="gray",lwd=1)
title(xlab="Index", cex.lab=1.5, line = 3.5)
title(ylab=expression(s[j]^u), cex.lab=1.5, line=3.5)
legend("topleft", legend=c("AAD", "WBPCA", "WAAD"), 
       lty=c(1,2,3), col=c(1,1,1), text.col = c(1,1,1), cex = 1)
dev.off()

################################################################

# png(file="Fig10.png")
# par(mar=c(6, 6, 4, 2))
# hist(BTdelta,freq=F,xlim = c(1,7),ylim=c(0,.6),xlab="", ylab="", cex.axis=1.5, las=1, main="")
# abline(v=delta0,col="gray",lwd=2)
# title(xlab=bquote(paste("max(",r[j]^u,")")), cex.lab=1.5, line = 3.5)
# title(ylab="Density", cex.lab=1.5, line=3.5)
# dev.off()
###############################################################################################
# bquote(paste("(",n[1],", ",n[2], ")", " = ","(",.(n1),", ",.(n2),")",", ",rho," = ","1/4" ))
# T1 <- table(log(ERw)>thetaw,attlabels)
# T2 <- table(log(ERaffected)>thetaaffected,attlabels)
# T3 <- table(log(ERWBPACA)>log(thetaWBPACA),attlabels)
# 
# c(c(T1[2,2]/(T1[2,1]+T1[2,2]),T2[2,2]/(T2[2,1]+T2[2,2]),T3[2,2]/(T3[2,1]+T3[2,2])),
#   c(T1[1,2]/(T1[1,1]+T1[1,2]),T2[1,2]/(T2[1,1]+T2[1,2]),T3[1,2]/(T3[1,1]+T3[1,2])))


#############################################

table(ADFA_NB15_dat$attack_cat)
levels(ADFA_NB15_dat$attack_cat)

attackdata <- ADFA_NB15_dat[!ADFA_NB15_dat$attack_cat=="",-c(27,28)]
# attackdata[is.na(attackdata)]=0
# attackdata <- ADFA_NB15_dat[-c(200001:1000000),-c(8,20,21,27,28)]

nout=100
nin=9900
n=nin+nout

res <- replicate(1000,{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=F),],
                       attackdata[sample(1:dim(attackdata)[1],nout,replace=F),]), 
                 center = attr(s_cleandata, 'scaled:center'),
                 scale=attr(s_cleandata, 'scaled:scale'))
  
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  sdsvds2 > delta0[2,]
})

# boxplot(res~row(res),outline=F)
# abline(h=1)
barplot(apply(res,1,sum)/1000)

barplot(table(apply(res,2,sum)))

#######################################################
png(file="barplotncUNSWall.png")
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

png(file="barplotwhichcUNSWall.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(apply(res,1,sum)/1000*100,ylim=c(0,105),col="gray",
                   las=1,ylab="Percentage",cex.lab=1.5,cex.axis=1.5)
text(x = barCenters[seq(1,26,2)], y = par("usr")[3] - 3, srt = 45,
     adj = 1, labels = paste0("Comp-",seq(1,26,2)), xpd = TRUE,cex.lab=2.5)
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
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=F),],
                       attackdata[sample(1:dim(attackdata)[1],nout,replace=F),]), 
                 center = attr(s_cleandata, 'scaled:center'),
                 scale=attr(s_cleandata, 'scaled:scale'))
  
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

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=2,xlim=c(0,.01))
lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)

#############################################

png(file="kdd99rocUNSWall.png")
par(mar=c(6, 7, 4, 2))

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=1.5,xlim=c(0,.0001),las=1)

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

#############################################

table(ADFA_NB15_dat$attack_cat)
levels(ADFA_NB15_dat$attack_cat)

attackdata <- ADFA_NB15_dat[ADFA_NB15_dat$attack_cat=="Generic",-c(27,28)]
# attackdata[is.na(attackdata)]=0
# attackdata <- ADFA_NB15_dat[-c(200001:1000000),-c(8,20,21,27,28)]

nout=100
nin=9900
n=nin+nout

res <- replicate(1000,{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=F),],
                       attackdata[sample(1:dim(attackdata)[1],nout,replace=F),]), 
                 center = attr(s_cleandata, 'scaled:center'),
                 scale=attr(s_cleandata, 'scaled:scale'))
  
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  sdsvds2 > delta0[2,]
})

# boxplot(res~row(res),outline=F)
# abline(h=1)
barplot(apply(res,1,sum)/1000)

barplot(table(apply(res,2,sum)))

#######################################################
png(file="barplotncUNSWGeneric.png")
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

png(file="barplotwhichcUNSWGeneric.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(apply(res,1,sum)/1000*100,ylim=c(0,105),col="gray",
                   las=1,ylab="Percentage",cex.lab=1.5,cex.axis=1.5)
text(x = barCenters[seq(1,26,2)], y = par("usr")[3] - 3, srt = 45,
     adj = 1, labels = paste0("Comp-",seq(1,26,2)), xpd = TRUE,cex.lab=2.5)
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
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=F),],
                       attackdata[sample(1:dim(attackdata)[1],nout,replace=F),]), 
                 center = attr(s_cleandata, 'scaled:center'),
                 scale=attr(s_cleandata, 'scaled:scale'))
  
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

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=2,xlim=c(0,.01))
lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)

#############################################

png(file="kdd99rocUNSWGeneric.png")
par(mar=c(6, 7, 4, 2))

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=1.5,xlim=c(0,.0001),las=1)

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

#############################################
#############################################
#############################################

table(ADFA_NB15_dat$attack_cat)
levels(ADFA_NB15_dat$attack_cat)

attackdata <- ADFA_NB15_dat[ADFA_NB15_dat$attack_cat=="Fuzzers",-c(27,28)]
# attackdata[is.na(attackdata)]=0
# attackdata <- ADFA_NB15_dat[-c(200001:1000000),-c(8,20,21,27,28)]

nout=100
nin=9900
n=nin+nout

res <- replicate(1000,{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=F),],
                       attackdata[sample(1:dim(attackdata)[1],nout,replace=F),]), 
                 center = attr(s_cleandata, 'scaled:center'),
                 scale=attr(s_cleandata, 'scaled:scale'))
  
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  sdsvds2 > delta0[2,]
})

# boxplot(res~row(res),outline=F)
# abline(h=1)
barplot(apply(res,1,sum)/1000)

barplot(table(apply(res,2,sum)))

#######################################################
png(file="barplotncUNSWFuzzers.png")
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

png(file="barplotwhichcUNSWFuzzers.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(apply(res,1,sum)/1000*100,ylim=c(0,105),col="gray",
                   las=1,ylab="Percentage",cex.lab=1.5,cex.axis=1.5)
text(x = barCenters[seq(1,26,2)], y = par("usr")[3] - 3, srt = 45,
     adj = 1, labels = paste0("Comp-",seq(1,26,2)), xpd = TRUE,cex.lab=2.5)
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
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=F),],
                       attackdata[sample(1:dim(attackdata)[1],nout,replace=F),]), 
                 center = attr(s_cleandata, 'scaled:center'),
                 scale=attr(s_cleandata, 'scaled:scale'))
  
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

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=2,xlim=c(0,.01))
lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)

#############################################

png(file="kdd99rocUNSWFuzzers.png")
par(mar=c(6, 7, 4, 2))

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=1.5,xlim=c(0,.0001),las=1)

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


#############################################
#############################################
#############################################

table(ADFA_NB15_dat$attack_cat)
levels(ADFA_NB15_dat$attack_cat)

attackdata <- ADFA_NB15_dat[ADFA_NB15_dat$attack_cat=="DoS",-c(27,28)]
# attackdata[is.na(attackdata)]=0
# attackdata <- ADFA_NB15_dat[-c(200001:1000000),-c(8,20,21,27,28)]

nout=100
nin=9900
n=nin+nout

res <- replicate(1000,{
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=F),],
                       attackdata[sample(1:dim(attackdata)[1],nout,replace=F),]), 
                 center = attr(s_cleandata, 'scaled:center'),
                 scale=attr(s_cleandata, 'scaled:scale'))
  
  svds2 <- NData%*%V%*%solve(D)
  sdsvds2 <- apply(svds2*sqrt(dim(s_cleandata)[1]-1),2,sd)
  
  sdsvds2 > delta0[2,]
})

# boxplot(res~row(res),outline=F)
# abline(h=1)
barplot(apply(res,1,sum)/1000)

barplot(table(apply(res,2,sum)))

#######################################################
png(file="barplotncUNSWDoS.png")
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

png(file="barplotwhichcUNSWDoS.png")
par(mar=c(6, 7, 4, 2))
barCenters=barplot(apply(res,1,sum)/1000*100,ylim=c(0,105),col="gray",
                   las=1,ylab="Percentage",cex.lab=1.5,cex.axis=1.5)
text(x = barCenters[seq(1,26,2)], y = par("usr")[3] - 3, srt = 45,
     adj = 1, labels = paste0("Comp-",seq(1,26,2)), xpd = TRUE,cex.lab=2.5)
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
  
  NData <- scale(rbind(cleandata[sample(1:dim(cleandata)[1],nin,replace=F),],
                       attackdata[sample(1:dim(attackdata)[1],nout,replace=F),]), 
                 center = attr(s_cleandata, 'scaled:center'),
                 scale=attr(s_cleandata, 'scaled:scale'))
  
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

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=2,xlim=c(0,.0001))
lines(sum0[,5],sum0[,4],type="l",col=1,lwd=2,lty=3)
lines(sum0[,8],sum0[,7],type="l",col=1,lwd=2,lty=2)

#############################################

png(file="kdd99rocUNSWDoS.png")
par(mar=c(6, 7, 4, 2))

plot(sum0[,2],sum0[,1],type="l",col=1,lwd=2,xlab="", ylab="", cex.axis=1.5,xlim=c(0,.0001),las=1)

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


