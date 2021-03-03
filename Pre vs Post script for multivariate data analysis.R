#remove peaks not present in at least 80% of samples in each group
setwd("D:/MEGA/Publication 2")
library(readr)
UF_PILOT4_3_withexcl <- read_csv("ENP_final_3_new_FINAL_peaks_excluded_new_PVP_fm.csv")
UF_PILOT4_3_withexcl$`2.022404917`=NULL
#UF_PILOT4_3_withexcl$`1.457467423`=NULL
#UF_PILOT4_3_withexcl$`1.471949712`=NULL
#UF_PILOT4_3_withexcl$`1.872658199`=NULL
UF_PILOT4_3_withexcl$`1.903997967`=NULL
#UF_PILOT4_3_withexcl$`2.144890749`=NULL
UF_PILOT4_3_withexcl$`2.214349772`=NULL
#UF_PILOT4_3_withexcl$`2.253153549`=NULL
#UF_PILOT4_3_withexcl$`2.324851518`=NULL
#UF_PILOT4_3_withexcl$`3.126405303`=NULL
#UF_PILOT4_3_withexcl$`3.136864883`=NULL
UF_PILOT4_3_withexcl$`3.185692061`=NULL
#UF_PILOT4_3_withexcl$`3.262367198`=NULL
#UF_PILOT4_3_withexcl$`3.321208614`=NULL
#UF_PILOT4_3_withexcl$`3.554511486`=NULL
#UF_PILOT4_3_withexcl$`3.567271011`=NULL
#UF_PILOT4_3_withexcl$`3.612111183`=NULL
#UF_PILOT4_3_withexcl$`3.939434677`=NULL
#UF_PILOT4_3_withexcl$`3.95085602`=NULL
#UF_PILOT4_3_withexcl$`3.962700713`=NULL
#UF_PILOT4_3_withexcl$`3.972674294`=NULL
#UF_PILOT4_3_withexcl$`4.107654959`=NULL
jack=UF_PILOT4_3_withexcl
require("plyr")
count(jack, c("GROUP"))
jack=jack[,-which(colSums(is.na(jack[which(jack$GROUP=="QC"),]))>=26  | colSums(is.na(jack[which(jack$GROUP=="HV"),]))>=19)]
pearl=jack#[,-which(colSums(is.na(jack[which(jack$GROUP=="POST"),]))>=7 & colSums(is.na(jack[which(jack$GROUP=="PPGL"),]))>=7)]
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN.csv", pearl)

#scale to MA peak (5.99)
captain=pearl
captain[is.na(captain)]=0
captain=captain[,order(colnames(captain))]
captain=cbind(captain$`ENSAT-HT ID`, captain$GROUP, captain[,1:(as.numeric(ncol(captain))-2)])
sparrow=cbind(captain[,1:2], captain[,3:(as.numeric(ncol(captain)))]/captain$`5.996102403`)
sparrow$`5.996102403`<- NULL
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN_MA.csv", sparrow)

#PQN
james=sparrow[,3:ncol(sparrow)]
james[james==0]=NA
black=sparrow[sparrow$`captain$GROUP`=="HV",3:ncol(sparrow)]
dauntless=apply(black, 2, median, na.rm=T)
check=t(t(james)/dauntless)
check2=apply(check, 1, median, na.rm=T)
checkmate=james/check2
checkmate[is.na(checkmate)]=0
beckett=cbind(sparrow[,1:2], checkmate)
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN_MA_PQN.csv", beckett)

#remove peaks with a RSD>0.3 in replicate samples (in this case both PA, PPGL)
cutler=beckett
library(sjstats)
hay=function(x) {
  cv(cutler[which(cutler$`captain$GROUP`=="QC"),3:(as.numeric(ncol(cutler)))][,x])
}
cvQC=sapply((1:(as.numeric(ncol(cutler))-2)), hay)
norrington=checkmate[,which(cvQC<0.3)]
commodore=cbind(cutler[,1:2], norrington)
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN_MA_PQN_0.3.csv", commodore)
#before exclusion: "2.030244923" "2.102557302" "2.13235951"  "2.139543856" "2.443586856" "2.448745467" "2.62883786"  "2.651707926" "2.705949063" "3.679241009" "4.102385996" "4.115225666"
#after exclusion: "2.030244923" "2.102557302" "2.13235951"  "2.139543856" "2.443586856" "2.448745467" "2.62883786"  "2.651707926" "2.705949063" "3.679241009" "4.102385996" "4.115225666"

#kNN missing value estimation
library(impute)
ncheck=norrington
ncheck[ncheck==0]=NA
ncheckO.imputed <- impute.knn(as.matrix(ncheck), k=10, rowmax = 1, colmax = 1, maxp=nrow(ncheck))
commodore=cbind(commodore[,1:2], ncheckO.imputed$data)
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN_MA_PQN_0.3_knn10.csv", commodore)

#glog
library("LMGene")
library(Biobase)
library(tools)
library(readr)

ENP_vanilla_forGLOG=commodore
PLASMA_QC_PQN=ENP_vanilla_forGLOG[which(ENP_vanilla_forGLOG$`captain$GROUP`=="QC"),]
#PLASMA_QC_PQN=PLASMA_QC_PQN[-c(92:105, 148:166),]
QCled_PQN=as.matrix(t(PLASMA_QC_PQN))
QCled_PQN=QCled_PQN[-2,]
colnames(QCled_PQN)=QCled_PQN[1,]
QCled_PQN=QCled_PQN[-1,]
QCr=apply(QCled_PQN, 1,as.numeric)
QCr=t(QCr)
QCmonster=as.factor(c(rep(1:2, each=(as.numeric(nrow(PLASMA_QC_PQN))/2)),1))
QCdose=as.numeric(c(rep(1, times=as.numeric(nrow(PLASMA_QC_PQN)))))
QCled_list=list("monster"=QCmonster, "dose"=QCdose)
QCled.eS=neweS(QCr, QCled_list)
tranpar <- tranest(QCled.eS)
tranpar

PLASMA_PQN=ENP_vanilla_forGLOG
led_PQN=as.matrix(t(PLASMA_PQN))
led_PQN=led_PQN[-2,]
colnames(led_PQN)=led_PQN[1,]
led_PQN=led_PQN[-1,]
r=apply(led_PQN, 1,as.numeric)
r=t(r)
monster=as.factor(c(1:as.numeric(nrow(PLASMA_PQN))))
dose=as.numeric(c(rep(1, times=as.numeric(nrow(PLASMA_PQN)))))
led_list=list("monster"=monster, "dose"=dose)
led.eS=neweS(r, led_list)
trled.eS <- transeS(led.eS, tranpar$lambda, tranpar$alpha)
kostakis=exprs(trled.eS)
colnames(kostakis)=as.factor(colnames(kostakis))
colnames(kostakis)=colnames(led_PQN)
final=t(kostakis)
final=cbind(cutler[,1:2], final[,1:as.numeric(ncol(final))])
final=final[,-1]
#write.csv(final, file="ENP_final_3_bothQCexp_16072019_excl_final_initial_sample_exclusions_20_MA_PQN_0.3_knn10_GLOG.csv")

#robust PCA
library(rospca)
X=final[which(final$`captain$GROUP`=='QC'),-1]
res=robpca(X, k=0, skew = F)
diagPlot((res))
which(res$flag.all==F)
EXCL_QC=which(res$flag.all==F)

library(rospca)
X=final[which(final$`captain$GROUP`=='HV'),-1]
res=robpca(X, k=0, skew = F)
diagPlot((res))
which(res$flag.all==F)
EXCL_HV=which(res$flag.all==F)

EXCL=c(EXCL_QC, EXCL_HV)


ix <- which(UF_PILOT4_3_withexcl$`ENSAT-HT ID` %in% c(names(EXCL)))
jack=UF_PILOT4_3_withexcl[-ix,]
require("plyr")
count(jack, c("GROUP"))
jack=jack[,-which(colSums(is.na(jack[which(jack$GROUP=="QC"),]))>=22  | colSums(is.na(jack[which(jack$GROUP=="HV"),]))>=14)]
pearl=jack#[,-which(colSums(is.na(jack[which(jack$GROUP=="POST"),]))>=7 & colSums(is.na(jack[which(jack$GROUP=="PPGL"),]))>=7)]
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN.csv", pearl)

#scale to MA peak (5.99)
captain=pearl
captain[is.na(captain)]=0
captain=captain[,order(colnames(captain))]
captain=cbind(captain$`ENSAT-HT ID`, captain$GROUP, captain[,1:(as.numeric(ncol(captain))-2)])
sparrow=cbind(captain[,1:2], captain[,3:(as.numeric(ncol(captain)))]/captain$`5.996102403`)
sparrow$`5.996102403`<- NULL

#PQN
james=sparrow[,3:ncol(sparrow)]
james[james==0]=NA
black=sparrow[sparrow$`captain$GROUP`=="HV",3:ncol(sparrow)]
dauntless=apply(black, 2, median, na.rm=T)
check=t(t(james)/dauntless)
check2=apply(check, 1, median, na.rm=T)
checkmate=james/check2
checkmate[is.na(checkmate)]=0
beckett=cbind(sparrow[,1:2], checkmate)

#remove peaks with a RSD>0.3 in replicate samples (in this case both PA, PPGL)
cutler=beckett
library(sjstats)
hay=function(x) {
  cv(cutler[which(cutler$`captain$GROUP`=="QC"),3:(as.numeric(ncol(cutler)))][,x])
}
cvQC=sapply((1:(as.numeric(ncol(cutler))-2)), hay)
norrington=checkmate[,which(cvQC<0.3)]
commodore=cbind(cutler[,1:2], norrington)
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN_MA_PQN_0.3.csv", commodore)
#before exclusion: "2.030244923" "2.102557302" "2.13235951"  "2.139543856" "2.443586856" "2.448745467" "2.62883786"  "2.651707926" "2.705949063" "3.679241009" "4.102385996" "4.115225666"
#after exclusion: "1.860408043" "1.884843241" "1.9553774"   "2.030244923" "2.102557302" "2.13235951"  "2.139543856" "2.443586856" "2.448745467" "2.62883786"  "2.651707926" "2.705949063" "3.679241009" "3.995651938" "4.102385996" "4.115225666" "5.171586437"

#kNN missing value estimation
library(impute)
ncheck=norrington
ncheck[ncheck==0]=NA
ncheckO.imputed <- impute.knn(as.matrix(ncheck), k=10, rowmax = 1, colmax = 1, maxp=nrow(ncheck))
commodore=cbind(commodore[,1:2], ncheckO.imputed$data)
#write.csv(file="ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_20_new_nokb_outliers_NOCEN_MA_PQN_0.3_knn10.csv", commodore)

#glog
library("LMGene")
library(Biobase)
library(tools)
library(readr)

ENP_vanilla_forGLOG=commodore
PLASMA_QC_PQN=ENP_vanilla_forGLOG[which(ENP_vanilla_forGLOG$`captain$GROUP`=="QC"),]
#PLASMA_QC_PQN=PLASMA_QC_PQN[-c(92:105, 148:166),]
QCled_PQN=as.matrix(t(PLASMA_QC_PQN))
QCled_PQN=QCled_PQN[-2,]
colnames(QCled_PQN)=QCled_PQN[1,]
QCled_PQN=QCled_PQN[-1,]
QCr=apply(QCled_PQN, 1,as.numeric)
QCr=t(QCr)
QCmonster=as.factor(c(rep(1:2, each=(as.numeric(nrow(PLASMA_QC_PQN))/2))))
QCdose=as.numeric(c(rep(1, times=as.numeric(nrow(PLASMA_QC_PQN)))))
QCled_list=list("monster"=QCmonster, "dose"=QCdose)
QCled.eS=neweS(QCr, QCled_list)
tranpar <- tranest(QCled.eS)
tranpar

PLASMA_PQN=ENP_vanilla_forGLOG
led_PQN=as.matrix(t(PLASMA_PQN))
led_PQN=led_PQN[-2,]
colnames(led_PQN)=led_PQN[1,]
led_PQN=led_PQN[-1,]
r=apply(led_PQN, 1,as.numeric)
r=t(r)
monster=as.factor(c(1:as.numeric(nrow(PLASMA_PQN))))
dose=as.numeric(c(rep(1, times=as.numeric(nrow(PLASMA_PQN)))))
led_list=list("monster"=monster, "dose"=dose)
led.eS=neweS(r, led_list)
trled.eS <- transeS(led.eS, tranpar$lambda, tranpar$alpha)
kostakis=exprs(trled.eS)
colnames(kostakis)=as.factor(colnames(kostakis))
colnames(kostakis)=colnames(led_PQN)
final=t(kostakis)
final=cbind(cutler[,1:2], final[,1:as.numeric(ncol(final))])
final=final[,-1]
row.names(final)=commodore[,1]
#write.csv(final, file="PVP_final_no2.02_noFRPAHVs.csv")

PVP2=final[-c(1:74, 147:256),-1]

PVP1=read_csv("D:/MEGA/Publication 2/PVP NB 19-02-2020_KL_NB2.csv")

C1C2=c(NA, NA, 1,1,1,NA,2,NA,1,2,NA,2,NA,2,NA,NA,2,2,1,NA,NA,1,NA,2,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
PNG=c(0,0,1,1,1,0,1,0,1,1,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0)

PVP=cbind(PVP1, PVP2)
PVP$`Plasma N (pg/ml)`[2]=NA
PVP$`MTY (pg/ml)`[2]=NA

X=final[,-1]
pca.ENP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.ENP, group = final[,1], title = "ALL SAMPLES", ind.names = FALSE,legend = T,cex=4, pch=c(15,16,17,18))

PVP_final_no2_02_noFRPAHVs <- final[which(final$`captain$GROUP`=="HV"),]
ENP_final_3_bothQCexp_24062020_aligned_exclusions_ALL_final_HV <- read_csv("ENP_final_3_bothQCexp_24062020_aligned_exclusions_ALL_final_allHV.csv")
SEQ=ENP_final_3_bothQCexp_24062020_aligned_exclusions_ALL_final_HV$`SEQUENCE NUMBER`[which(ENP_final_3_bothQCexp_24062020_aligned_exclusions_ALL_final_HV$`ENSAT-HT ID` %in% row.names(PVP_final_no2_02_noFRPAHVs))]
SEQ=ENP_final_3_bothQCexp_24062020_aligned_exclusions_ALL_final_HV$RUN[which(ENP_final_3_bothQCexp_24062020_aligned_exclusions_ALL_final_HV$`ENSAT-HT ID` %in% row.names(PVP_final_no2_02_noFRPAHVs))]


X=PVP_final_no2_02_noFRPAHVs[,-(1)]
Y=SEQ
Y[which(Y<median(Y))]=1
Y[which(!Y==1)]=2
pca.ENP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plotIndiv(pca.ENP, group = Y, title = "HV SAMPLES", ind.names = FALSE,legend = T,cex=4)

hay=function(x) {
  {
    samp12=x
    test <- samp12
    train <- setdiff(1:nrow(X), test)
    #design <- data.frame(patient=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ENSAT-HT ID`[train])
    #design2 <- data.frame(patient=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ENSAT-HT ID`[test])
    plsda.train <- plsda(X[train, ], Y[train], ncomp = 10, scale = F)
    perf.plsda.train <- perf(plsda.train, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 25, dist = "mahalanobis.dist")
    plsda.train <- plsda(X[train, ], Y[train], ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), scale = F)
    test.predict <- predict(plsda.train, X[test, ], dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1])))]                         #number of components
    well=Y[test]==Prediction
    rich=vip(plsda.train)
    wellrich=list(well, rich)
  }
  return(wellrich)
}

system.time( PREvPOST <- mclapply.hack( 1:(as.numeric(nrow(X))),
                                        hay)) 

#PREvPOST=lapply(1:(as.numeric(nrow(X))), hay)          #number of repeats

MCs=function(x) {
  return(PREvPOST[[x]][[1]])
}

MCPVP=lapply((1:(as.numeric(nrow(X)))), MCs)          #number of repeats
prevpost=unlist(MCPVP)
sum(prevpost==F)
jim=SEQ[which(prevpost==TRUE)]
jim[which(jim<median(SEQ))]
jim[which(jim>=median(SEQ))]
length(jim[which(jim<median(SEQ))])
length(jim[which(jim>=median(SEQ))])
((length(jim[which(jim<median(SEQ))])/length(which(SEQ<median(SEQ))))+(length(jim[which(jim>=median(SEQ))])/length(which(SEQ>=median(SEQ)))))/2
#result: 73% accuracy, after exclusion: 66% SEQ, 70% RUN

is.odd <- function(x) x %% 2 != 0
jake=cbind(SEQ, PVP_final_no2_02_noFRPAHVs)
jake=PVP_final_no2_02_noFRPAHVs[order(row.names(PVP_final_no2_02_noFRPAHVs)),]
row.names(jake)[which(prevpost==FALSE)]


lord=function(x) {
  {
    kel=PREvPOST[[x]][[2]][,as.numeric(ncol(PREvPOST[[x]][[2]]))]
  }
  return(kel)
}

kel=lapply((1:(length(prevpost))), lord)

output <- matrix(unlist(kel), ncol = (length(prevpost)), byrow = F)
row.names(output)=names(kel[[1]])

arth=function(x) {
  {
    impo=which(output[x,]>1)
  }
  return(length(impo))
}

as=sapply(1:as.numeric(nrow(output)), arth)
arthas=cbind(names(kel[[1]]), as)
arthas=apply(arthas, 2, as.numeric)

menet=function(x){
  median(output[x,])
}

hil=function(x){
  mad(output[x,], constant = 1)
}

menets=sapply(1:as.numeric(nrow(output)), menet)
hils=sapply(1:as.numeric(nrow(output)), hil)
arthas=cbind(arthas, menets, hils)
arthas=arthas[order(arthas[,3]),]

HV_SEQ=list(PREvPOST, sum(prevpost==F), arthas)
HV_RUN=list(PREvPOST, sum(prevpost==F), arthas)

comp.test=function(x){
  cor.test(X[,x], SEQ, method = "spearman")$p.value
}

all.cor=sapply(1:91, comp.test)
all.comp.p=all.cor
jason=p.adjust(all.cor, method = "fdr")
all.comp.fdr=jason
colnames(X)[which(all.cor<0.05)]
colnames(X)[which(jason<0.05)]
#"1.872658199" "1.903997967" "2.214349772" "2.324851518" "3.126405303" "3.185692061" "3.612111183"


#At baseline
library(mixOmics)
PRE=PVP[37:72,]
ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG=PRE
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`CENTER ID`)
levels(enp.fac)
levels(enp.fac)=c("1", "2","3", "4", "5", "6")
center=as.matrix(enp.fac)
center=as.numeric(center)
center=as.matrix(center)
colnames(center)=c("CENTER ID")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Sample Age in days`)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
SA=as.matrix(enp.fac)
SA=as.numeric(SA)
SA=as.matrix(SA)
colnames(SA)=c("Sample Age")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$GENDER)
levels(enp.fac)
levels(enp.fac)=c("0","1")
levels(enp.fac)=as.numeric(levels(enp.fac))
gender=as.matrix(enp.fac)
gender=as.numeric(gender)
gender=as.matrix(gender)
colnames(gender)=c("SEX")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Patient age at sampling`)
levels(enp.fac)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
YOB=as.matrix(enp.fac)
YOB=as.numeric(YOB)
YOB=as.matrix(YOB)
colnames(YOB)=c("Patient Age")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$RUN)
levels(enp.fac)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
RUN=as.matrix(enp.fac)
RUN=as.numeric(RUN)
RUN=as.matrix(RUN)
colnames(RUN)=c("RUN")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$SEQ)
levels(enp.fac)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
SEQ=as.matrix(enp.fac)
SEQ=as.numeric(SEQ)
SEQ=as.matrix(SEQ)
colnames(SEQ)=c("SEQ")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$BMI)
levels(enp.fac)=c(levels(enp.fac), "na")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
BMI=as.matrix(enp.fac)
BMI=as.numeric(BMI)
BMI=as.matrix(BMI)
colnames(BMI)=c("BMI")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ARTERIAL HYPERTENSION`)
levels(enp.fac)=c("0", "1")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
HT=as.matrix(enp.fac)
HT=as.numeric(HT)
HT=as.matrix(HT)
colnames(HT)=c("Arterial Hypertension")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$DM)
levels(enp.fac)=c("0", "1")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
DM=as.matrix(enp.fac)
DM=as.numeric(DM)
DM=as.matrix(DM)
colnames(DM)=c("Diabetes Mellitus")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`TUMOR LOCALIZATION`)
levels(enp.fac)=c("0", "1")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
TL=as.matrix(enp.fac)
TL=as.numeric(TL)
TL=as.matrix(TL)
colnames(TL)=c("Tumor Location")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`SECRETORY PHENOTYPE`)
levels(enp.fac)=c("1", "2")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
SP=as.matrix(enp.fac)
SP=as.numeric(SP)
SP=as.matrix(SP)
colnames(SP)=c("SP")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`SECRETORY PHENOTYPE2`)
levels(enp.fac)=c("0", "1")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
SP2=as.matrix(enp.fac)
SP2=as.numeric(SP2)
SP2=as.matrix(SP2)
colnames(SP2)=c("Secretory Phenotype")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$GENOTYPE)
levels(enp.fac)=c("1","1","2", "2", "2","2","2", "3")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
GENO=as.matrix(enp.fac)
GENO=as.numeric(GENO)
GENO=as.matrix(GENO)
colnames(GENO)=c("GENO")
GENO[GENO==3]=NA
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$GENOTYPE)
levels(enp.fac)=c("1","1","2", "2", "3","3","3", "4")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
GENO2=as.matrix(enp.fac)
GENO2=as.numeric(GENO2)
GENO2=as.matrix(GENO2)
colnames(GENO2)=c("GENO2")
GENO2[GENO2==4]=NA
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`TUMOR DIMENSIONS`)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
TS=as.matrix(enp.fac)
TS=as.numeric(TS)
TS=as.matrix(TS)
colnames(TS)=c("Tumor Size")
Plasma_total=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Plasma N (pg/ml)`+ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Plasma M1 (pg/ml)`+ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`MTY (pg/ml)`
Urine_total=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Urine free DA (ug/day)`+ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Urine Free EPI (ug/day)`+ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Urine Free NE (ug/day)`
enp.fac=as.factor(Plasma_total)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
Plasma=as.matrix(enp.fac)
Plasma=as.numeric(Plasma)
Plasma=as.matrix(Plasma)
colnames(Plasma)=c("Total Plasma Metanephrines")
Urine_total=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Urine free DA (ug/day)`+ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Urine Free EPI (ug/day)`+ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Urine Free NE (ug/day)`
enp.fac=as.factor(Urine_total)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
Urine=as.matrix(enp.fac)
Urine=as.numeric(Urine)
Urine=as.matrix(Urine)
colnames(Urine)=c("Total Urine Catecholamines")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Days between pre and post sampling`)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
dPP=as.matrix(enp.fac)
dPP=as.numeric(dPP)
dPP=as.matrix(dPP)
colnames(dPP)=c("days between Pre and Post sampling")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Collection time from surgery in days`)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
dSur=as.matrix(enp.fac)
dSur=as.numeric(dSur)
dSur=as.matrix(dSur)
colnames(dSur)=c("Days from surgery")
GYDR=center
GYDR[which(GYDR==1)]=1000
GYDR[which(GYDR<1000)]=0
GYDR[which(GYDR==1000)]=1
colnames(GYDR)="GYDR"
GYLU=center
GYLU[which(GYLU==2)]=1000
GYLU[which(GYLU<1000)]=0
GYLU[which(GYLU==1000)]=1
colnames(GYLU)="GYLU"
GYMU=center
GYMU[which(GYMU==3)]=1000
GYMU[which(GYMU<1000)]=0
GYMU[which(GYMU==1000)]=1
colnames(GYMU)="GYMU"
GYWU=center
GYWU[which(GYWU==4)]=1000
GYWU[which(GYWU<1000)]=0
GYWU[which(GYWU==1000)]=1
colnames(GYWU)="GYWU"
NLNI=center
NLNI[which(NLNI==5)]=1000
NLNI[which(NLNI<1000)]=0
NLNI[which(NLNI==1000)]=1
colnames(NLNI)="NLNI"
PLWW=center
PLWW[which(PLWW==6)]=1000
PLWW[which(PLWW<1000)]=0
PLWW[which(PLWW==1000)]=1
colnames(PLWW)="PLWW"
pos=GENO
pos[which(pos==2)]=1000
pos[which(pos<1000)]=0
pos[which(pos==1000)]=1
colnames(pos)="POSITIVE"
neg=GENO
neg[which(neg==1)]=1000
neg[which(neg<1000)]=0
neg[which(neg==1000)]=1
colnames(neg)="NEGATIVE"
C1=GENO2
C1[which(C1==2)]=1000
C1[which(C1<1000)]=0
C1[which(C1==1000)]=1
colnames(C1)="CLUSTER 1"
C2=GENO2
C2[which(C2==3)]=1000
C2[which(C2<1000)]=0
C2[which(C2==1000)]=1
colnames(C2)="CLUSTER 2"
WHYF=cbind(SA, gender, YOB, TS, Plasma, Urine, BMI, HT, DM, TL, SP2, pos, C1, C2, dSur, GYDR, GYLU, GYMU, GYWU, NLNI, PLWW, RUN, SEQ)
GENO2[which(GENO2==3)]=NA
C1=unmap(C1C2)[,1]
C2=unmap(C1C2)[,2]
SDH=c(0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0)
WHYF=cbind(SA, gender, YOB, TS, Plasma, Urine, BMI, HT, DM, TL, SP2, pos, GENO2, dSur, GYDR, GYLU, GYMU, GYWU, NLNI, PLWW, RUN, SEQ, C1, C2, PNG)
colnames(WHYF)=c("sample age", "sex", "patient age", "tumor size", "plasma metanephrines", "urine catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor location", "Secretory Phenotype", "Presence of mutation", "Mutation Cluster", "days before surgery", "GYDR", "GYLU", "GYMU", "GYWU", "NLNI", "PLWW", "RUN", "SEQ", "C1", "C2", "PNG")

X=PRE[,-(1:37)]
pca.ENP = mixOmics::pca(X, ncomp = 10, center = TRUE, scale = FALSE)
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
sex=gender
sex[which(sex==0)]="FEMALE"
sex[which(sex==1)]="MALE"
plotIndiv(pca.ENP, group = sex, title = "PPGL gender", ind.names = FALSE,legend = T,cex=7, comp = c(1,2), pch = c(17,19))
Y=SEQ
Y[which(Y<median(Y))]=1
Y[which(!Y==1)]=2
Y=SP2
plotIndiv(pca.ENP, group = Y, title = "PPGL Secretory Phenotype", ind.names = FALSE,legend = T,cex=7, comp = c(1,2), pch = c(17,19))

X=PRE[,-(1:37)]
X=as.matrix(X[,1:(as.numeric(ncol(X)))])
colnames(X)=round(as.numeric(colnames(X)), 3)

Y=WHYF[,-c(12,13,15,16,17,18,19,25)] #all
Y=WHYF[,c(1,14,20,21,22)] #tech
Y=WHYF[,-c(1,12,13,14,15,16,17,18,19,20,21,22,25)] #bio
Y=WHYF[,-c(1,2,3,12,13,14,15,16,17,18,19,20,21,22,25)] #clin
Y=WHYF[,-c(6,7,8,9,12,13,14,15,16,17,18,19,25)] #all no NA
Y=WHYF[,c(1,20,21,22)] #tech
Y=WHYF[,-c(1,6,7,8,9,12,13,14,15,16,17,18,19,20,21,22,25)] #bio no NA
Y=WHYF[,-c(1,2,3,6,7,8,9,12,13,14,15,16,17,18,19,20,21,22,25)] #clin no NA

Y2=Y[!(rowSums(is.na(Y))>0),]
X=X[!(rowSums(is.na(Y))>0),]
Y=Y2
pls.PvP <- pls(Y, X, ncomp = ncol(Y)-1, scale = T)
perf.pls.PvP <- perf(pls.PvP, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 50)
pls.PvP <- pls(Y, X, ncomp = min(which(perf.pls.PvP$Q2.total>0.0975)), scale = T)
vipsPvP_pls=vip(pls.PvP)
plotLoadings(pls.PvP, ndisplay = 30)
cim(pls.PvP, comp = 1, xlab = "NMR peaks", ylab = "Factors", margins = c(7, 7))
tune = tune.spls(Y, X, ncomp=ncol(Y)-1, test.keepX = rep(1:as.numeric(ncol(Y))), progressBar = TRUE, scale=T, validation="loo")
plot(tune)
tune$choice.keepX
toxicity.spls <- spls(Y, X, ncomp = ncol(Y)-1, keepX = tune$choice.keepX)
perf.spls.PvP <- perf(toxicity.spls, validation = "loo", auc = F, progressBar = T)
perf.spls.PvP$Q2.total

hay=function(x) {
  {
    samp12=x
    samp1234=samp12
    test <- samp1234
    train <- setdiff(1:nrow(X), test)
    #pls.train <- pls(Y[train,], X[train,], ncomp = 10, scale = T)
    pls.train <- pls(Y[train,], X[train,], ncomp = ncol(Y)-1, scale = T)
    perf.pls.train <- perf(pls.train, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 50)
    pls.train <- pls(Y[train, ], X[train,], ncomp = which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total)), scale = T)
    #pls.train <- pls(Y[train, ], X[train,], ncomp = min(which(perf.pls.train$Q2.total>0.0975)), scale = T)
    if (which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total))==1){
      RSS=rep(ncol(X)*(length(train)), ncol(X))
      RSS=colSums(do.call("rbind", replicate(length(train), RSS, simplify=F)))
    } else {
      test.predict.train <- predict(pls.train, Y[train,, drop = F])$predict[,,which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total))-1]
      rich = X[train,] - test.predict.train
      RSS=colSums(rich^2)
    }
    test.predict.test <- predict(pls.train, Y[test,, drop = F])$predict[,,which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total))]
    well = X[test, ] - test.predict.test
    PRESS=well^2
    PRESS1=colSums(do.call("rbind", replicate(length(train), PRESS, simplify=F)))
    Q2.total=1 -sum(PRESS1)/sum(RSS)
    wellrich=list(PRESS, RSS, perf.pls.train$Q2.total, Q2.total)
    #wellrich=perf.pls.train$Q2.total
  }
  return(wellrich)
}

PREvPOST=lapply(1:nrow(Y), hay)

allcolmean=do.call(rbind, replicate(nrow(Y), colMeans(X), simplify=FALSE))
TSS1=X-allcolmean

TSS = apply(TSS1, 2, function(x) {
  base::norm(x, type = "2")^2
})


PVP_Yres=function(x){
  PREvPOST[[x]][[1]]
}

PVP_Yres1=t(sapply(1:nrow(X),PVP_Yres))
PRESS1=PVP_Yres1
PRESS = apply(PRESS1, 2, function(x) {
  base::norm(x, type = "2")
})

Q2.total2 = 1 - sum(PRESS)/sum(TSS)

each.Q2=function(x){
  PREvPOST[[x]][[4]]
}

Q2.total=median(sapply(1:nrow(X),each.Q2))
Q2.total_mad=mad(sapply(1:nrow(X),each.Q2))

ea.Q2=function(x){
  PREvPOST[[x]][[3]]
}

Q2.ea=sapply(1:nrow(X),ea.Q2)
which(Q2.ea>0)


#permutation test

heave=function(x){
  X1=sample(nrow(X), replace=F)
  X=X[X1,]
}

set.seed(1)
newX=lapply(1:1000, heave)

Sheev=function(y){
  {
    X=newX[[y]]
    hay=function(x) {
      {
        samp12=x
        samp1234=samp12
        test <- samp1234
        train <- setdiff(1:nrow(X), test)
        pls.train <- pls(Y[train,], X[train,], ncomp = ncol(Y)-1, scale = T)
        perf.pls.train <- perf(pls.train, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 50)
        pls.train <- pls(Y[train, ], X[train,], ncomp = which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total)), scale = T)
        #pls.train <- pls(Y[train, ], X[train,], ncomp = min(which(perf.pls.train$Q2.total>0.0975)), scale = T)
        if (which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total))==1){
          RSS=rep(ncol(X)*(length(train)), ncol(X))
          RSS=colSums(do.call("rbind", replicate(length(train), RSS, simplify=F)))
        } else {
          test.predict.train <- predict(pls.train, Y[train,, drop = F])$predict[,,which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total))-1]
          rich = X[train,] - test.predict.train
          RSS=colSums(rich^2)
        }
        test.predict.test <- predict(pls.train, Y[test,, drop = F])$predict[,,which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total))]
        well = X[test, ] - test.predict.test
        PRESS=well^2
        PRESS1=colSums(do.call("rbind", replicate(length(train), PRESS, simplify=F)))
        Q2.total=1 -sum(PRESS1)/sum(RSS)
        wellrich=list(PRESS, RSS, perf.pls.train$Q2.total, Q2.total)
        #wellrich=perf.pls.train$Q2.total
      }
      return(wellrich)
    }
    
    PREvPOST=lapply(1:nrow(Y), hay)
    
    allcolmean=do.call(rbind, replicate(nrow(Y), colMeans(X), simplify=FALSE))
    TSS1=X-allcolmean
    
    TSS = apply(TSS1, 2, function(x) {
      base::norm(x, type = "2")^2
    })
    
    
    PVP_Yres=function(x){
      PREvPOST[[x]][[1]]
    }
    
    PVP_Yres1=t(sapply(1:nrow(X),PVP_Yres))
    PRESS1=PVP_Yres1
    PRESS = apply(PRESS1, 2, function(x) {
      base::norm(x, type = "2")
    })
    
    Q2.total2 = 1 - sum(PRESS)/sum(TSS)
    
    each.Q2=function(x){
      PREvPOST[[x]][[4]]
    }
    
    Q2.total=median(sapply(1:nrow(X),each.Q2))
    Q2.total_mad=mad(sapply(1:nrow(X),each.Q2))
    fbf1=list(Q2.total, Q2.total_mad, Q2.total2)
  }
  return(fbf1)
}


system.time( perm.all <- mclapply.hack( 1:1000 , Sheev ) )

perm.each=function(x) {
  perm.all[[x]][[3]]
}

perm.ea=sapply(1:1000, perm.each)
length(which(perm.ea>=Q2.total2))

atb_total=list(Q2.total2, perm.ea)
atb_tech=list(Q2.total2, perm.ea)
atb_bio=list(Q2.total2, perm.ea)
atb_clin=list(Q2.total2, perm.ea)
atb_total_noNA=list(Q2.total2, perm.ea)
atb_tech_noNA=list(Q2.total2, perm.ea)
atb_bio_noNA=list(Q2.total2, perm.ea)
atb_clin_noNA=list(Q2.total2, perm.ea)

#FACTOR PLSDA
#center: GYDR v PLWW
ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG=PRE
CEN=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG[which(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`CENTER ID`=="GYDR" | ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`CENTER ID`=="PLWW"),]
X=CEN[,-(1:37)]
Y=CEN$`CENTER ID`
test=CEN
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

SAm=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=SAm[,-(1:37)]
Y=SAm$`Sample Age in days`
Y=as.numeric(Y)
Y[which(Y<median(Y))]=1
Y[which(!Y==1)]=2
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)


SEX=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=SEX[,-(1:37)]
Y=SEX$GENDER
test=SEX
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

AGE=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=AGE[,-(1:37)]
Y=AGE$`Patient age at sampling`
Y[which(Y<45)]=1
Y[which(Y>=45)]=2
test=AGE
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

BMI25=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG[-which(is.na(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$BMI)),]
X=BMI25[,-(1:37)]
Y=BMI25$BMI
#X=X[-which(is.na(BMI25$BMI)),]
#Y=Y[-which(is.na(BMI25$BMI))]
Y[which(Y<25)]=1
Y[which(Y>=25)]=2
test=BMI25
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

AHT=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG[-which(is.na(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ARTERIAL HYPERTENSION`)),]
X=AHT[,-(1:37)]
Y=AHT$`ARTERIAL HYPERTENSION`
#X=X[-which(is.na(AHT$`ARTERIAL HYPERTENSION`)),]
#Y=Y[-which(is.na(AHT$`ARTERIAL HYPERTENSION`))]
test=AHT
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

Diab=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG[-which(is.na(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$DM)),]
X=Diab[,-(1:37)]
Y=Diab$DM
#X=X[-which(is.na(Diab$DM)),]
#Y=Y[-which(is.na(Diab$DM))]
test=Diab
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

TLoc=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=TLoc[,-(1:37)]
Y=TLoc$`TUMOR LOCALIZATION`
#X=X[-which(is.na(Diab$DM)),]
#Y=Y[-which(is.na(Diab$DM))]
test=TLoc
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)



SecPh2=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=SecPh2[,-(1:37)]
Y=SecPh2$`SECRETORY PHENOTYPE2`
#X=X[-which(is.na(Diab$DM)),]
test=SecPh2
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)
perf.plsda.PvP <- perf(plsda.PvP, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 25, dist = "mahalanobis.dist")
plsda.PvP <- plsda(X, Y, ncomp = min(which(perf.plsda.PvP$error.rate$BER[,1]==min(perf.plsda.PvP$error.rate$BER[,1]))), scale = F)
plotIndiv(plsda.PvP)
plotLoadings(plsda.PvP, ndisplay = 30)
plotVar(plsda.PvP)

Y=as.factor(Y)
levels(Y)=c("1", "2")
FC=rep(1,91)
delta=PRE[,-(1:37)]
pre_factor=WHYF[,12]
delta_1=delta[-which(is.na(pre_factor)),]
pre_factor_1=pre_factor[-which(is.na(pre_factor))]
delta_1=delta
pre_factor_1=pre_factor
X=delta_1
Y=pre_factor_1

comp.test=rep(1, 91)

for(i in 1:91){
  if(shapiro.test(X[which(Y==2),i])$p.value<0.05){
    big_av=median(X[which(Y==2),i])
  } else {
    big_av=mean(X[which(Y==2),i])
  }
  if(shapiro.test(X[which(Y==1),i])$p.value<0.05){
    sm_av=median(X[which(Y==1),i])
  } else {
    sm_av=mean(X[which(Y==1),i])
  }
  FC[i]=big_av/sm_av
}


GENO=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG[-which(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$GENOTYPE=="Unknown"),]
X=GENO[,-(1:37)]
Y=GENO$GENOTYPE
Y[which(Y=="negative")]="Negative"
#X=X[-which(Y=="Unknown"),]
#Y=Y[-which(Y=="Unknown")]
Y[-which(Y=="Negative")]="Positive"
test=GENO
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)


TSm5=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=TSm5[,-(1:37)]
Y=TSm5$`TUMOR DIMENSIONS`
Y=as.numeric(Y)
Y[which(Y<5)]=1
Y[which(Y>=5)]=2
test=TSm5
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)


PM=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=PM[,-(1:37)]
Y=Plasma_total
Y=as.numeric(Y)
Y[which(Y<median(Y))]=1
Y[which(!Y==1)]=2
test=PM
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

UC=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG[-which(is.na(Urine_total)),]
X=UC[,-(1:37)]
Y=Urine_total
#X=X[-which(is.na(Urine_total)),]
Y=Y[-which(is.na(Urine_total))]
Y=as.numeric(Y)
Y[which(Y<median(Y))]=1
Y[which(!Y==1)]=2
test=UC
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

RUNPRE=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=RUNPRE[,-(1:37)]
Y=RUN
Y=as.numeric(Y)
Y[which(Y<median(Y))]=1
Y[which(!Y==1)]=2
test=RUNPRE
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

SEQPRE=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=SEQPRE[,-(1:37)]
Y=SEQ
Y=as.numeric(Y)
Y[which(Y<median(Y))]=1
Y[which(!Y==1)]=2
test=SEQPRE
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

surPRE=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG[-which(is.na(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Collection time from surgery in days`)),]
X=surPRE[,-(1:37)]
Y=surPRE$`Collection time from surgery in days`
#X=X[-which(is.na(Y)),]
#Y=Y[-which(is.na(Y))]
Y=as.numeric(Y)
Y[which(Y<median(Y))]=1
Y[which(!Y==1)]=2
test=surPRE
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

c1c2PRE=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG[-which(is.na(C1C2)),]
X=c1c2PRE[,-(1:37)]
Y=C1C2
#X=X[-which(is.na(Y)),]
Y=Y[-which(is.na(Y))]
Y=as.numeric(Y)
test=c1c2PRE
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

PNGPRE=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=PNGPRE[,-(1:37)]
Y=PNG
test=PNGPRE
plsda.PvP = plsda(X,Y,ncomp=10, scale=F)

c1PRE=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=c1PRE[,-(1:37)]
Y=C1
#X=X[-which(is.na(Y)),]
#Y=Y[-which(is.na(Y))]
Y=as.numeric(Y)
test=c1PRE
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

c2PRE=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=c2PRE[,-(1:37)]
Y=C2
#X=X[-which(is.na(Y)),]
#Y=Y[-which(is.na(Y))]
Y=as.numeric(Y)
test=c2PRE
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

SDHPRE=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG
X=SDHPRE[,-(1:37)]
Y=SDH
#X=X[-which(is.na(Y)),]
#Y=Y[-which(is.na(Y))]
Y=as.numeric(Y)
test=SDHPRE
plsda.PvP <- plsda(X, Y, ncomp = 10, scale = F)

hay=function(x) {
  {
    samp12=x
    test <- samp12
    train <- setdiff(1:nrow(X), test)
    #design <- data.frame(patient=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ENSAT-HT ID`[train])
    #design2 <- data.frame(patient=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ENSAT-HT ID`[test])
    plsda.train <- plsda(X[train, ], Y[train], ncomp = 10, scale = F)
    perf.plsda.train <- perf(plsda.train, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 25, dist = "mahalanobis.dist")
    plsda.train <- plsda(X[train, ], Y[train], ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), scale = F)
    test.predict <- predict(plsda.train, X[test, ], dist = "mahalanobis.dist")     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1])))]                         #number of components
    well=Y[test]==Prediction
    rich=vip(plsda.train)
    wellrich=list(well, rich)
  }
  return(wellrich)
}

system.time( PREvPOST <- mclapply.hack( 1:(as.numeric(nrow(X))),
                                        hay)) 

#PREvPOST=lapply(1:(as.numeric(nrow(X))), hay)          #number of repeats

MCs=function(x) {
  return(PREvPOST[[x]][[1]])
}

MCPVP=lapply((1:(as.numeric(nrow(X)))), MCs)          #number of repeats
prevpost=unlist(MCPVP)
sum(prevpost==F)


#Acc:
(sum(prevpost==T)/length(prevpost==F))
(length(Y[which(Y[which(prevpost==TRUE)]=="no")]) + length(Y[which(Y[which(prevpost==TRUE)]=="yes")])) / (length(Y[which(Y=="no")]) + length(Y[which(Y=="yes")]))
(length(Y[which(Y[which(prevpost==TRUE)]=="GYDR")]) + length(Y[which(Y[which(prevpost==TRUE)]=="PLWW")])) / (length(Y[which(Y=="GYDR")]) + length(Y[which(Y=="PLWW")]))
(length(Y[which(Y[which(prevpost==TRUE)]=="FEMALE")]) + length(Y[which(Y[which(prevpost==TRUE)]=="MALE")])) / (length(Y[which(Y=="FEMALE")]) + length(Y[which(Y=="MALE")]))
(length(Y[which(Y[which(prevpost==TRUE)]==1)]) + length(Y[which(Y[which(prevpost==TRUE)]==2)])) / (length(Y[which(Y==1)]) + length(Y[which(Y==2)]))
(length(Y[which(Y[which(prevpost==TRUE)]=="Adrenal")]) + length(Y[which(Y[which(prevpost==TRUE)]=="Extra-adrenal")])) / (length(Y[which(Y=="Adrenal")]) + length(Y[which(Y=="Extra-adrenal")]))
(length(Y[which(Y[which(prevpost==TRUE)]=="Negative")]) + length(Y[which(Y[which(prevpost==TRUE)]=="Positive")])) / (length(Y[which(Y=="Negative")]) + length(Y[which(Y=="Positive")]))
(length(Y[which(Y[which(prevpost==TRUE)]=="Adrenergic")]) + length(Y[which(Y[which(prevpost==TRUE)]=="Nonadrenergic")])) / (length(Y[which(Y=="Adrenergic")]) + length(Y[which(Y=="Nonadrenergic")]))


#Bacc:
((length(Y[which(Y[which(prevpost==TRUE)]=="no")]) / length(Y[which(Y=="no")])) + (length(Y[which(Y[which(prevpost==TRUE)]=="yes")]) / length(Y[which(Y=="yes")])))/2
((length(Y[which(Y[which(prevpost==TRUE)]=="GYDR")]) / length(Y[which(Y=="GYDR")])) + (length(Y[which(Y[which(prevpost==TRUE)]=="PLWW")]) / length(Y[which(Y=="PLWW")])))/2
((length(Y[which(Y[which(prevpost==TRUE)]=="FEMALE")]) / length(Y[which(Y=="FEMALE")])) + (length(Y[which(Y[which(prevpost==TRUE)]=="MALE")]) / length(Y[which(Y=="MALE")])))/2
((length(Y[which(Y[which(prevpost==TRUE)]==1)]) / length(Y[which(Y==1)])) + (length(Y[which(Y[which(prevpost==TRUE)]==2)]) / length(Y[which(Y==2)])))/2
((length(Y[which(Y[which(prevpost==TRUE)]=="Adrenal")]) / length(Y[which(Y=="Adrenal")])) + (length(Y[which(Y[which(prevpost==TRUE)]=="Extra-adrenal")]) / length(Y[which(Y=="Extra-adrenal")])))/2
((length(Y[which(Y[which(prevpost==TRUE)]=="Negative")]) / length(Y[which(Y=="Negative")])) + (length(Y[which(Y[which(prevpost==TRUE)]=="Positive")]) / length(Y[which(Y=="Positive")])))/2
((length(Y[which(Y[which(prevpost==TRUE)]=="Adrenergic")]) / length(Y[which(Y=="Adrenergic")])) + (length(Y[which(Y[which(prevpost==TRUE)]=="Nonadrenergic")]) / length(Y[which(Y=="Nonadrenergic")])))/2
((length(Y[which(Y[which(prevpost==TRUE)]==0)]) / length(Y[which(Y==0)])) + (length(Y[which(Y[which(prevpost==TRUE)]==1)]) / length(Y[which(Y==1)])))/2



lord=function(x) {
  {
    kel=PREvPOST[[x]][[2]][,as.numeric(ncol(PREvPOST[[x]][[2]]))]
  }
  return(kel)
}

kel=lapply((1:(length(prevpost))), lord)

output <- matrix(unlist(kel), ncol = (length(prevpost)), byrow = F)
row.names(output)=names(kel[[1]])

arth=function(x) {
  {
    impo=which(output[x,]>1)
  }
  return(length(impo))
}

as=sapply(1:as.numeric(nrow(output)), arth)
arthas=cbind(names(kel[[1]]), as)
arthas=apply(arthas, 2, as.numeric)

menet=function(x){
  median(output[x,])
}

hil=function(x){
  mad(output[x,], constant = 1)
}

menets=sapply(1:as.numeric(nrow(output)), menet)
hils=sapply(1:as.numeric(nrow(output)), hil)
arthas=cbind(arthas, menets, hils)
arthas=arthas[order(arthas[,3]),]

CENTER1=list(PREvPOST, prevpost, arthas, plsda.PvP)
SA1=list(PREvPOST, prevpost, arthas, plsda.PvP)
SEX1=list(PREvPOST, prevpost, arthas, plsda.PvP)
AGE1.1=list(PREvPOST, prevpost, arthas, plsda.PvP)
BMI251=list(PREvPOST, prevpost, arthas, plsda.PvP)
AHT1=list(PREvPOST, prevpost, arthas, plsda.PvP)
DM1.1=list(PREvPOST, prevpost, arthas, plsda.PvP)
TL1=list(PREvPOST, prevpost, arthas, plsda.PvP)
SP2.1=list(PREvPOST, prevpost, arthas, plsda.PvP)
GENO1=list(PREvPOST, prevpost, arthas, plsda.PvP)
TS1.5=list(PREvPOST, prevpost, arthas, plsda.PvP)
PLASMA1=list(PREvPOST, prevpost, arthas, plsda.PvP)
URINE1=list(PREvPOST, prevpost, arthas, plsda.PvP)
RUNNUM1=list(PREvPOST, prevpost, arthas, plsda.PvP)
SEQNUM1=list(PREvPOST, prevpost, arthas, plsda.PvP)
surPRE1=list(PREvPOST, prevpost, arthas, plsda.PvP)
C1C21=list(PREvPOST, prevpost, arthas, plsda.PvP)
PNG1=list(PREvPOST, prevpost, arthas, plsda.PvP)
C111=list(PREvPOST, prevpost, arthas, plsda.PvP)
C211=list(PREvPOST, prevpost, arthas, plsda.PvP)
SDH1=list(PREvPOST, prevpost, arthas, plsda.PvP)

#PRE V POST
#PVP$BMI[which(PVP$`Time between sampling and weight in days`>10 | PVP$`Time between sampling and weight in days`<(-10))]=NA
PVP$BMI[which(is.na(PVP$BMI) & PVP$GROUP=="PPGL")]=NA
BMI1PRE=PVP[which(PVP$BMI<25 & PVP$GROUP=="PPGL"),]
BMI1POST=PVP[(which(PVP$BMI<25 & PVP$GROUP=="PPGL")-36),]
BMI1=rbind(BMI1POST, BMI1PRE)
BMI2PRE=PVP[which(PVP$BMI>=25 & PVP$GROUP=="PPGL"),]
BMI2POST=PVP[(which(PVP$BMI>=25 & PVP$GROUP=="PPGL")-36),]
BMI2=rbind(BMI2POST, BMI2PRE)
MALEPRE=PVP[which(PVP$GENDER=="MALE" & PVP$GROUP=="PPGL"),]
MALEPOST=PVP[(which(PVP$GENDER=="MALE" & PVP$GROUP=="PPGL")-36),]
MALE=rbind(MALEPOST, MALEPRE)
FEMALEPRE=PVP[which(PVP$GENDER=="FEMALE" & PVP$GROUP=="PPGL"),]
FEMALEPOST=PVP[(which(PVP$GENDER=="FEMALE" & PVP$GROUP=="PPGL")-36),]
FEMALE=rbind(FEMALEPOST, FEMALEPRE)
AGE1PRE=PVP[which(PVP$`Patient age at sampling`<45 & PVP$GROUP=="PPGL"),]
AGE1POST=PVP[(which(PVP$`Patient age at sampling`<45 & PVP$GROUP=="PPGL")-36),]
AGE1=rbind(AGE1POST, AGE1PRE)
AGE2PRE=PVP[which(PVP$`Patient age at sampling`>=45 & PVP$GROUP=="PPGL"),]
AGE2POST=PVP[(which(PVP$`Patient age at sampling`>=45 & PVP$GROUP=="PPGL")-36),]
AGE2=rbind(AGE2POST, AGE2PRE)
ADRPRE=PVP[which(PVP$`SECRETORY PHENOTYPE`=="Adrenergic" & PVP$GROUP=="PPGL"),]
ADRPOST=PVP[(which(PVP$`SECRETORY PHENOTYPE`=="Adrenergic" & PVP$GROUP=="PPGL")-36),]
ADR=rbind(ADRPOST, ADRPRE)
NADRPRE=PVP[which(PVP$`SECRETORY PHENOTYPE`=="Nonadrenergic" & PVP$GROUP=="PPGL"),]
NADRPOST=PVP[(which(PVP$`SECRETORY PHENOTYPE`=="Nonadrenergic" & PVP$GROUP=="PPGL")-36),]
NADR=rbind(NADRPOST, NADRPRE)
ADRPRE2=PVP[which(PVP$`SECRETORY PHENOTYPE2`=="Adrenergic" & PVP$GROUP=="PPGL"),]
ADRPOST2=PVP[(which(PVP$`SECRETORY PHENOTYPE2`=="Adrenergic" & PVP$GROUP=="PPGL")-36),]
ADR2=rbind(ADRPOST2, ADRPRE2)
NADRPRE2=PVP[which(PVP$`SECRETORY PHENOTYPE2`=="Nonadrenergic" & PVP$GROUP=="PPGL"),]
NADRPOST2=PVP[(which(PVP$`SECRETORY PHENOTYPE2`=="Nonadrenergic" & PVP$GROUP=="PPGL")-36),]
NADR2=rbind(NADRPOST2, NADRPRE2)
ADRENALPRE=PVP[which(PVP$`TUMOR LOCALIZATION`=="Adrenal" & PVP$GROUP=="PPGL"),]
ADRENALPOST=PVP[(which(PVP$`TUMOR LOCALIZATION`=="Adrenal" & PVP$GROUP=="PPGL")-36),]
ADRENAL=rbind(ADRENALPOST, ADRENALPRE)
TS1PRE=PVP[which(as.numeric(PVP$`TUMOR DIMENSIONS`)<4 & PVP$GROUP=="PPGL"),]
TS1POST=PVP[(which(as.numeric(PVP$`TUMOR DIMENSIONS`)<4 & PVP$GROUP=="PPGL")-36),]
TS1=rbind(TS1POST, TS1PRE)
TS2PRE=PVP[which(as.numeric(PVP$`TUMOR DIMENSIONS`)>=4 & PVP$GROUP=="PPGL"),]
TS2POST=PVP[(which(as.numeric(PVP$`TUMOR DIMENSIONS`)>=4 & PVP$GROUP=="PPGL")-36),]
TS2=rbind(TS2POST, TS2PRE)
TS15PRE=PVP[which(as.numeric(PVP$`TUMOR DIMENSIONS`)<5 & PVP$GROUP=="PPGL"),]
TS15POST=PVP[(which(as.numeric(PVP$`TUMOR DIMENSIONS`)<5 & PVP$GROUP=="PPGL")-36),]
TS15=rbind(TS15POST, TS15PRE)
TS25PRE=PVP[which(as.numeric(PVP$`TUMOR DIMENSIONS`)>=5 & PVP$GROUP=="PPGL"),]
TS25POST=PVP[(which(as.numeric(PVP$`TUMOR DIMENSIONS`)>=5 & PVP$GROUP=="PPGL")-36),]
TS25=rbind(TS25POST, TS25PRE)
HT1POST=PVP[which(PVP$`ARTERIAL HYPERTENSION`[37:72]=="yes" & PVP$`ARTERIAL HYPERTENSION`[1:36]=="no"),]
HT1PRE=PVP[which(PVP$`ARTERIAL HYPERTENSION`[37:72]=="yes" & PVP$`ARTERIAL HYPERTENSION`[1:36]=="no")+36,]
HT1=rbind(HT1POST, HT1PRE)
HT2POST=PVP[which(PVP$`ARTERIAL HYPERTENSION`[37:72]=="yes" & PVP$`ARTERIAL HYPERTENSION`[1:36]=="yes"),]
HT2PRE=PVP[which(PVP$`ARTERIAL HYPERTENSION`[37:72]=="yes" & PVP$`ARTERIAL HYPERTENSION`[1:36]=="yes")+36,]
HT2=rbind(HT2POST, HT2PRE)
HT3POST=PVP[which(PVP$`ARTERIAL HYPERTENSION`[37:72]=="no" & PVP$`ARTERIAL HYPERTENSION`[1:36]=="no"),]
HT3PRE=PVP[which(PVP$`ARTERIAL HYPERTENSION`[37:72]=="no" & PVP$`ARTERIAL HYPERTENSION`[1:36]=="no")+36,]
HT3=rbind(HT3POST, HT3PRE)
DM1POST=PVP[which(PVP$DM[37:72]=="no" & PVP$DM[1:36]=="no"),]
DM1PRE=PVP[which(PVP$DM[37:72]=="no" & PVP$DM[1:36]=="no")+36,]
DM1=rbind(DM1POST, DM1PRE)
DM2POST=PVP[which(PVP$DM[37:72]=="yes" & PVP$DM[1:36]=="no"),]
DM2PRE=PVP[which(PVP$DM[37:72]=="yes" & PVP$DM[1:36]=="no")+36,]
DM2=rbind(DM2POST, DM2PRE)
DM3POST=PVP[which(PVP$DM[37:72]=="yes" & PVP$DM[1:36]=="yes"),]
DM3PRE=PVP[which(PVP$DM[37:72]=="yes" & PVP$DM[1:36]=="yes")+36,]
DM3=rbind(DM3POST, DM3PRE)
SA1PRE=PVP[which(PRE$`Sample Age in days`<median(PRE$`Sample Age in days`) & PVP$GROUP=="PPGL"),]
SA1POST=PVP[(which(PRE$`Sample Age in days`<median(PRE$`Sample Age in days`) & PVP$GROUP=="PPGL")-36),]
SA1=rbind(SA1POST, SA1PRE)
SA2PRE=PVP[which(PRE$`Sample Age in days`>=median(PRE$`Sample Age in days`) & PVP$GROUP=="PPGL"),]
SA2POST=PVP[(which(PRE$`Sample Age in days`>=median(PRE$`Sample Age in days`) & PVP$GROUP=="PPGL")-36),]
SA2=rbind(SA2POST,SA2PRE)
TPP1PRE=PVP[which(PRE$`Days between pre and post sampling`<median(PRE$`Days between pre and post sampling`) & PVP$GROUP=="PPGL"),]
TPP1POST=PVP[(which(PRE$`Days between pre and post sampling`<median(PRE$`Days between pre and post sampling`) & PVP$GROUP=="PPGL")-36),]
TPP1=rbind(TPP1POST, TPP1PRE)
TPP2PRE=PVP[which(PRE$`Days between pre and post sampling`>=median(PRE$`Days between pre and post sampling`) & PVP$GROUP=="PPGL"),]
TPP2POST=PVP[(which(PRE$`Days between pre and post sampling`>=median(PRE$`Days between pre and post sampling`) & PVP$GROUP=="PPGL")-36),]
TPP2=rbind(TPP2POST,TPP2PRE)
GYDRPRE=PVP[which(PVP$`CENTER ID`=="GYDR" & PVP$GROUP=="PPGL"),]
GYDRPOST=PVP[(which(PVP$`CENTER ID`=="GYDR" & PVP$GROUP=="PPGL")-36),]
GYDR1=rbind(GYDRPOST, GYDRPRE)
PLWWPRE=PVP[which(PVP$`CENTER ID`=="PLWW" & PVP$GROUP=="PPGL"),]
PLWWPOST=PVP[(which(PVP$`CENTER ID`=="PLWW" & PVP$GROUP=="PPGL")-36),]
PLWW1=rbind(PLWWPOST, PLWWPRE)
C1PRE=PRE[which(C1C2==1),]
C1POST=PVP[which(C1C2==1),]
C11=rbind(C1POST, C1PRE)
C2PRE=PRE[which(C1C2==2),]
C2POST=PVP[which(C1C2==2),]
C21=rbind(C2POST, C2PRE)
POSPRE=PRE[which(PNG==1),]
POSPOST=PVP[which(PNG==1),]
POS1=rbind(POSPOST, POSPRE)
NEGPRE=PRE[which(PNG==0),]
NEGPOST=PVP[which(PNG==0),]
NEG1=rbind(NEGPOST, NEGPRE)

ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG=PVP
x <- ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG[,-(1:37)]
SEQ.fac=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$SEQ
SEQ.fac[SEQ.fac>=median(SEQ.fac)]=1
SEQ.fac[SEQ.fac>1]=0
y1=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$GROUP)
X <- x
enp_names=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ENSAT-HT ID`
row.names(X)=make.unique(enp_names)
Y <- as.factor(y1)
#levels(Y)=c("POST", "PPGL")

design <- data.frame(patient=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ENSAT-HT ID`)
plsda.PvP <- plsda(X, Y, ncomp = 10, multilevel = design, scale = F)
perf.plsda.PvP <- perf(plsda.PvP, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 25, dist = "mahalanobis.dist")
plsda.PvP <- plsda(X, Y, ncomp = min(which(perf.plsda.PvP$error.rate$BER[,1]==min(perf.plsda.PvP$error.rate$BER[,1]))), multilevel = design, scale = F)
vipsPvP=vip(plsda.PvP)
vipsPvP=vipsPvP[order(vipsPvP[,min(which(perf.plsda.PvP$error.rate$BER[,1]==min(perf.plsda.PvP$error.rate$BER[,1])))]),]
#write.csv(file="PvPvipsML.csv", vipsPvP)

tune.pca(X, ncomp=10, center = T, scale = F, multilevel = design)
pca.ENP = mixOmics::pca(X, ncomp = 3, center = TRUE, scale = FALSE, multilevel = design)
plot3dML=plotIndiv(pca.ENP, group = Y, ind.names = FALSE,legend = T,title = 'Pre vs. Post', style = '3d', cex=0.5)
levels(Y)=c("POST", "PRE")
plotIndiv(pca.ENP, group = Y, ind.names = F,legend = T,title = 'Pre vs. Post', cex=4, pch = c(15, 19), col.per.group = c("blue","darkorange2"))

hay=function(x) {
  {
    samp12=x
    samp34=samp12+(as.numeric(nrow(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG))/2)
    samp1234=c(samp12, samp34)    
    test <- samp1234
    train <- setdiff(1:nrow(X), test)
    design <- data.frame(patient=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ENSAT-HT ID`[train])
    design2 <- data.frame(patient=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ENSAT-HT ID`[test])
    plsda.train <- plsda(X[train, ], Y[train], ncomp = 10, multilevel = design, scale = F)
    perf.plsda.train <- perf(plsda.train, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 25, dist = "mahalanobis.dist")
    plsda.train <- plsda(X[train, ], Y[train], ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), multilevel = design, scale = F)
    test.predict <- predict(plsda.train, X[test, ], dist = "mahalanobis.dist", multilevel = design2)     #change distance
    Prediction <- test.predict$class$mahalanobis.dist[, min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1])))]                         #number of components
    well=Y[test]==Prediction
    rich=vip(plsda.train)
    wellrich=list(well, rich)
  }
  return(wellrich)
}

system.time( PREvPOST <- mclapply.hack( 1:(as.numeric(nrow(X))/2),
                                        hay)) 

#PREvPOST=lapply(1:(as.numeric(nrow(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG)/2)), hay)          #number of repeats

MCs=function(x) {
  return(PREvPOST[[x]][[1]])
}

MCPVP=lapply((1:(as.numeric(nrow(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG))/2)), MCs)          #number of repeats
prevpost=unlist(MCPVP)
sum(prevpost==F)
is.odd <- function(x) x %% 2 != 0
jake=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG[order(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ENSAT-HT ID`),]
jake$`ENSAT-HT ID`[which(prevpost==FALSE)[is.odd(which(prevpost==FALSE))]]
jake$SEQ[which(prevpost==FALSE)]
jake$`ENSAT-HT ID`[which(prevpost==FALSE)]
jake$GENDER[which(prevpost==FALSE)]
jake$`Sample Age in days`[which(prevpost==FALSE)]
jake$`Days between pre and post sampling`[which(prevpost==FALSE)]
(length(which(jake$GENDER[which(prevpost==FALSE)]=="FEMALE"))/length(which(jake$GENDER=="FEMALE")))/(length(which(jake$GENDER[which(prevpost==FALSE)]=="MALE"))/length(which(jake$GENDER=="MALE")))
(length(which(jake$`Sample Age in days`[which(prevpost==FALSE)]<median(jake$`Sample Age in days`)))/length(which(jake$`Sample Age in days`<median(jake$`Sample Age in days`))))/(length(which(jake$`Sample Age in days`[which(prevpost==FALSE)]>=median(jake$`Sample Age in days`)))/length(which(jake$`Sample Age in days`>=median(jake$`Sample Age in days`))))
(length(which(jake$`Days between pre and post sampling`[which(prevpost==FALSE)]<median(jake$`Days between pre and post sampling`)))/length(which(jake$`Days between pre and post sampling`<median(jake$`Days between pre and post sampling`))))/(length(which(jake$`Days between pre and post sampling`[which(prevpost==FALSE)]>=median(jake$`Days between pre and post sampling`)))/length(which(jake$`Days between pre and post sampling`>=median(jake$`Days between pre and post sampling`))))
median(jake$`Days between pre and post sampling`[which(prevpost==FALSE)])/median(jake$`Days between pre and post sampling`[which(prevpost==TRUE)])

jakeP=jake[which(jake$GROUP=="PPGL"),]
prevpostP=prevpost[seq(2,72,2)]
cen=as.factor(jakeP$`CENTER ID`)
levels(cen)=c("1", "2", "3", "4", "5", "6")
#levels(cen)=c("1", "2", "3", "4")
gen=as.factor(jakeP$GENDER)
levels(gen)=c("0", "1")
gen=as.numeric(as.character(gen))
age=jakeP$`Patient age at sampling`
run=jakeP$RUN
seq=jakeP$SEQ
sa=jakeP$`Sample Age in days`
dpp=jake$`Days between pre and post sampling`
dppP=dpp[seq(2,72,2)]
tl=as.factor(jakeP$`TUMOR LOCALIZATION`)
levels(tl)=c("0", "1")
tl=as.numeric(as.character(tl))
ts=as.numeric(jakeP$`TUMOR DIMENSIONS`)
sp=as.factor(jakeP$`SECRETORY PHENOTYPE2`)
levels(sp)=c("0", "1")
sp=as.numeric(as.character(sp))
pvp=as.factor(prevpostP)
levels(pvp)=c("0", "1")
pvp=as.numeric(as.character(pvp))
cc=as.factor(C1C2)
levels(cc)=c("0", "1")
cc=as.numeric(as.character(cc))
pn=as.factor(PNG)
levels(pn)=c("0", "1")
pn=as.numeric(as.character(pn))
c1=as.factor(C1)
levels(c1)=c("0","1")
c1=as.numeric(as.character(c1))
c2=as.factor(C2)
levels(c2)=c("0","1")
c2=as.numeric(as.character(c2))
mylogit <- glm(pvp ~ cen + gen + age + run + seq + sa + dppP + tl + ts + sp + c1 + c2, family = "binomial")
summary(mylogit)

bmi=as.numeric(jakeP$BMI)
dm=as.factor(jakeP$DM)
levels(dm)=c("0", "1")
aht=as.factor(jakeP$`ARTERIAL HYPERTENSION`)
levels(aht)=c("0", "1")
pvp=pvp[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
cen=cen[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
gen=gen[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
age=age[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
run=run[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
seq=seq[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
sa=sa[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
dppP=dppP[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
tl=tl[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
ts=ts[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
sp=sp[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
bmi2=bmi[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
dm2=dm[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
aht2=aht[-which(is.na(bmi) | is.na(dm) | is.na(aht) | is.na(cc))]
mylogit <- glm(pvp ~ cen + gen + age + run + seq + sa + dppP + tl + ts + sp + bmi2 + dm2 + aht2, family = "binomial")
summary(mylogit)

jakePO=jake[which(jake$GROUP=="POST"),]
prevpostPO=prevpost[seq(1,72,2)]
cen=as.factor(jakePO$`CENTER ID`)
levels(cen)=c("1", "2", "3", "4", "5", "6")
gen=as.factor(jakePO$GENDER)
levels(gen)=c("0", "1")
gen=as.numeric(as.character(gen))
age=jakePO$`Patient age at sampling`
run=jakePO$RUN
seq=jakePO$SEQ
sa=jakePO$`Sample Age in days`
dpp=jake$`Days between pre and post sampling`
dppP=dpp[seq(1,72,2)]
tl=as.factor(jakePO$`TUMOR LOCALIZATION`)
levels(tl)=c("0", "1")
tl=as.numeric(as.character(tl))
ts=as.numeric(jakePO$`TUMOR DIMENSIONS`)
sp=as.factor(jakePO$`SECRETORY PHENOTYPE2`)
levels(sp)=c("0", "1")
sp=as.numeric(as.character(sp))
pvp=as.factor(prevpostPO)
levels(pvp)=c("0", "1")
pvp=as.numeric(as.character(pvp))
mylogit <- glm(pvp ~ cen + gen + age + run + seq + sa + dppP, family = "binomial")
summary(mylogit)
bmi=as.numeric(jakePO$BMI)
dm=as.factor(jakePO$DM)
levels(dm)=c("0", "1")
aht=as.factor(jakePO$`ARTERIAL HYPERTENSION`)
levels(aht)=c("0", "1")
pvp=pvp[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
cen=cen[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
gen=gen[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
age=age[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
run=run[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
seq=seq[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
sa=sa[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
dppP=dppP[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
tl=tl[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
ts=ts[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
sp=sp[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
bmi2=bmi[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
dm2=dm[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
aht2=aht[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
mylogit <- glm(pvp ~ cen + gen + age + run + seq + sa + dppP + bmi2 + dm2 + aht2, family = "binomial")
summary(mylogit)

cen=as.factor(jake$`CENTER ID`)
levels(cen)=c("1", "2", "3", "4", "5", "6")
gen=as.factor(jake$GENDER)
levels(gen)=c("0", "1")
gen=as.numeric(as.character(gen))
age=jake$`Patient age at sampling`
run=jake$RUN
seq=jake$SEQ
sa=jake$`Sample Age in days`
dpp=jake$`Days between pre and post sampling`
tl=as.factor(jakeP$`TUMOR LOCALIZATION`)
levels(tl)=c("0", "1")
tl=as.numeric(as.character(tl))
tl=c(tl,tl)
ts=as.numeric(jakeP$`TUMOR DIMENSIONS`)
ts=c(ts,ts)
sp=as.factor(jakeP$`SECRETORY PHENOTYPE2`)
levels(sp)=c("0", "1")
sp=as.numeric(as.character(sp))
sp=c(sp,sp)
pvp=as.factor(prevpost)
levels(pvp)=c("0", "1")
pvp=as.numeric(as.character(pvp))
mylogit <- glm(pvp ~ cen + gen + age + run + seq + sa + dpp + tl + ts + sp, family = "binomial")
summary(mylogit)
bmi=as.numeric(jake$BMI)
dm=as.factor(jake$DM)
levels(dm)=c("0", "1")
aht=as.factor(jake$`ARTERIAL HYPERTENSION`)
levels(aht)=c("0", "1")
pvp=pvp[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
cen=cen[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
gen=gen[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
age=age[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
run=run[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
seq=seq[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
sa=sa[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
dpp=dpp[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
tl=tl[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
ts=ts[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
sp=sp[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
bmi2=bmi[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
dm2=dm[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
aht2=aht[-which(is.na(bmi) | is.na(dm) | is.na(aht))]
mylogit <- glm(pvp ~ cen + gen + age + run + seq + sa + dpp + tl + ts + sp + bmi2 + dm2 + aht2, family = "binomial")
summary(mylogit)

lord=function(x) {
  {
    kel=PREvPOST[[x]][[2]][,as.numeric(ncol(PREvPOST[[x]][[2]]))]
  }
  return(kel)
}

kel=lapply((1:(length(prevpost)/2)), lord)

output <- matrix(unlist(kel), ncol = (length(prevpost)/2), byrow = F)
row.names(output)=names(kel[[1]])

arth=function(x) {
  {
    impo=which(output[x,]>1)
  }
  return(length(impo))
}

as=sapply(1:as.numeric(nrow(output)), arth)
arthas=cbind(names(kel[[1]]), as)
arthas=apply(arthas, 2, as.numeric)

menet=function(x){
  median(output[x,])
}

hil=function(x){
  mad(output[x,], constant = 1)
}

menets=sapply(1:as.numeric(nrow(output)), menet)
hils=sapply(1:as.numeric(nrow(output)), hil)
arthas=cbind(arthas, menets, hils)
arthas=arthas[order(arthas[,3]),]

PVP1=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
BMI11=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
BMI21=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
MALE1=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
FEMALE1=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
AGE11=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
AGE21=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
ADR21=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
NADR21=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
ADRENAL1=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
TS115=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
TS215=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
HT11=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
HT21=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
DM11=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
SA11=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
SA21=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
TPP11=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
TPP21=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
PLWW11=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
GYDR11=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
C111=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
C211=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
POS11=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)
NEG11=list(PREvPOST, prevpost, arthas, vipsPvP, plsda.PvP)


#DELTA
delta=PVP[37:72,-(1:37)]-PVP[1:36,-(1:37)]


ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG=PRE
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`CENTER ID`)
levels(enp.fac)
levels(enp.fac)=c("1", "2","3", "4", "5", "6")
center=as.matrix(enp.fac)
center=as.numeric(center)
center=as.matrix(center)
colnames(center)=c("CENTER ID")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Sample Age in days`)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
SA=as.matrix(enp.fac)
SA=as.numeric(SA)
SA=as.matrix(SA)
colnames(SA)=c("Sample Age")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$GENDER)
levels(enp.fac)
levels(enp.fac)=c("0","1")
levels(enp.fac)=as.numeric(levels(enp.fac))
gender=as.matrix(enp.fac)
gender=as.numeric(gender)
gender=as.matrix(gender)
colnames(gender)=c("SEX")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Patient age at sampling`)
levels(enp.fac)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
YOB=as.matrix(enp.fac)
YOB=as.numeric(YOB)
YOB=as.matrix(YOB)
colnames(YOB)=c("Patient Age")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$RUN)
levels(enp.fac)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
RUN=as.matrix(enp.fac)
RUN=as.numeric(RUN)
RUN=as.matrix(RUN)
colnames(RUN)=c("RUN")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$SEQ)
levels(enp.fac)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
SEQ=as.matrix(enp.fac)
SEQ=as.numeric(SEQ)
SEQ=as.matrix(SEQ)
colnames(SEQ)=c("SEQ")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$BMI)
levels(enp.fac)=c(levels(enp.fac), "na")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
BMI=as.matrix(enp.fac)
BMI=as.numeric(BMI)
BMI=as.matrix(BMI)
colnames(BMI)=c("BMI")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ARTERIAL HYPERTENSION`)
levels(enp.fac)=c("0", "1")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
HT=as.matrix(enp.fac)
HT=as.numeric(HT)
HT=as.matrix(HT)
colnames(HT)=c("Arterial Hypertension")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$DM)
levels(enp.fac)=c("0", "1")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
DM=as.matrix(enp.fac)
DM=as.numeric(DM)
DM=as.matrix(DM)
colnames(DM)=c("Diabetes Mellitus")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`TUMOR LOCALIZATION`)
levels(enp.fac)=c("0", "1")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
TL=as.matrix(enp.fac)
TL=as.numeric(TL)
TL=as.matrix(TL)
colnames(TL)=c("Tumor Location")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`SECRETORY PHENOTYPE`)
levels(enp.fac)=c("1", "2")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
SP=as.matrix(enp.fac)
SP=as.numeric(SP)
SP=as.matrix(SP)
colnames(SP)=c("SP")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`SECRETORY PHENOTYPE2`)
levels(enp.fac)=c("0", "1")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
SP2=as.matrix(enp.fac)
SP2=as.numeric(SP2)
SP2=as.matrix(SP2)
colnames(SP2)=c("Secretory Phenotype")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$GENOTYPE)
levels(enp.fac)=c("1","1","2", "2", "2","2","2", "3")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
GENO=as.matrix(enp.fac)
GENO=as.numeric(GENO)
GENO=as.matrix(GENO)
colnames(GENO)=c("GENO")
GENO[GENO==3]=NA
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$GENOTYPE)
levels(enp.fac)=c("1","1","2", "2", "3","3","3", "4")
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
GENO2=as.matrix(enp.fac)
GENO2=as.numeric(GENO2)
GENO2=as.matrix(GENO2)
colnames(GENO2)=c("GENO2")
GENO2[GENO2==4]=NA
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`TUMOR DIMENSIONS`)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
TS=as.matrix(enp.fac)
TS=as.numeric(TS)
TS=as.matrix(TS)
colnames(TS)=c("Tumor Size")
Plasma_total=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Plasma N (pg/ml)`+ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Plasma M1 (pg/ml)`+ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`MTY (pg/ml)`
Urine_total=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Urine free DA (ug/day)`+ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Urine Free EPI (ug/day)`+ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Urine Free NE (ug/day)`
enp.fac=as.factor(Plasma_total)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
Plasma=as.matrix(enp.fac)
Plasma=as.numeric(Plasma)
Plasma=as.matrix(Plasma)
colnames(Plasma)=c("Total Plasma Metanephrines")
Urine_total=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Urine free DA (ug/day)`+ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Urine Free EPI (ug/day)`+ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Urine Free NE (ug/day)`
enp.fac=as.factor(Urine_total)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
Urine=as.matrix(enp.fac)
Urine=as.numeric(Urine)
Urine=as.matrix(Urine)
colnames(Urine)=c("Total Urine Catecholamines")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Days between pre and post sampling`)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
dPP=as.matrix(enp.fac)
dPP=as.numeric(dPP)
dPP=as.matrix(dPP)
colnames(dPP)=c("days between Pre and Post sampling")
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`Collection time from surgery in days`)
meisio=as.numeric(levels(enp.fac))
levels(enp.fac)=meisio
levels(enp.fac)=as.numeric(levels(enp.fac))
dSur=as.matrix(enp.fac)
dSur=as.numeric(dSur)
dSur=as.matrix(dSur)
colnames(dSur)=c("Days from surgery")
GYDR=center
GYDR[which(GYDR==1)]=1000
GYDR[which(GYDR<1000)]=0
GYDR[which(GYDR==1000)]=1
colnames(GYDR)="GYDR"
GYLU=center
GYLU[which(GYLU==2)]=1000
GYLU[which(GYLU<1000)]=0
GYLU[which(GYLU==1000)]=1
colnames(GYLU)="GYLU"
GYMU=center
GYMU[which(GYMU==3)]=1000
GYMU[which(GYMU<1000)]=0
GYMU[which(GYMU==1000)]=1
colnames(GYMU)="GYMU"
GYWU=center
GYWU[which(GYWU==4)]=1000
GYWU[which(GYWU<1000)]=0
GYWU[which(GYWU==1000)]=1
colnames(GYWU)="GYWU"
NLNI=center
NLNI[which(NLNI==5)]=1000
NLNI[which(NLNI<1000)]=0
NLNI[which(NLNI==1000)]=1
colnames(NLNI)="NLNI"
PLWW=center
PLWW[which(PLWW==6)]=1000
PLWW[which(PLWW<1000)]=0
PLWW[which(PLWW==1000)]=1
colnames(PLWW)="PLWW"
pos=GENO
pos[which(pos==2)]=1000
pos[which(pos<1000)]=0
pos[which(pos==1000)]=1
colnames(pos)="POSITIVE"
neg=GENO
neg[which(neg==1)]=1000
neg[which(neg<1000)]=0
neg[which(neg==1000)]=1
colnames(neg)="NEGATIVE"
WHYF=cbind(SA, gender, YOB, TS, Plasma, Urine, BMI, HT, DM, TL, SP2, pos, C1, C2, dPP, GYDR, GYLU, GYMU, GYWU, NLNI, PLWW)
GENO2[which(GENO2==3)]=NA
WHYF=cbind(SA, gender, YOB, TS, Plasma, Urine, BMI, HT, DM, TL, SP2, pos, GENO2, dPP, GYDR, GYLU, GYMU, GYWU, NLNI, PLWW, C1, C2, PNG)
colnames(WHYF)=c("sample age", "sex", "patient age", "tumor size", "plasma metanephrines", "urine catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor location", "Secretory Phenotype", "Presence of mutation", "Mutation Cluster", "days between sampling", "GYDR", "GYLU", "GYMU", "GYWU", "NLNI", "PLWW", "C1", "C2", "PNG")

X=data.frame(delta$`2.313307518`, delta$`2.37023004`, delta$`2.388646009`, delta$`4.132772125`, delta$`2.262080201`, delta$`3.136864883`, delta$`5.219707859`, delta$`5.227224754`, delta$`3.548012672`, delta$`3.126405303`, delta$`3.161934796`, delta$`3.984713634`, delta$`4.080120359`, delta$`4.093821199`, delta$`4.107654959`, delta$`4.120725314`, delta$`3.012515074`, delta$`3.346299278`, delta$`3.041061161`, delta$`1.996338747`, delta$`2.059835912`, delta$`3.31178287`, delta$`2.355572113`, delta$`6.891572492`, delta$`7.185187396`, delta$`3.177317058`, delta$`3.283913914`, delta$`3.670008318`)
X=as.matrix(X[,1:(as.numeric(ncol(X)))])
colnames(X)=c("2.313", "2.370", "2.389", "4.133", "2.262", "3.137","5.220", "5.227", "3.548", "3.126", "3.162", "3.985", "4.080", "4.094", "4.108", "4.121", "3.013","3.346", "3.041", "1.996", "2.060", "3.312", "2.356", "6.892", "7.185", "3.177", "3.284", "3.670")
#colnames(X)=c("3-hydroxybutyrate1", "3-hydroxybutyrate2", "3-hydroxybutyrate3", "3-OH-bu/Proline", "Acetoacetate", "Dimethyl sulfone", "Glucose1", "Glucose2","Glycine", "Histidine", "Lactate", "Lactate/Proline/3-OH-bu", "Lysine", "Methanol", "Ornithine", "Proline1", "Proline2", "Proline3","Pyruvate", "Tyrosine1", "Tyrosine2", "Unknown1", "Unknown2", "Unknown3")

Y=WHYF[,-c(1,12,13,15,16,17,18,19,23)]  #ALL
Y=WHYF[,-c(1,12,13,14,15,16,17,18,19,20,23)]  #BIOLOGICAL
Y=WHYF[,-c(1,2,3,12,13,14,15,16,17,18,19,20,23)]  #CLINICAL
Y=WHYF[,-c(1,6,7,8,9,12,13,15,16,17,18,19,23)] #ALL NONA
Y=WHYF[,-c(1,6,7,8,9,12,13,14,15,16,17,18,19,20,23)]  #BIOLOGICAL NONA
Y=WHYF[,-c(1,6,7,8,9,2,3,12,13,14,15,16,17,18,19,20,23)]  #CLINICAL NONA

Y2=Y[!(rowSums(is.na(Y))>0),]
X=X[!(rowSums(is.na(Y))>0),]
Y=Y2
pls.PvP <- pls(Y, X, ncomp = ncol(Y)-1, scale = T)
perf.pls.PvP <- perf(pls.PvP, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 50)
pls.PvP <- pls(Y, X, ncomp = min(which(perf.pls.PvP$Q2.total>0.0975)), scale = T)
vipsPvP_pls=vip(pls.PvP)
plotLoadings(pls.PvP, ndisplay = 30)
plotVar(pls.PvP)
plotIndiv(pls.PvP)
plotArrow(pls.PvP)
cim(pls.PvP, comp = 1, xlab = "NMR peaks", ylab = "Factors", margins = c(7, 7))
tune = tune.spls(Y, X, ncomp=ncol(Y)-1, test.keepX = rep(1:12), progressBar = TRUE, scale=T, validation="loo")
plot(tune)
tune$choice.keepX
toxicity.spls <- spls(Y, X, ncomp = ncol(Y)-1, keepX = tune$choice.keepX)
perf.spls.PvP <- perf(toxicity.spls, validation = "loo", auc = F, progressBar = T)


hay=function(x) {
  {
    samp12=x
    samp1234=samp12
    test <- samp1234
    train <- setdiff(1:nrow(X), test)
    #pls.train <- pls(Y[train,], X[train,], ncomp = 10, scale = T)
    pls.train <- pls(Y[train,], X[train,], ncomp = ncol(Y)-1, scale = T)
    perf.pls.train <- perf(pls.train, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 50)
    pls.train <- pls(Y[train, ], X[train,], ncomp = which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total)), scale = T)
    #pls.train <- pls(Y[train, ], X[train,], ncomp = min(which(perf.pls.train$Q2.total>0.0975)), scale = T)
    if (which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total))==1){
      RSS=rep(ncol(X)*(length(train)), ncol(X))
      RSS=colSums(do.call("rbind", replicate(length(train), RSS, simplify=F)))
    } else {
      test.predict.train <- predict(pls.train, Y[train,, drop = F])$predict[,,which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total))-1]
      rich = X[train,] - test.predict.train
      RSS=colSums(rich^2)
    }
    test.predict.test <- predict(pls.train, Y[test,, drop = F])$predict[,,which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total))]
    well = X[test, ] - test.predict.test
    PRESS=well^2
    PRESS1=colSums(do.call("rbind", replicate(length(train), PRESS, simplify=F)))
    Q2.total=1 -sum(PRESS1)/sum(RSS)
    wellrich=list(PRESS, RSS, perf.pls.train$Q2.total, Q2.total)
    #wellrich=perf.pls.train$Q2.total
  }
  return(wellrich)
}

PREvPOST=lapply(1:nrow(Y), hay)

allcolmean=do.call(rbind, replicate(nrow(Y), colMeans(X), simplify=FALSE))
TSS1=X-allcolmean

TSS = apply(TSS1, 2, function(x) {
  base::norm(x, type = "2")^2
})


PVP_Yres=function(x){
  PREvPOST[[x]][[1]]
}

PVP_Yres1=t(sapply(1:nrow(X),PVP_Yres))
PRESS1=PVP_Yres1
PRESS = apply(PRESS1, 2, function(x) {
  base::norm(x, type = "2")
})

Q2.total2 = 1 - sum(PRESS)/sum(TSS)

each.Q2=function(x){
  PREvPOST[[x]][[4]]
}

Q2.total=median(sapply(1:nrow(X),each.Q2))
Q2.total_mad=mad(sapply(1:nrow(X),each.Q2))

ea.Q2=function(x){
  PREvPOST[[x]][[3]]
}

Q2.ea=sapply(1:nrow(X),ea.Q2)
which(Q2.ea>0)

#permutation test

heave=function(x){
  X1=sample(nrow(X), replace=F)
  X=X[X1,]
}

set.seed(1)
newX=lapply(1:1000, heave)

Sheev=function(y){
  {
    X=newX[[y]]
    hay=function(x) {
      {
        samp12=x
        samp1234=samp12
        test <- samp1234
        train <- setdiff(1:nrow(X), test)
        pls.train <- pls(Y[train,], X[train,], ncomp = ncol(Y)-1, scale = T)
        perf.pls.train <- perf(pls.train, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 50)
        pls.train <- pls(Y[train, ], X[train,], ncomp = which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total)), scale = T)
        #pls.train <- pls(Y[train, ], X[train,], ncomp = min(which(perf.pls.train$Q2.total>0.0975)), scale = T)
        if (which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total))==1){
          RSS=rep(ncol(X)*(length(train)), ncol(X))
          RSS=colSums(do.call("rbind", replicate(length(train), RSS, simplify=F)))
        } else {
          test.predict.train <- predict(pls.train, Y[train,, drop = F])$predict[,,which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total))-1]
          rich = X[train,] - test.predict.train
          RSS=colSums(rich^2)
        }
        test.predict.test <- predict(pls.train, Y[test,, drop = F])$predict[,,which(perf.pls.train$Q2.total==max(perf.pls.train$Q2.total))]
        well = X[test, ] - test.predict.test
        PRESS=well^2
        PRESS1=colSums(do.call("rbind", replicate(length(train), PRESS, simplify=F)))
        Q2.total=1 -sum(PRESS1)/sum(RSS)
        wellrich=list(PRESS, RSS, perf.pls.train$Q2.total, Q2.total)
        #wellrich=perf.pls.train$Q2.total
      }
      return(wellrich)
    }
    
    PREvPOST=lapply(1:nrow(Y), hay)
    
    allcolmean=do.call(rbind, replicate(nrow(Y), colMeans(X), simplify=FALSE))
    TSS1=X-allcolmean
    
    TSS = apply(TSS1, 2, function(x) {
      base::norm(x, type = "2")^2
    })
    
    
    PVP_Yres=function(x){
      PREvPOST[[x]][[1]]
    }
    
    PVP_Yres1=t(sapply(1:nrow(X),PVP_Yres))
    PRESS1=PVP_Yres1
    PRESS = apply(PRESS1, 2, function(x) {
      base::norm(x, type = "2")
    })
    
    Q2.total2 = 1 - sum(PRESS)/sum(TSS)
    
    each.Q2=function(x){
      PREvPOST[[x]][[4]]
    }
    
    Q2.total=median(sapply(1:nrow(X),each.Q2))
    Q2.total_mad=mad(sapply(1:nrow(X),each.Q2))
    fbf1=list(Q2.total, Q2.total_mad, Q2.total2)
  }
  return(fbf1)
}


#perm.all=lapply(1:100, Sheev)
set.seed(1)
system.time( perm.all <- mclapply.hack( 1:100 , Sheev ) )

perm.each=function(x) {
  perm.all[[x]][[3]]
}

perm.ea=sapply(1:100, perm.each)
length(which(perm.ea>=Q2.total2))

total=list(Q2.total2, perm.ea)
bio=list(Q2.total2, perm.ea)
clin=list(Q2.total2, perm.ea)
total_noNA=list(Q2.total2, perm.ea)
bio_noNA=list(Q2.total2, perm.ea)
clin_noNA=list(Q2.total2, perm.ea)
