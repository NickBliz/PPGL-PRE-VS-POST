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

row.names(commodore)=commodore[,1]
commodore=commodore[,-1]
final=commodore

library(sjstats)
hay=function(x) {
  cv(final[which(final$`captain$GROUP`=="QC"),2:(as.numeric(ncol(final)))][,x])
}
cvQC=sapply((1:(as.numeric(ncol(final))-1)), hay)
cvQC=2*cvQC
hay=function(x) {
  cv(final[which(final$`captain$GROUP`=="HV"),2:(as.numeric(ncol(final)))][,x])
}
cvHV=sapply((1:(as.numeric(ncol(final))-1)), hay)

library(PASWR2)
difff=cvHV-cvQC
z.test(difff, sigma.x = sqrt( sd(difff)^2/length(difff) ), mu=0)


hay=function(x) {
  cv(final[which(final$`captain$GROUP`=="PPGL" | final$`captain$GROUP`=="POST"),2:(as.numeric(ncol(final)))][,x])
}
cvP=sapply((1:(as.numeric(ncol(final))-1)), hay)
hay=function(x) {
  cv(final[which(final$`captain$GROUP`=="PPGL"),2:(as.numeric(ncol(final)))][,x])
}
cvPP=sapply((1:(as.numeric(ncol(final))-1)), hay)
hay=function(x) {
  cv(final[which(final$`captain$GROUP`=="POST"),2:(as.numeric(ncol(final)))][,x])
}
cvPO=sapply((1:(as.numeric(ncol(final))-1)), hay)


PVP2=final[-c(1:74, 147:256),-1]

PVP1=read_csv("D:/MEGA/Publication 2/PVP NB 19-02-2020_KL_NB2.csv")

PVP=cbind(PVP1, PVP2)
PVP$`Plasma N (pg/ml)`[2]=NA
PVP$`MTY (pg/ml)`[2]=NA

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
C1C2=c(NA, NA, 1,1,1,NA,2,NA,1,2,NA,2,NA,2,NA,NA,2,2,1,NA,NA,1,NA,2,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
PNG=c(0,0,1,1,1,0,1,0,1,1,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0)
C1=unmap(C1C2)[,1]
C2=unmap(C1C2)[,2]
WHYF=cbind(SA, gender, YOB, TS, Plasma, Urine, BMI, HT, DM, TL, SP2, pos, C1, C2, dSur, GYDR, GYLU, GYMU, GYWU, NLNI, PLWW, RUN, SEQ, C1, C2, PNG)


#UVA T-TEST
delta=PRE[,-(1:37)]
WHYF=cbind(center, SA, gender, YOB, TS, Plasma, Urine, BMI, HT, DM, TL, SP2, pos, GENO2, dSur, RUN, SEQ, C1, C2, PNG)
all.comp.p=matrix(1, nrow=20, ncol=91)
colnames(all.comp.p)=colnames(X)
row.names(all.comp.p)=c("center", "sample age", "sex", "patient age", "tumor size", "plasma metanephrines", "urine catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor location", "Secretory Phenotype", "Presence of mutation", "Mutation Cluster", "days before surgery", "RUN", "SEQ", "C1", "C2", "PNG")
all.comp.fdr=all.comp.p
all.comp.fc=all.comp.p
comp.test=rep(1, 91)
FC=rep(1,91)
pre_factor=WHYF[,1]
#delta_1=delta[-which(is.na(pre_factor)),]
#pre_factor_1=pre_factor[-which(is.na(pre_factor))]
delta_1=delta
pre_factor_1=pre_factor
X=delta_1
Y=pre_factor_1

for(i in 1:91){
  comp.test[i]=wilcox.test(X[which(Y=="1"),i], X[which(Y=="6"),i])$p.value
    big_av=median(X[which(Y=="6"),i], na.rm = T)
    sm_av=median(X[which(Y=="1"),i], na.rm = T)
    FC[i]=big_av/sm_av
}

jason=p.adjust(comp.test, method = "fdr")
all.comp.p[1,]=comp.test
all.comp.fdr[1,]=jason
all.comp.fc[1,]=FC

which(jason<0.05)


delta=PRE[,-(1:37)]
pre_factor=WHYF[,7]
delta_1=delta[-which(is.na(pre_factor)),]
pre_factor_1=pre_factor[-which(is.na(pre_factor))]
#delta_1=delta
#pre_factor_1=pre_factor
X=delta_1
Y=pre_factor_1

comp.test=rep(1, 91)

for(i in 1:91){
  comp.test[i]=wilcox.test(X[which(Y<median(Y)),i], X[which(Y>=median(Y)),i])$p.value
    big_av=median(X[which(Y>=median(Y)),i], na.rm = T)
    sm_av=median(X[which(Y<median(Y)),i], na.rm = T)
  FC[i]=big_av/sm_av
}

jason=p.adjust(comp.test, method = "fdr")
all.comp.p[7,]=comp.test
all.comp.fdr[7,]=jason
all.comp.fc[7,]=FC

which(jason<0.05)

delta=PRE[,-(1:37)]
pre_factor=WHYF[,5]
delta_1=delta[-which(is.na(pre_factor)),]
pre_factor_1=pre_factor[-which(is.na(pre_factor))]
delta_1=delta
pre_factor_1=pre_factor
X=delta_1
Y=pre_factor_1

comp.test=rep(1, 91)

for(i in 1:91){
  comp.test[i]=wilcox.test(X[which(Y<5),i], X[which(Y>=5),i])$p.value
    big_av=median(X[which(Y>=5),i], na.rm = T)
    sm_av=median(X[which(Y<5),i], na.rm = T)
  FC[i]=big_av/sm_av
}

jason=p.adjust(comp.test, method = "fdr")
all.comp.p[5,]=comp.test
all.comp.fdr[5,]=jason
all.comp.fc[5,]=FC

which(jason<0.05)

delta=PRE[,-(1:37)]
pre_factor=WHYF[,4]
delta_1=delta[-which(is.na(pre_factor)),]
pre_factor_1=pre_factor[-which(is.na(pre_factor))]
delta_1=delta
pre_factor_1=pre_factor
X=delta_1
Y=pre_factor_1

comp.test=rep(1, 91)

for(i in 1:91){
  comp.test[i]=wilcox.test(X[which(Y<45),i], X[which(Y>=45),i])$p.value
    big_av=median(X[which(Y>=45),i], na.rm = T)
    sm_av=median(X[which(Y<45),i], na.rm = T)
  FC[i]=big_av/sm_av
}

jason=p.adjust(comp.test, method = "fdr")
all.comp.p[4,]=comp.test
all.comp.fdr[4,]=jason
all.comp.fc[4,]=FC

which(jason<0.05)

delta=PRE[,-(1:37)]
pre_factor=WHYF[,8]
delta_1=delta[-which(is.na(pre_factor)),]
pre_factor_1=pre_factor[-which(is.na(pre_factor))]
#delta_1=delta
#pre_factor_1=pre_factor
X=delta_1
Y=pre_factor_1

comp.test=rep(1, 91)

for(i in 1:91){
  comp.test[i]=wilcox.test(X[which(Y<25),i], X[which(Y>=25),i])$p.value
    big_av=median(X[which(Y>=25),i], na.rm = T)
    sm_av=median(X[which(Y<25),i], na.rm = T)
  FC[i]=big_av/sm_av
}

jason=p.adjust(comp.test, method = "fdr")
all.comp.p[8,]=comp.test
all.comp.fdr[8,]=jason
all.comp.fc[8,]=FC

which(jason<0.05)

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
  comp.test[i]=wilcox.test(X[which(Y==0),i], X[which(Y==1),i])$p.value
    big_av=median(X[which(Y==1),i], na.rm = T)
    sm_av=median(X[which(Y==0),i], na.rm = T)
  FC[i]=big_av/sm_av
}

jason=p.adjust(comp.test, method = "fdr")
all.comp.p[10,]=comp.test
all.comp.fdr[10,]=jason
all.comp.fc[10,]=FC

which(jason<0.05)

delta=PRE[,-(1:37)]
pre_factor=WHYF[,14]
delta_1=delta[-which(is.na(pre_factor)),]
pre_factor_1=pre_factor[-which(is.na(pre_factor))]
#delta_1=delta
#pre_factor_1=pre_factor
X=delta_1
Y=pre_factor_1

comp.test=rep(1, 91)

for(i in 1:91){
  comp.test[i]=wilcox.test(X[which(Y==2),i], X[which(Y==3),i])$p.value
    big_av=median(X[which(Y==3),i], na.rm = T)
    sm_av=median(X[which(Y==2),i], na.rm = T)
  FC[i]=big_av/sm_av
}

jason=p.adjust(comp.test, method = "fdr")
all.comp.p[14,]=comp.test
all.comp.fdr[14,]=jason
all.comp.fc[14,]=FC

which(jason<0.05)

AT_BASELINE_T_TEST=list(all.comp.p, all.comp.fdr, all.comp.fc)
write.csv(file="atb_t_fc_new.csv", all.comp.fc)
write.csv(file="atb_t_p_new.csv", all.comp.p)
write.csv(file="atb_t_fdr_new.csv", all.comp.fdr)

#UVA SPEARMAN
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
WHYF=cbind(SA, gender, YOB, TS, Plasma, Urine, BMI, HT, DM, TL, SP2, pos, GENO2, dSur, GYDR, GYLU, GYMU, GYWU, NLNI, PLWW, RUN, SEQ)

all.comp.p=matrix(1, nrow=22, ncol=91)
colnames(all.comp.p)=colnames(X)
row.names(all.comp.p)=c("sample age", "sex", "patient age", "tumor size", "plasma metanephrines", "urine catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor location", "Secretory Phenotype", "Presence of mutation", "Mutation Cluster", "days before surgery", "GYDR", "GYLU", "GYMU", "GYWU", "NLNI", "PLWW", "RUN", "SEQ")
all.comp.fdr=all.comp.p
all.comp.rho=all.comp.p

abcd=14

delta=PRE[,-(1:37)]
pre_factor=WHYF[,abcd]
delta_1=delta[-which(is.na(pre_factor)),]
pre_factor_1=pre_factor[-which(is.na(pre_factor))]
#delta_1=delta
#pre_factor_1=pre_factor
X=delta_1
Y=pre_factor_1

comp.test=function(x){
  cor.test(X[,x], Y, method = "spearman")$p.value
}

all.cor=sapply(1:91, comp.test)
all.comp.p[abcd,]=all.cor
jason=p.adjust(all.cor, method = "fdr")
all.comp.fdr[abcd,]=jason

comp.test=function(x){
  cor.test(X[,x], Y, method = "spearman")$estimate
}

all.cor=sapply(1:91, comp.test)

all.comp.rho[abcd,]=all.cor

which(jason<0.05)

AT_BASELINE_Spearman=list(all.comp.p, all.comp.fdr, all.comp.rho)
write.csv(file="atb_s_rho_new.csv", all.comp.rho)
write.csv(file="atb_s_p_new.csv", all.comp.p)
write.csv(file="atb_s_fdr_new.csv", all.comp.fdr)



#PRE V POST
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


all.comp.p=matrix(1, nrow=24, ncol=91)
colnames(all.comp.p)=colnames(X)
row.names(all.comp.p)=c("ALL", "BMI < 25 KG/M2", "BMI >= 25 KG/M2", "MALE", "FEMALE", "AGE < 45 YRS", "AGE >= 45 YRS", "ADRENERGIC", "NONADRENERGIC", "ADRENAL", "TUMOR SIZE < 5CM", "TUMOR SIZE >= 5CM", "HYPERTENSION YES PRE NO POST", "HYPERTENSION YES PRE YES POST", "DIABETES NO PRE NO POST", "GYDR", "PLWW", "SA1", "SA2", "TPP1", "TPP2", "C1", "C2", "POS")
all.comp.fdr=all.comp.p
all.comp.fc=all.comp.p
comp.test=rep(1, 91)
AFC=rep(1,91)

delta=C21[,-(1:37)]
pre_factor=C21$GROUP
abcd=23
delta_1=delta
pre_factor_1=pre_factor
X=delta_1
Y=pre_factor_1


for(i in 1:91){
  comp.test[i]=wilcox.test(X[which(Y=="POST"),i], X[which(Y=="PPGL"),i], paired = T)$p.value
}

FC=X[which(Y=="PPGL"),]/X[which(Y=="POST"),]
FC[(FC==0)]=NA
FC[(FC==Inf)]=NA

for(i in 1:91){
    AFC[i]=median(FC[,i], na.rm = T)
}

jason=p.adjust(comp.test, method = "fdr")
which(jason<0.05)
all.comp.p[abcd,]=comp.test
all.comp.fdr[abcd,]=jason
all.comp.fc[abcd,]=AFC

PVP_T_TEST=list(all.comp.p, all.comp.fdr, all.comp.fc)
write.csv(file="PVP_T_TEST_fc_new.csv", PVP_T_TEST[[3]])
write.csv(file="PVP_T_TEST_p_new.csv", PVP_T_TEST[[1]])
write.csv(file="PVP_T_TEST_fdr_new.csv", PVP_T_TEST[[2]])

#DELTA

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
C1C2=c(NA, NA, 1,1,1,NA,2,NA,1,2,NA,2,NA,2,NA,NA,2,2,1,NA,NA,1,NA,2,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
PNG=c(0,0,1,1,1,0,1,0,1,1,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0)
C1=unmap(C1C2)[,1]
C2=unmap(C1C2)[,2]
dPP=PVP$`Days between pre and post sampling`[1:36]


delta=PVP[37:72,-(1:37)]-PVP[1:36,-(1:37)]
WHYF=cbind(SA, PLWW, gender, YOB, TS, Plasma, Urine, BMI, HT, DM, TL, SP2, dPP, C1, C2, PNGF)
#colnames(WHYF)=c("Patient Age", "Tumor Size", "Total Plasma Metanephrines", "Total Urine Catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor Location", "Secretory Phenotype", "Days between Pre and Post sampling", "C1", "C2", "Presence of mutation")
colnames(WHYF)=c("Sample Age", "PLWW", "Sex", "Patient Age", "Tumor Size", "Total Plasma Metanephrines", "Total Urine Catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor Location", "Secretory Phenotype", "Days between Pre and Post sampling", "C1", "C2", "Presence of mutation")
X=data.frame(delta$`2.313307518`, delta$`2.37023004`, delta$`4.132772125`, delta$`2.262080201`, delta$`3.136864883`, delta$`5.219707859`, delta$`5.227224754`, delta$`3.548012672`, delta$`3.126405303`, delta$`3.984713634`, delta$`4.080120359`, delta$`4.093821199`, delta$`4.107654959`, delta$`2.997167774`, delta$`3.012515074`, delta$`3.346299278`, delta$`3.041061161`, delta$`1.996338747`, delta$`2.059835912`, delta$`3.31178287`, delta$`2.355572113`, delta$`2.388646009`, delta$`6.891572492`, delta$`7.185187396`, delta$`3.161934796`, delta$`3.283913914`, delta$`3.670008318`)
X=as.matrix(X[,1:(as.numeric(ncol(X)))])
#colnames(X)=c("2.313", "2.370", "2.389", "4.133", "2.262", "3.137","5.220", "5.227", "3.548", "3.126", "3.162", "3.985", "4.080", "4.094", "4.108", "4.121", "3.013","3.346", "3.041", "1.996", "2.060", "3.312", "2.356", "6.892", "7.185", "3.177", "3.284", "3.670")
colnames(X)=c("3-Hydroxybutyrate1", "3-Hydroxybutyrate2",  "3-Hydroxybutyrate/Proline", "Acetoacetate", "Dimethyl sulfone", "Glucose1", "Glucose2","Glycine", "Histidine/Phenylalanine", "Histidine/Phenylalanine/Serine", "Lactate1", "Lactate2", "Lactate3",  "Lysine1", "Lysine2", "Methanol", "Ornithine", "Proline1", "Proline2", "Proline3","Pyruvate", "Succinate/3-hydroxybutyrate", "Tyrosine1", "Tyrosine2", "Unknown1", "Unknown2", "Unknown3")
#X=X[,order(colnames(X))]
rquery.cormat3(X)

inspect_cor_plot=rquery.cormat3(X)
all.comp.fdr=inspect_cor_plot$p
all.comp.rho=inspect_cor_plot$r
colnames(all.comp.fdr)=row.names(all.comp.fdr)=colnames(all.comp.rho)=row.names(all.comp.rho)=colnames(X)

DELTA_Spearman=list(all.comp.p, all.comp.fdr, all.comp.rho)
write.csv(file="DELTA_Spearman_rho_corplot.csv", all.comp.rho)
write.csv(file="DELTA_Spearman_fdr_corplot.csv", all.comp.fdr)

X_1=X

X=X[!rowSums(is.na(WHYF)),]
WHYF=WHYF[!rowSums(is.na(WHYF)),]
test.all=cbind(WHYF, X)
rquery.cormat2(test.all)
inspect_cor_plot=rquery.cormat2(test.all)
all.comp.fdr=inspect_cor_plot$p
all.comp.rho=inspect_cor_plot$r
colnames(all.comp.fdr)=row.names(all.comp.fdr)=colnames(all.comp.rho)=row.names(all.comp.rho)=colnames(test.all)

DELTA_Spearman=list(all.comp.p, all.comp.fdr, all.comp.rho)
write.csv(file="DELTA_Spearman_rho_corplot2.csv", all.comp.rho)


all.comp.p=matrix(1, nrow=16, ncol=27)
colnames(all.comp.p)=colnames(X)
row.names(all.comp.p)=c("Sample Age", "PLWW", "Sex", "Patient Age", "Tumor Size", "Total Plasma Metanephrines", "Total Urine Catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor Location", "Secretory Phenotype", "Days between Pre and Post sampling", "C1", "C2", "Presence of mutation")
all.comp.fdr=all.comp.p
all.comp.rho=all.comp.p

abcd=10

pre_factor=WHYF[,abcd]
X=X_1[-which(is.na(pre_factor)),]
pre_factor_1=pre_factor[-which(is.na(pre_factor))]
Y=pre_factor_1

#pre_factor=WHYF[,abcd]
#pre_factor_1=pre_factor
#Y=pre_factor_1
#X=X_1

comp.test=function(x){
  cor.test(X[,x], Y, method = "spearman")$p.value
}

all.cor=sapply(1:27, comp.test)
all.comp.p[abcd,]=all.cor
jason=p.adjust(all.cor, method = "fdr")
all.comp.fdr[abcd,]=jason

comp.test=function(x){
  cor.test(X[,x], Y, method = "spearman")$estimate
}

all.cor=sapply(1:27, comp.test)

all.comp.rho[abcd,]=all.cor

which(jason<0.05)

test4=all.comp.rho
test5=all.comp.p
test6=all.comp.fdr

write.csv(file="DELTA_Spearman_rho.csv", test4)
write.csv(file="DELTA_Spearman_fdr.csv", test6)



ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG=PRE[which(PRE$GENDER=="FEMALE"),]
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`CENTER ID`)
levels(enp.fac)
levels(enp.fac)=c("1", "2","3", "4", "5")
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
PLWW=center
PLWW[which(PLWW==5)]=1000
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
C1C2=c(NA, NA, 1,1,1,NA,2,NA,1,2,NA,2,NA,2,NA,NA,2,2,1,NA,NA,1,NA,2,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
C1C2F=C1C2[which(PRE$GENDER=="FEMALE")]
PNG=c(0,0,1,1,1,0,1,0,1,1,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0)
PNGF=PNG[which(PRE$GENDER=="FEMALE")]
C1=unmap(C1C2F)[,1]
C2=unmap(C1C2F)[,2]
dPP=PRE$`Days between pre and post sampling`[which(PRE$GENDER=="FEMALE")]



delta=FEMALE[28:54,-(1:37)]-FEMALE[1:27,-(1:37)]
WHYF=cbind(SA, PLWW, YOB, TS, Plasma, Urine, BMI, HT, DM, TL, SP2, dPP, C1, C2, PNGF)
#colnames(WHYF)=c("Patient Age", "Tumor Size", "Total Plasma Metanephrines", "Total Urine Catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor Location", "Secretory Phenotype", "Days between Pre and Post sampling", "C1", "C2", "Presence of mutation")
colnames(WHYF)=c("Sample Age", "PLWW", "Patient Age", "Tumor Size", "Total Plasma Metanephrines", "Total Urine Catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor Location", "Secretory Phenotype", "Days between Pre and Post sampling", "C1", "C2", "Presence of mutation")
X=data.frame(delta$`2.313307518`, delta$`2.37023004`, delta$`4.132772125`, delta$`2.262080201`, delta$`3.136864883`, delta$`5.219707859`, delta$`5.227224754`, delta$`3.548012672`, delta$`3.126405303`, delta$`3.984713634`, delta$`4.080120359`, delta$`4.093821199`, delta$`4.107654959`, delta$`2.997167774`, delta$`3.012515074`, delta$`3.346299278`, delta$`3.041061161`, delta$`1.996338747`, delta$`2.059835912`, delta$`3.31178287`, delta$`2.355572113`, delta$`2.388646009`, delta$`6.891572492`, delta$`7.185187396`, delta$`3.161934796`, delta$`3.283913914`, delta$`3.670008318`)
X=as.matrix(X[,1:(as.numeric(ncol(X)))])
X_1=X
#colnames(X)=c("2.313", "2.370", "2.389", "4.133", "2.262", "3.137","5.220", "5.227", "3.548", "3.126", "3.162", "3.985", "4.080", "4.094", "4.108", "4.121", "3.013","3.346", "3.041", "1.996", "2.060", "3.312", "2.356", "6.892", "7.185", "3.177", "3.284", "3.670")
colnames(X)=c("3-Hydroxybutyrate1", "3-Hydroxybutyrate2",  "3-Hydroxybutyrate/Proline", "Acetoacetate", "Dimethyl sulfone", "Glucose1", "Glucose2","Glycine", "Histidine/Phenylalanine", "Histidine/Phenylalanine/Serine", "Lactate1", "Lactate2", "Lactate3",  "Lysine1", "Lysine2", "Methanol", "Ornithine", "Proline1", "Proline2", "Proline3","Pyruvate", "Succinate/3-hydroxybutyrate", "Tyrosine1", "Tyrosine2", "Unknown1", "Unknown2", "Unknown3")
#X=X[,order(colnames(X))]
rquery.cormat3(X)

all.comp.p=matrix(1, nrow=15, ncol=27)
colnames(all.comp.p)=colnames(X)
row.names(all.comp.p)=c("Sample Age", "PLWW", "Patient Age", "Tumor Size", "Total Plasma Metanephrines", "Total Urine Catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor Location", "Secretory Phenotype", "Days between Pre and Post sampling", "C1", "C2", "Presence of mutation")
all.comp.fdr=all.comp.p
all.comp.rho=all.comp.p

abcd=8

pre_factor=WHYF[,abcd]
X=X_1[-which(is.na(pre_factor)),]
pre_factor_1=pre_factor[-which(is.na(pre_factor))]
Y=pre_factor_1

pre_factor=WHYF[,abcd]
pre_factor_1=pre_factor
Y=pre_factor_1
X=X_1

comp.test=function(x){
  cor.test(X[,x], Y, method = "spearman")$p.value
}

all.cor=sapply(1:27, comp.test)
all.comp.p[abcd,]=all.cor
jason=p.adjust(all.cor, method = "fdr")
all.comp.fdr[abcd,]=jason

comp.test=function(x){
  cor.test(X[,x], Y, method = "spearman")$estimate
}

all.cor=sapply(1:27, comp.test)

all.comp.rho[abcd,]=all.cor

which(jason<0.05)

test4=all.comp.rho
test5=all.comp.p
test6=all.comp.fdr

write.csv(file="DELTA_Spearman_rho_FEMALE.csv", test4)
write.csv(file="DELTA_Spearman_fdr_FEMALE.csv", test6)


ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG=PRE[which(PRE$BMI<25),]
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`CENTER ID`)
levels(enp.fac)
levels(enp.fac)=c("1", "2","3", "4", "5")
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
GYMU=center
GYMU[which(GYMU==2)]=1000
GYMU[which(GYMU<1000)]=0
GYMU[which(GYMU==1000)]=1
colnames(GYMU)="GYMU"
GYWU=center
GYWU[which(GYWU==3)]=1000
GYWU[which(GYWU<1000)]=0
GYWU[which(GYWU==1000)]=1
colnames(GYWU)="GYWU"
NLNI=center
NLNI[which(NLNI==4)]=1000
NLNI[which(NLNI<1000)]=0
NLNI[which(NLNI==1000)]=1
colnames(NLNI)="NLNI"
PLWW=center
PLWW[which(PLWW==5)]=1000
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
C1C2=c(NA, NA, 1,1,1,NA,2,NA,1,2,NA,2,NA,2,NA,NA,2,2,1,NA,NA,1,NA,2,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
C1C2=C1C2[which(PRE$BMI<25)]
PNG=c(0,0,1,1,1,0,1,0,1,1,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0)
PNG=PNG[which(PRE$BMI<25)]
C1=unmap(C1C2)[,1]
C2=unmap(C1C2)[,2]
dPP=PVP$`Days between pre and post sampling`[which(PRE$BMI<25)]


delta=BMI1[16:30,-(1:37)]-BMI1[1:15,-(1:37)]
WHYF=cbind(PLWW, gender, YOB, TS, Plasma, Urine, BMI, HT, DM, TL, SP2, dPP, C1, C2, PNG)
#colnames(WHYF)=c("Patient Age", "Tumor Size", "Total Plasma Metanephrines", "Total Urine Catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor Location", "Secretory Phenotype", "Days between Pre and Post sampling", "C1", "C2", "Presence of mutation")
colnames(WHYF)=c("PLWW", "Sex", "Patient Age", "Tumor Size", "Total Plasma Metanephrines", "Total Urine Catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor Location", "Secretory Phenotype", "Days between Pre and Post sampling", "C1", "C2", "Presence of mutation")
X=data.frame(delta$`2.313307518`, delta$`2.37023004`, delta$`4.132772125`, delta$`2.262080201`, delta$`3.136864883`, delta$`5.219707859`, delta$`5.227224754`, delta$`3.548012672`, delta$`3.126405303`, delta$`3.984713634`, delta$`4.080120359`, delta$`4.093821199`, delta$`4.107654959`, delta$`2.997167774`, delta$`3.012515074`, delta$`3.346299278`, delta$`3.041061161`, delta$`1.996338747`, delta$`2.059835912`, delta$`3.31178287`, delta$`2.355572113`, delta$`2.388646009`, delta$`6.891572492`, delta$`7.185187396`, delta$`3.161934796`, delta$`3.283913914`, delta$`3.670008318`)
X=as.matrix(X[,1:(as.numeric(ncol(X)))])
X_1=X
#colnames(X)=c("2.313", "2.370", "2.389", "4.133", "2.262", "3.137","5.220", "5.227", "3.548", "3.126", "3.162", "3.985", "4.080", "4.094", "4.108", "4.121", "3.013","3.346", "3.041", "1.996", "2.060", "3.312", "2.356", "6.892", "7.185", "3.177", "3.284", "3.670")
colnames(X)=c("3-Hydroxybutyrate1", "3-Hydroxybutyrate2",  "3-Hydroxybutyrate/Proline", "Acetoacetate", "Dimethyl sulfone", "Glucose1", "Glucose2","Glycine", "Histidine/Phenylalanine", "Histidine/Phenylalanine/Serine", "Lactate1", "Lactate2", "Lactate3",  "Lysine1", "Lysine2", "Methanol", "Ornithine", "Proline1", "Proline2", "Proline3","Pyruvate", "Succinate/3-hydroxybutyrate", "Tyrosine1", "Tyrosine2", "Unknown1", "Unknown2", "Unknown3")
#X=X[,order(colnames(X))]
rquery.cormat3(X)



all.comp.p=matrix(1, nrow=15, ncol=27)
colnames(all.comp.p)=colnames(X)
row.names(all.comp.p)=c("PLWW", "Sex", "Patient Age", "Tumor Size", "Total Plasma Metanephrines", "Total Urine Catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor Location", "Secretory Phenotype", "Days between Pre and Post sampling", "C1", "C2", "Presence of mutation")
all.comp.fdr=all.comp.p
all.comp.rho=all.comp.p

abcd=15

pre_factor=WHYF[,abcd]
X=X_1[-which(is.na(pre_factor)),]
pre_factor_1=pre_factor[-which(is.na(pre_factor))]
Y=pre_factor_1

pre_factor=WHYF[,abcd]
pre_factor_1=pre_factor
Y=pre_factor_1
X=X_1

comp.test=function(x){
  cor.test(X[,x], Y, method = "spearman")$p.value
}

all.cor=sapply(1:27, comp.test)
all.comp.p[abcd,]=all.cor
jason=p.adjust(all.cor, method = "fdr")
all.comp.fdr[abcd,]=jason

comp.test=function(x){
  cor.test(X[,x], Y, method = "spearman")$estimate
}

all.cor=sapply(1:27, comp.test)

all.comp.rho[abcd,]=all.cor

which(jason<0.05)

test4=all.comp.rho
test5=all.comp.p
test6=all.comp.fdr

write.csv(file="DELTA_Spearman_rho_BMI1.csv", test4)
write.csv(file="DELTA_Spearman_fdr_BMI1.csv", test6)


ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG=PRE[which(PRE$`Days between pre and post sampling`>=median(PRE$`Days between pre and post sampling`)),]
enp.fac=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`CENTER ID`)
levels(enp.fac)
levels(enp.fac)=c("1", "2","3", "4")
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
GYWU=center
GYWU[which(GYWU==3)]=1000
GYWU[which(GYWU<1000)]=0
GYWU[which(GYWU==1000)]=1
colnames(GYWU)="GYWU"
PLWW=center
PLWW[which(PLWW==4)]=1000
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
C1C2=c(NA, NA, 1,1,1,NA,2,NA,1,2,NA,2,NA,2,NA,NA,2,2,1,NA,NA,1,NA,2,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
C1C2=C1C2[which(PRE$`Days between pre and post sampling`>=median(PRE$`Days between pre and post sampling`))]
PNG=c(0,0,1,1,1,0,1,0,1,1,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0)
PNG=PNG[which(PRE$`Days between pre and post sampling`>=median(PRE$`Days between pre and post sampling`))]
C1=unmap(C1C2)[,1]
C2=unmap(C1C2)[,2]
dPP=PVP$`Days between pre and post sampling`[which(PRE$`Days between pre and post sampling`>=median(PRE$`Days between pre and post sampling`))]


delta=TPP2[19:36,-(1:37)]-TPP2[1:18,-(1:37)]
WHYF=cbind(PLWW, gender, YOB, TS, Plasma, Urine, BMI, HT, DM, TL, SP2, dPP, C1, C2, PNG)
#colnames(WHYF)=c("Patient Age", "Tumor Size", "Total Plasma Metanephrines", "Total Urine Catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor Location", "Secretory Phenotype", "Days between Pre and Post sampling", "C1", "C2", "Presence of mutation")
colnames(WHYF)=c("PLWW", "Sex", "Patient Age", "Tumor Size", "Total Plasma Metanephrines", "Total Urine Catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor Location", "Secretory Phenotype", "Days between Pre and Post sampling", "C1", "C2", "Presence of mutation")
X=data.frame(delta$`2.313307518`, delta$`2.37023004`, delta$`4.132772125`, delta$`2.262080201`, delta$`3.136864883`, delta$`5.219707859`, delta$`5.227224754`, delta$`3.548012672`, delta$`3.126405303`, delta$`3.984713634`, delta$`4.080120359`, delta$`4.093821199`, delta$`4.107654959`, delta$`2.997167774`, delta$`3.012515074`, delta$`3.346299278`, delta$`3.041061161`, delta$`1.996338747`, delta$`2.059835912`, delta$`3.31178287`, delta$`2.355572113`, delta$`2.388646009`, delta$`6.891572492`, delta$`7.185187396`, delta$`3.161934796`, delta$`3.283913914`, delta$`3.670008318`)
X=as.matrix(X[,1:(as.numeric(ncol(X)))])
X_1=X
#colnames(X)=c("2.313", "2.370", "2.389", "4.133", "2.262", "3.137","5.220", "5.227", "3.548", "3.126", "3.162", "3.985", "4.080", "4.094", "4.108", "4.121", "3.013","3.346", "3.041", "1.996", "2.060", "3.312", "2.356", "6.892", "7.185", "3.177", "3.284", "3.670")
colnames(X)=c("3-Hydroxybutyrate1", "3-Hydroxybutyrate2",  "3-Hydroxybutyrate/Proline", "Acetoacetate", "Dimethyl sulfone", "Glucose1", "Glucose2","Glycine", "Histidine/Phenylalanine", "Histidine/Phenylalanine/Serine", "Lactate1", "Lactate2", "Lactate3",  "Lysine1", "Lysine2", "Methanol", "Ornithine", "Proline1", "Proline2", "Proline3","Pyruvate", "Succinate/3-hydroxybutyrate", "Tyrosine1", "Tyrosine2", "Unknown1", "Unknown2", "Unknown3")
#X=X[,order(colnames(X))]
rquery.cormat3(X)



all.comp.p=matrix(1, nrow=15, ncol=27)
colnames(all.comp.p)=colnames(X)
row.names(all.comp.p)=c("PLWW", "Sex", "Patient Age", "Tumor Size", "Total Plasma Metanephrines", "Total Urine Catecholamines", "BMI", "Hypertension", "Diabetes", "Tumor Location", "Secretory Phenotype", "Days between Pre and Post sampling", "C1", "C2", "Presence of mutation")
all.comp.fdr=all.comp.p
all.comp.rho=all.comp.p

abcd=9

pre_factor=WHYF[,abcd]
X=X_1[-which(is.na(pre_factor)),]
pre_factor_1=pre_factor[-which(is.na(pre_factor))]
Y=pre_factor_1

#pre_factor=WHYF[,abcd]
#pre_factor_1=pre_factor
#Y=pre_factor_1
#X=X_1

comp.test=function(x){
  cor.test(X[,x], Y, method = "spearman")$p.value
}

all.cor=sapply(1:27, comp.test)
all.comp.p[abcd,]=all.cor
jason=p.adjust(all.cor, method = "fdr")
all.comp.fdr[abcd,]=jason

comp.test=function(x){
  cor.test(X[,x], Y, method = "spearman")$estimate
}

all.cor=sapply(1:27, comp.test)

all.comp.rho[abcd,]=all.cor

which(jason<0.05)

test4=all.comp.rho
test5=all.comp.p
test6=all.comp.fdr

write.csv(file="DELTA_Spearman_rho_TPP2.csv", test4)
write.csv(file="DELTA_Spearman_fdr_TPP2.csv", test6)
