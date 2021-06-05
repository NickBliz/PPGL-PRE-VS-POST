library(readr)
ENU1 <- read_csv("C:/Users/z072108/Downloads/U_R/ENP_features1_3.csv")     #open datasets
ENU2 <- read_csv("C:/Users/z072108/Downloads/U_R/ENP_features2_3.csv")

setwd("C:/Users/z072108/Downloads/U_R")
ENU1.d=as.data.frame(ENU1)          #change name                                   #convert to data frames
ENU2.d=as.data.frame(ENU2)          #change name
ENU1.d=ENU1.d[,-1]                                                     #remove 1st column with sample names
ENU2.d=ENU2.d[,-1]
ENU1.d$names=c(rep(1:as.numeric(nrow(ENU1.d))))                        #add a column with a number for each sample
ENU2.d$names=c(rep((as.numeric(nrow(ENU1.d))+1):(as.numeric(nrow(ENU1.d))+as.numeric(nrow(ENU2.d)))))
library("data.table", lib.loc="~/R/win-library/3.4")
library("reshape2", lib.loc="~/R/win-library/3.4")
ENU1.dm=melt(ENU1.d, id="names", na.rm = T)                            # melt data frames
ENU2.dm=melt(ENU2.d, id="names", na.rm = T)
ENU1.2dm=rbind(ENU1.dm, ENU2.dm)                                       #bind melted data frames
ENU1.2dmdc=dcast(ENU1.2dm, names~variable, value.var = 'value')        #un-melt combined data frame
#ENU1.2dmdc=ENU1.2dmdc[, !grepl("_", colnames(ENU1.2dmdc))]             #remove columns with 0_
ENU1.2dmdc=ENU1.2dmdc[, !colnames(ENU1.2dmdc)<0.5]                     #remove columns with ppm values < 0.5
UR1=ENU1.2dmdc[1:as.numeric(nrow(ENU1.d)),]
UR2=ENU1.2dmdc[((as.numeric(nrow(ENU1.d))+1):(as.numeric(nrow(ENU1.d))+as.numeric(nrow(ENU2.d)))),]
UR1=UR1[, order(colnames(UR1))]
UR2=UR2[, order(colnames(UR2))]

UR1.1m=UR1
UR1.1m=rbind(UR1.1m, (as.numeric(nrow(UR1.1m))+1))
checkUR1=function(x) {
  abs(as.numeric(colnames(UR1.1m)[x])-as.numeric(colnames(UR1.1m)[x+1]))
}
UR1.1mm=sapply(c(1:(as.numeric(ncol(UR1.1m))-1)), checkUR1)
UR1.1m[as.numeric(nrow(UR1.1m)), ]=cbind(UR1.1mm, as.numeric(ncol(UR1.1m)))
UR1.1m=rbind(UR1.1m, (as.numeric(nrow(UR1.1m))+1))
UR1.1m[as.numeric(nrow(UR1.1m)),]= UR1.1m[(as.numeric(nrow(UR1.1m))-1),]
for (i in seq_along(UR1.1m)) {
  if((((abs(as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i+1])))-(abs(as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i+2])))>0) | ((abs(as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i+1])))-(abs(as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i])))>0)) | (((colSums(is.na(UR1.1m[i+2]))>((as.numeric(nrow(UR1)))-1)) & (colSums(is.na(UR1.1m[i+1]))>((as.numeric(nrow(UR1)))-1))) | ((colSums(!is.na(UR1.1m[i+2]))>((as.numeric(nrow(UR1)))-1)) & (colSums(!is.na(UR1.1m[i+1]))>((as.numeric(nrow(UR1)))-1))))) {UR1.1m[(as.numeric(nrow(UR1.1m))),i+1]=0}
}
UR1.1m=rbind(UR1.1m, (as.numeric(nrow(UR1.1m))+1))
UR1.1m[as.numeric(nrow(UR1.1m)),]= UR1.1m[(as.numeric(nrow(UR1.1m))-1),]
for (i in seq_along(UR1.1m)) {
  if((((colSums(is.na(UR1.1m[i+2]))>((as.numeric(nrow(UR1)))-1)) & (colSums(!is.na(UR1.1m[i+1]))>((as.numeric(nrow(UR1)))-1))) | ((colSums(!is.na(UR1.1m[i+2]))>((as.numeric(nrow(UR1)))-1)) & (colSums(is.na(UR1.1m[i+1]))>((as.numeric(nrow(UR1)))-1)))) & (as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i+1]==0)) & (as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i]==0)) & (as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i+2]==0))) {UR1.1m[(as.numeric(nrow(UR1.1m))),i+1]=(as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-2),i+1]))}
}
UR1.1m=rbind(UR1.1m, (as.numeric(nrow(UR1.1m))+1))
UR1.1m[as.numeric(nrow(UR1.1m)),]= UR1.1m[(as.numeric(nrow(UR1.1m))-1),]
for (i in seq_along(UR1.1m)) {
  if((as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i]!=0)) & (as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i+1]!=0)) & (as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i+1])>as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i]))) {UR1.1m[(as.numeric(nrow(UR1.1m))),i+1]=0}
}
for (i in seq_along(UR1.1m)) {
  if((as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i]!=0)) & (as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i+1]!=0)) & (as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i])>as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))-1),i+1]))) {UR1.1m[(as.numeric(nrow(UR1.1m))),i]=0}
}
pam=UR1.1m[,-as.numeric(ncol(UR1.1m))]
bay=max(pam[(as.numeric(nrow(UR1.1m))-1),], na.rm=T)
for (i in seq_along(UR1.1m)) {
  if((as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))),i])>=as.numeric(bay))) {UR1.1m[(as.numeric(nrow(UR1.1m))),i]=0}
}
for (i in seq_along(colnames(UR1.1m))) {
  if(as.numeric(UR1.1m[(as.numeric(nrow(UR1.1m))),i])>0) {colnames(UR1.1m)[i]=colnames(UR1.1m)[i+1]}
}
UR1.1=UR1.1m
UR1.1=UR1.1[-c((as.numeric(nrow(UR1.1))-1),(as.numeric(nrow(UR1.1))), (as.numeric(nrow(UR1.1))-3), (as.numeric(nrow(UR1.1))-2)),]
UR1.1=rbind(UR1.1, (as.numeric(nrow(UR1.1))+1))
UR1.1[(as.numeric(nrow(UR1.1))),]=colnames(UR1.1)
UR1.1=UR1.1[,-(as.numeric(ncol(UR1.1)))]
UR1.1=UR1.1[,order(UR1.1[1,])]
colnames(UR1.1)=UR1.1[(as.numeric(nrow(UR1.1))),]
UR1.5=UR1.1[(!duplicated(colnames(UR1.1), fromLast=F))]
UR1.5=UR1.5[-(as.numeric(nrow(UR1.1))),]
row.names(UR1.5)=ENU1$X1                  #change name
UR1.6=UR1.5[,order(colnames(UR1.5))]

UR2.1=UR2
colnames(UR2.1)=colnames(UR1.1m)
UR2.1.1=UR2.1
UR2.1.1=rbind(UR2.1.1, (as.numeric(nrow(UR2.1.1))+1))
UR2.1.1[as.numeric(nrow(UR2.1.1)),]=colnames(UR2.1.1)
UR2.1.1=UR2.1.1[,-as.numeric(ncol(UR2.1.1))]
UR2.1.1=UR2.1.1[,order(UR2.1.1[1,])]
colnames(UR2.1.1)=UR2.1.1[as.numeric(nrow(UR2.1.1)),]
UR2.1.5=UR2.1.1[(!duplicated(colnames(UR2.1.1), fromLast=F))]
UR2.1.5=UR2.1.5[-as.numeric(nrow(UR2.1.5)),]
ENU2$X1[299]="QC 3 RUN 45"
ENU2$X1[300]="QC 3 RUN 45b"
ENU2$X1[314]="QC 3 RUN 46"
row.names(UR2.1.5)=ENU2$X1                #change name
UR2.1.6=UR2.1.5[,order(colnames(UR2.1.5))]
UR1_2=rbind(UR1.6, UR2.1.6)               #change name
will=UR1_2                                #change name
for (i in seq_along(will)) {
  if((((colSums(is.na(will[1:as.numeric(nrow(UR1)),][i]))>((as.numeric(nrow(UR1))-1))) & (colSums(is.na(will[(as.numeric(nrow(UR1))+1):as.numeric(nrow(ENU1.2dmdc)),][i+1]))>((as.numeric(nrow(UR2))-1)))) | ((colSums(is.na(will[1:as.numeric(nrow(UR1)),][i+1]))>((as.numeric(nrow(UR1))-1))) & (colSums(is.na(will[(as.numeric(nrow(UR1))+1):as.numeric(nrow(ENU1.2dmdc)),][i]))>((as.numeric(nrow(UR2))-1)))))) {colnames(will)[i]=colnames(will)[i+1]}
}
bill=cbind(will, ncol(will)+1)
bill=rbind(bill, nrow(bill)+1)
bill[nrow(bill),]=colnames(bill)
bill[,ncol(bill)]=c(1:(as.numeric(nrow(ENU1.2dmdc))+1))
bill=bill[-nrow(bill),]
colnames(bill)[colnames(bill)=="ncol(will) + 1"] = "names"
davy=bill[(!duplicated(colnames(bill), fromLast=T))]
jones=bill[(duplicated(colnames(bill), fromLast=T))]
jones=cbind(jones, ncol(jones)+1)
jones[,ncol(jones)]=davy$names
colnames(jones)[colnames(jones)=="ncol(jones) + 1"] = "names"
davym=melt(davy, id="names", na.rm = T)
jonesm=melt(jones, id="names", na.rm = T)
bootstrapm=rbind(davym, jonesm)
bootstrap=dcast(bootstrapm, names~variable, value.var = 'value')
row.names(bootstrap)=row.names(UR1_2)     #change name
bootstrap=bootstrap[,-1]
bootstrap[is.na(bootstrap)]=0
#bootstrap=bootstrap[, -c(141, 53, 52, 59, 64, 65, 82, 104, 119, 73, 74, 77, 76, 75, 103)]   #remove blank peaks and histidine:3.346, 2.056, 1.904, 2.638-2.6886, 2.50-5.545, 7.055-7.0758, 7.79-7.868
write.csv(file="ENP_final_3_new.csv", bootstrap)   #change name

#mitch=bootstrap
#mitch[mitch==0]<-NA
#mitch=mitch[,-which(colSums(is.na(mitch[1:348,]))==348)]
#mitch=mitch[,-which(colSums(is.na(mitch[349:684,]))==336)]

