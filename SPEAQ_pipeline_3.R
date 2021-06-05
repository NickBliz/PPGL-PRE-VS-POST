library("batman", lib.loc="~/R/win-library/3.4")
brukerdata<-readBruker("C:/Users/z072108/Downloads/U_R")
new_bd=brukerdata[-c(rep(1:31586), rep(61825:62540), rep(63126:68002), rep(71189:72684), rep(73595:74765), rep(75416:75870), rep(88162:89461), rep(90763:91737), rep(96615:131072)),]
spectra=as.matrix(t(new_bd))
colnames(spectra)=spectra[1,]
spectra=spectra[-1,]
garrosh=as.data.frame(row.names(spectra))
duplicated(garrosh)
spectra1=spectra[1:348,]
spectra2=spectra[349:684,]
ppm1=as.numeric(colnames(spectra1))
indx1 <- grepl('QC', row.names(spectra1))
grom1=as.numeric(indx1)
grom1=as.factor(grom1)
ENP1=list("spectra"=spectra1, "ppm"=ppm1, "QCs"=grom1)
ppm2=as.numeric(colnames(spectra2))
indx2 <- grepl('QC', row.names(spectra2))
grom2=as.numeric(indx2)
grom2=as.factor(grom2)
ENP2=list("spectra"=spectra2, "ppm"=ppm2, "QCs"=grom2)
ppm=as.numeric(colnames(spectra))
indx <- grepl('QC', row.names(spectra))
grom=as.numeric(indx)
grom=as.factor(grom)
ENP=list("spectra"=spectra, "ppm"=ppm, "QCs"=grom)

library("speaq", lib.loc="~/R/win-library/3.4")

#peakList1=detectSpecPeaks(ENP1$spectra, nDivRange = c(128), scales = seq(1,16,2), baselineThresh = 0, SNR.Th = -1, verbose = TRUE)
#Y1=dohCluster(ENP1$spectra, peakList = peakList1, refInd = 233, maxShift = 32, acceptLostPeak = TRUE, verbose = TRUE)
Y1.peaks=getWaveletPeaks(Y.spec = ENP1$spectra, X.ppm = ENP1$ppm, baselineThresh = 0, SNR.Th = 3, nCPU = -1, include_nearbyPeaks = TRUE)
Y1.grouped=PeakGrouper(Y.peaks = Y1.peaks, min.samp.grp = 1, grouping.window.width = 200)
Y1.grouped_3=Y1.grouped[-which(Y1.grouped$peakSNR<3),]
Y1.filled=PeakFilling(Y.grouped = Y1.grouped_3, Y.spec=ENP1$spectra, max.index.shift = 8, nCPU = -1)
ENP_features_new1=BuildFeatureMatrix(Y1.filled)
ENP_features_newppm1=BuildFeatureMatrix(Y1.filled, var="peakPPM")
nzmean <- function(x) {
  zvals <- x==0
  if(all(zvals)) 0 else mean(x[!zvals])
}
PEAKSHIFT=ENP_features_newppm1
PEAKSHIFT2=apply(PEAKSHIFT, 2, nzmean)
PEAKSHIFT=rbind(PEAKSHIFT, (as.numeric(nrow(PEAKSHIFT)))+1)
PEAKSHIFT[nrow(PEAKSHIFT),]=PEAKSHIFT2
colnames(ENP_features_new1)=PEAKSHIFT[nrow(PEAKSHIFT),]
row.names(ENP_features_new1)=row.names(spectra1)

#gibbs=ENP_features_new1
#gibbs=gibbs[,-c(61, 62, 63, 64, 65, 66, 67, 68, 114, 115, 116, 154, 155, 156, 157, 158, 159, 257)]

#peakList2=detectSpecPeaks(ENP2$spectra, nDivRange = c(128), scales = seq(1,16,2), baselineThresh = 0, SNR.Th = -1, verbose = TRUE)
#Y2=dohCluster(ENP2$spectra, peakList = peakList2, refInd = 1, maxShift = 8, acceptLostPeak = TRUE, verbose = TRUE)
Y2.peaks=getWaveletPeaks(Y.spec = ENP2$spectra, X.ppm = ENP2$ppm, baselineThresh = 0, SNR.Th = 3, nCPU = -1, include_nearbyPeaks = TRUE)
Y2.grouped=PeakGrouper(Y.peaks = Y2.peaks, min.samp.grp = 1, grouping.window.width = 200)
Y2.grouped_3=Y2.grouped[-which(Y2.grouped$peakSNR<3),]
#Y2.filled=PeakFilling(Y.grouped = Y2.grouped, Y.spec=Y2, max.index.shift = 200, nCPU = -1)
Y22.filled=PeakFilling(Y.grouped = Y2.grouped_3, Y.spec=ENP2$spectra, max.index.shift = 8, nCPU = -1)
#Y2.filled_3=Y2.filled[-which(Y2.filled$peakSNR<3),]
ENP_features_new2=BuildFeatureMatrix(Y22.filled)
ENP_features_newppm2=BuildFeatureMatrix(Y22.filled, var="peakPPM")
nzmean <- function(x) {
  zvals <- x==0
  if(all(zvals)) 0 else mean(x[!zvals])
}
PEAKSHIFT=ENP_features_newppm2
PEAKSHIFT2=apply(PEAKSHIFT, 2, nzmean)
PEAKSHIFT=rbind(PEAKSHIFT, (as.numeric(nrow(PEAKSHIFT)))+1)
PEAKSHIFT[nrow(PEAKSHIFT),]=PEAKSHIFT2
colnames(ENP_features_new2)=PEAKSHIFT[nrow(PEAKSHIFT),]
row.names(ENP_features_new2)=row.names(spectra2)

#jack=ENP_features_new2
#jack=jack[,-c(80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 182, 183, 133, 134, 135, 136, 282)]

Y.peaks=getWaveletPeaks(Y.spec = ENP$spectra, X.ppm = ENP$ppm, baselineThresh = 0, SNR.Th = -1, nCPU = -1, include_nearbyPeaks = TRUE)
Y.grouped=PeakGrouper(Y.peaks = Y.peaks, min.samp.grp = 1, grouping.window.width = 200)
Y.grouped_3=Y.grouped[-which(Y.grouped$peakSNR<3),]
Y.filled=PeakFilling(Y.grouped = Y.grouped_3, Y.spec=ENP$spectra, max.index.shift = 8, nCPU = -1)
ENP_features_new=BuildFeatureMatrix(Y.filled)
ENP_features_newppm=BuildFeatureMatrix(Y.filled, var="peakPPM")
nzmean <- function(x) {
  zvals <- x==0
  if(all(zvals)) 0 else mean(x[!zvals])
}
PEAKSHIFT=ENP_features_newppm
PEAKSHIFT2=apply(PEAKSHIFT, 2, nzmean)
PEAKSHIFT=rbind(PEAKSHIFT, (as.numeric(nrow(PEAKSHIFT)))+1)
PEAKSHIFT[nrow(PEAKSHIFT),]=PEAKSHIFT2
colnames(ENP_features_new)=PEAKSHIFT[nrow(PEAKSHIFT),]
row.names(ENP_features_new)=row.names(spectra)


setwd("C:/Users/z072108/Downloads/U_R")
write.csv(ENP_features_new1, file="ENP_features1_3.csv")
write.csv(ENP_features_new2, file="ENP_features2_3.csv")
