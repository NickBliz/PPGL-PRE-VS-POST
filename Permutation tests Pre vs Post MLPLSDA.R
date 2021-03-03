ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG=PVP
x <- ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG[,-(1:37)]
y1=as.factor(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$GROUP)
X <- x
enp_names=ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG$`ENSAT-HT ID`
row.names(X)=make.unique(enp_names)
Y <- as.factor(y1)
levels(Y)=c("1", "2")

say=function(x) {
  {
    Y=sample(Y, replace = F)
    Y1=as.integer(Y[1:(as.numeric(nrow(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG))/2)])
    Y2=Y1
    Y2[Y2==1]=3
    Y2[Y2==2]=1
    Y2[Y2==3]=2
    Y=c(Y1, Y2)
    Y=as.factor(Y)
  }
}

set.seed(1)
newY=lapply(1:1000, say)

way=function(y) {
  {
    Y=newY[[y]]
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
      }
      return(well)
    }
    
    PREvPOST=lapply((1:(as.numeric(nrow(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG))/2)), hay)          #number of repeats
    prevpost=unlist(PREvPOST)
    post=prevpost[seq(1,as.numeric(nrow(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG)), 2)]
    ppgl=prevpost[seq(2,as.numeric(nrow(ENP_final_3_bothQCexp_16072019_excl_blank_PHT_outliers_forPREvPOST_20_new_MA_PQN_min_0_3_GLOG)), 2)]
    post[which(post==FALSE)]="FN"
    post[which(post==TRUE)]="TN"
    ppgl[which(ppgl==FALSE)]="FP"
    ppgl[which(ppgl==TRUE)]="TP"
    preNpost=cbind(post, ppgl)
  }
  return(preNpost)
}

system.time( par.output <- mclapply.hack( 1:1000,
                                          way)) 

PREvPOSTpermute=par.output


#PREvPOSTpermute=lapply((1:1000), way)
p=do.call(cbind, PREvPOSTpermute)

ray=function(z){
  l=p[,c(z,z+1)]
  FP=length(which(l=="FP"))
  FN=length(which(l=="FN"))
  NMC=FP+FN
  return(NMC)
}

totNMC=sapply(as.numeric(seq(1, ncol(p), 2)), ray)

save.image("C:/Users/manz072108/Downloads/Permute_Res_PVP_PQN0_ML_complete.RData")