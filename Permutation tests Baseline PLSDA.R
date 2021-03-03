
Y=as.factor(Y)
levels(Y)=c("1", "2")

say=function(x) {
  {
    Y=sample(Y, replace = F)
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
        test <- samp12
        train <- setdiff(1:nrow(X), test)
        plsda.train <- plsda(X[train, ], Y[train], ncomp = 10, scale = F)
        perf.plsda.train <- perf(plsda.train, validation = "loo", folds = 10, auc = F, progressBar = T, nrepeat = 25, dist = "mahalanobis.dist")
        plsda.train <- plsda(X[train, ], Y[train], ncomp = min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1]))), scale = F)
        test.predict <- predict(plsda.train, X[test, ], dist = "mahalanobis.dist")     #change distance
        Prediction <- test.predict$class$mahalanobis.dist[, min(which(perf.plsda.train$error.rate$BER[,1]==min(perf.plsda.train$error.rate$BER[,1])))]                         #number of components
        well=Y[test]==Prediction
      }
      return(well)
    }
    
    PREvPOST=lapply(1:(as.numeric(nrow(X))), hay)          #number of repeats
    #PREvPOST=lapply(1:2, hay)
    prevpost=unlist(PREvPOST)
    sum(prevpost==F)
    Acc=(sum(prevpost==T)/length(prevpost==F))
    Acc2=(length(Y[which(Y[which(prevpost==TRUE)]=="1")]) + length(Y[which(Y[which(prevpost==TRUE)]=="2")])) / (length(Y[which(Y=="1")]) + length(Y[which(Y=="2")]))
    Bacc=((length(Y[which(Y[which(prevpost==TRUE)]=="1")]) / length(Y[which(Y=="1")])) + (length(Y[which(Y[which(prevpost==TRUE)]=="2")]) / length(Y[which(Y=="2")])))/2
    wellrich=list(PREvPOST, Acc, Bacc)
    #wellrich=Y
  }
  return(wellrich)
}


system.time( par.output <- mclapply.hack( 1:1000,
                                          way)) 

each_bacc=function(x){
  par.output[[x]][[3]]
}
all_bacc=lapply(1:1000, each_bacc)
totNMC=unlist(all_bacc)


#save.image("C:/Users/manz072108/Downloads/Permute_Res_PVP_PQN0_ML_complete.RData")