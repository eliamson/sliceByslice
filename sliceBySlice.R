#Slice by slice acquisition of global compactness and cross-sectional area.
#V0.3
#See accompanying ImageJ macro.
#Eli Amson. eli.amson1988@gmail.com

library(easycsv)# for choose_dir()

#### 1. Read and clean each profile ####
#Function to read a result .csv file from the accompanying ImageJ macro, and computes the mean of each parameter,
#either the raw measurements ("ResC_NoP" and "ResArea_NoP"), only excluding manually-excluded slices ("ResC" and "ResArea"),
#applying 'local outlier' ("ResCcomp"and "ResAreaComp"), or 'loess regression' ("ResCloess" and "ResAreaLoess") correction procedure.

vertCplot<- function(rowI,DoPlots = FALSE){ 
  Range <- 4  # for outlier detection
  dfFull <- data.frame(read.csv(DataMSCF_List[rowI],header = T,fileEncoding="latin1"))
  levels(dfFull$ToRem) <- c(levels(dfFull$ToRem),"3")
  Z1=dfFull$X[dfFull$ToRem=="Z1"]
  Z2=dfFull$X[dfFull$ToRem=="Z2"]
  
  #Local outlier detector: > Q_3 + 1.5 * IQR of a  -Range-3:+Range+3 subsample (excluding slice of interest and direct neighbours)
  for (j in Z1:Z2){ 
    if(j>Z1+Range && j<Z2-Range && dfFull$ToRem[j]==0){
        subSamp <- dfFull[c((j-Range):(j-2),(j+2):(j+Range)),]
        subSamp <- subset(subSamp,ToRem=="0")$ResC
      if (length(subSamp)!=0) {
          boxplotStats <- boxplot(subSamp,plot = F)$stats
          U_W <- boxplotStats[5, ]+(boxplotStats[4, ]-boxplotStats[2, ])
        if(dfFull$ResC[j]>U_W && dfFull$ToRem[j]==0){
            dfFull$ToRem[j]=3
            j=j+1
          while(dfFull$ResC[j]>U_W && dfFull$ToRem[j]!=1){
              dfFull$ToRem[j]=3
              j=j+1
          }
        } 
      }
    }
  }
    dfInt <- dfFull[Z1:Z2,]
    #No procedrue
    dfInt$ResC_NoP <- dfInt$ResC
    dfInt$ResArea_NoP <- dfInt$ResArea
    #'Outlier procedure"  
    dfInt$ResCcomp <- dfInt$ResC
    dfInt$ResAreaComp <- dfInt$ResArea
    
    #Replace values to exclude by a sequence based on neighbours.
  for (k in 1:nrow(dfInt)){
    if(dfInt$ToRem[k]=="1" && dfInt$ToRem[k-1]!="1"){
        nToRep <- 0
      while (dfInt$ToRem[k]!="0"){
          nToRep <- nToRep+1
          k=k+1
      }
        dfInt$ResCcomp[(k-nToRep):(k-1)] <- seq(from=dfInt$ResCcomp[k-1-nToRep],to=dfInt$ResCcomp[k],length.out=nToRep+2)[2:(nToRep+1)]
        dfInt$ResAreaComp[(k-nToRep):(k-1)] <- seq(from=dfInt$ResAreaComp[k-1-nToRep],to=dfInt$ResAreaComp[k],length.out=nToRep+2)[2:(nToRep+1)]
        
        dfInt$ResC_NoP[(k-nToRep):(k-1)] <- seq(from=dfInt$ResC_NoP[k-1-nToRep],to=dfInt$ResC_NoP[k],length.out=nToRep+2)[2:(nToRep+1)]
        dfInt$ResArea_NoP[(k-nToRep):(k-1)] <- seq(from=dfInt$ResArea_NoP[k-1-nToRep],to=dfInt$ResArea_NoP[k],length.out=nToRep+2)[2:(nToRep+1)]
    }
    if(dfInt$ToRem[k]=="3"){
        nToRep3 <- 0
      while (dfInt$ToRem[k]!="0"){
          nToRep3 <- nToRep3+1
          k=k+1
      }
        dfInt$ResCcomp[(k-nToRep3):(k-1)] <- seq(from=dfInt$ResCcomp[k-nToRep3-1],to=dfInt$ResCcomp[k],length.out=nToRep3+2)[2:(nToRep3+1)]
        dfInt$ResAreaComp[(k-nToRep3):(k-1)] <- seq(from=dfInt$ResAreaComp[k-nToRep3-1],to=dfInt$ResAreaComp[k],length.out=nToRep3+2)[2:(nToRep3+1)]
    } 
  }
    dfInt$ResCcomp <- signif(dfInt$ResCcomp,5); dfInt$ResAreaComp <- signif(dfInt$ResAreaComp,5)
    dfInt$ResC_NoP <- signif(dfInt$ResC_NoP,5); dfInt$ResArea_NoP <- signif(dfInt$ResArea_NoP,5)
    #Plots
  if (DoPlots==T){
      cols1  <- dfInt$ToRem; levels(cols1) <- c(levels(cols1),"black","coral"); cols1[cols1=="Z1"|cols1=="Z2"|cols1=="0"] <- "black"; cols1[cols1=="1"|cols1=="3"] <- "coral"; cols1 <- factor(cols1); cols1 <- c("black","coral")[cols1]
      plot(dfInt$ResCcomp,xaxt = "n",col=cols1)
      axis(1, at=round(seq(Z1, Z2, length.out=10)-Z1), labels=round(seq(Z1, Z2, length.out=10)))
      lines(1:nrow(dfInt), dfInt$ResCcomp)
      title(DataMSCF_List[rowI], outer=F)
  }
    #Comparison (of resuting mean) WITH LOESS (so not using "3")
      # i.e., non-parametric Local Polynomial Regression using a nearest neighbour approach
    dfInt$ResCno1_2 <- dfInt$ResC
    dfInt$ResAreaNo1_2 <- dfInt$ResArea
    for (l in 1:length(dfInt$ResCno1_2)){
      if (dfInt$ToRem[l] ==1)
        dfInt$ResCno1_2[l] <- dfInt$ResCcomp[l];
        dfInt$ResAreaNo1_2[l] <- dfInt$ResAreaComp[l];
  }
    Seq_No1_2 <- 1:length(dfInt$ResCno1_2)
    lo1_2 <- loess(dfInt$ResCno1_2~Seq_No1_2,span=1) #excluding only initial exclusions ("1")
    Arlo1_2 <- loess(dfInt$ResAreaNo1_2~Seq_No1_2,span=1) #excluding only initial exclusions ("1")
    dfInt$ResCloess<- predict(lo1_2) #append C values predicted by LOESS (to new variable)
    dfInt$ResAreaLoess<- predict(Arlo1_2) #append Area values predicted by LOESS (to new variable)
    lines(dfInt$ResCloess, col='red', lwd=2)
    
  # Append results
  return(dfInt)
} 

#### 2. Combine all profiles ####
DataComb <- data.frame()
#Both the 'local outlier' and the 'loess regression' procedures are included, respectively yielding the parameters MeanC and MeanArea, and MeanC_Loess and MeanA_Loess.

#Choose data folder
setwd(choose_dir())
#
DataMSCF_List <- list.files(full.names = F,recursive=F)
for (i in 1:length(DataMSCF_List)){
  df_profile <- vertCplot(i,DoPlots = T)
  LabelVertTemp <- gsub(".csv","",DataMSCF_List[i])
  DataComb[i,"LabelVertRaw"] <- gsub("-An-F","-An",LabelVertTemp)
  DataComb[i,"MeanAreaNoP"] <- mean(df_profile$ResArea_NoP)
  DataComb[i,"MeanAreaRaw"] <- mean(df_profile$ResArea)
  DataComb[i,"MeanArea"] <- mean(df_profile$ResAreaComp)
  DataComb[i,"MeanArea_SD"] <- sd(df_profile$ResAreaComp)
  DataComb[i,"MeanA_Loess"] <- mean(df_profile$ResAreaLoess)
  DataComb[i,"MeanCraw"] <- mean(df_profile$ResC)
  DataComb[i,"MeanCNoP"] <- mean(df_profile$ResC_NoP)
  DataComb[i,"MeanC"] <- mean(df_profile$ResCcomp)
  DataComb[i,"MeanC_SD"] <- sd(df_profile$ResCcomp)
  DataComb[i,"MeanC_Loess"] <- mean(df_profile$ResCloess)
  DataComb[i,"PropRem"] <- (length(which(df_profile$ToRem!=0))+2)/nrow(df_profile)  #Proportion of excluded slices
}
