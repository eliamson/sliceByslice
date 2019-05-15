library("ape")# read.nexus
library(geiger) # for name.check
library("phytools") # contMap

#Slice by slice acquisition of global compactness and cross-sectional area.
#V0.1
#See accompanying ImageJ macro.
#Eli Amson. eli.amson1988@gmail.com


#### 1. Read and clean each profile ####
#Function to read a result .csv file from the accompanying ImageJ macro, and computes the mean of each parameter,
#either the raw measurements ("ResC_NoP" and "ResArea_NoP"), only excluding manually-excluded slices ("ResC" and "ResArea"),
#applying 'local outlier' ("ResCcomp"and "ResAreaComp"), or 'loess regression' ("ResCloess" and "ResAreaLoess") correction procedure.

setwd("~/Documents/utilitaire/R/MammalBoneStruc/MSC-Full-lumbar-allRough")
DataMSCF_List <- list.files(full.names = F,recursive=F)

vertCplot<- function(rowI,DoPlots = FALSE){ 
  #rowI=2; DoPlots=T
  Range <- 4  # for outlier detection
  dfFull <- data.frame(read.csv(DataMSCF_List[rowI],header = T,fileEncoding="latin1"))
  names(dfFull)[names(dfFull)=="ResArea..mm."] <- "ResArea"
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

  #### 2.1. Accurate data acquisition ####  
setwd("~/Documents/utilitaire/R/MammalBoneStruc/MSC-Full-lumbar-allNew")
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
  #### 2.2. Rough data acquisition ####  
setwd("~/Documents/utilitaire/R/MammalBoneStruc/MSC-Full-lumbar-allRough")
for (i in 1:length(DataMSCF_List)){
  #i=1
  df_profile <- vertCplot(i,DoPlots = T)
  LabelVertTemp <- gsub(".csv","",DataMSCF_List[i])
  DataComb[i,"MeanAreaRawRough"] <- mean(df_profile$ResArea)
  DataComb[i,"MeanAreaNoPRough"] <- mean(df_profile$ResArea_NoP)
  DataComb[i,"MeanAreaRough"] <- mean(df_profile$ResAreaComp)
  DataComb[i,"MeanA_LoessRough"] <- mean(df_profile$ResAreaLoess)
  DataComb[i,"MeanCrawRough"] <- mean(df_profile$ResC)
  DataComb[i,"MeanCRough"] <- mean(df_profile$ResCcomp)
  DataComb[i,"MeanC_LoessRough"] <- mean(df_profile$ResCloess)
  DataComb[i,"MeanCNoPRough"] <- mean(df_profile$ResC_NoP)
  
}

#### 3.  Descriptive stats ####
  #### 3.1.  Compactness ####
    #Extremes
DataComb$LabelVertRaw[DataComb$MeanC==max(DataComb$MeanC,na.rm = T)];max(DataComb$MeanC,na.rm = T); 
vertCplot(20,DoPlots = T)
DataComb$LabelVertRaw[DataComb$MeanC==min(DataComb$MeanC,na.rm = T)];min(DataComb$MeanC,na.rm = T)
vertCplot(52,DoPlots = T)
    #Mean of whole dataset
mean(DataComb$MeanC,na.rm = T)
vertCplot(21,DoPlots = T)

  #### 3.2.  Area ####
    #Extremes
DataComb$LabelVertRaw[DataComb$MeanArea==max(DataComb$MeanArea,na.rm = T)];max(DataComb$MeanArea,na.rm = T); 
vertCplot(20,DoPlots = T)
DataComb$LabelVertRaw[DataComb$MeanArea==min(DataComb$MeanArea,na.rm = T)];min(DataComb$MeanArea,na.rm = T)
vertCplot(53,DoPlots = T)
    #Mean of whole dataset
mean(DataComb$MeanArea,na.rm = T)
vertCplot(63,DoPlots = T)


  #### 4. Difference between exclusion methods ####
    #### 4.1. Using accurate acquisition 
DataComb$CleanDiff <- (abs(DataComb$MeanC-DataComb$MeanC_Loess)/DataComb$MeanC)*100
TotalMeanCspan <- max(DataComb$MeanC, na.rm = T)-min(DataComb$MeanC, na.rm = T)
DataComb$CleanDiffToOverallVar <- (abs(DataComb$MeanC-DataComb$MeanC_Loess)/TotalMeanCspan)*100
mean(DataComb$CleanDiffToOverallVar,na.rm = T)
max(DataComb$CleanDiffToOverallVar,na.rm = T)
#=>The greatest difference between MeanC ('semi-automatic'outlier procedure') and MeanC_Loess ('loess regression procedure') represents < 1% of variation described by whole dataset. Mean difference <0.2%

  #### 4.2. Using "rough" acquisition ####
DataComb$DiffMeanCRough <- (abs(DataComb$MeanC-DataComb$MeanCRough)/DataComb$MeanC)*100
DataComb$DiffLoessCRough <- (abs(DataComb$MeanC-DataComb$MeanC_LoessRough)/DataComb$MeanC)*100
mean(DataComb$DiffMeanCRough)
mean(DataComb$DiffLoessCRough)
#=>with MeanC as reference, 0.7% difference between methods, MeanC better
DataComb$DiffMeanCRough <- (abs(DataComb$MeanC_Loess-DataComb$MeanCRough)/DataComb$MeanC)*100
DataComb$DiffLoessCRough <- (abs(DataComb$MeanC_Loess-DataComb$MeanC_LoessRough)/DataComb$MeanC)*100
mean(DataComb$DiffMeanCRough)
mean(DataComb$DiffLoessCRough)
#=>with MeanC_Loess as reference, 0.6% difference between methods, MeanC still better

#### 5. Phylo mapping ####
setwd("~/Documents/utilitaire/R/MammalBoneStruc")
KeyForSp <- data.frame(read.csv("KeyForSp.csv",header = T,fileEncoding="latin1",sep = ";"))
DataComb4tree <- DataComb[DataComb$LabelVertRaw!="ResultsMSC-Petaurista_sp_ZMB_MAM_78560_lumbar-reorVG-InvImJ-An-N",]
DataComb4tree$Taxon <- KeyForSp$Taxon
#DataComb4tree$LabelVert <- KeyForSp$LabelVert
  #TRee
treeSp <- read.nexus("timetreeMammalia_speciesFigtreeExp.nex")$mod1
    #correcting tree labels
treeSp$tip.label[treeSp$tip.label=="Aotus_azarai"] <- "Aotus_azarae"
treeSp$tip.label[treeSp$tip.label=="'Tatera_sp._KIK1704'"] <- "Tatera_indica"
    #swapping species
treeSp$tip.label[treeSp$tip.label=="Rhizomys_pruinosus"] <- "Rhizomys_sumatrensis"
    #rename species (should be valid but not found in TOL)
treeSp$tip.label[treeSp$tip.label=="Proechimys_hoplomyoides"] <- "Proechimys_mincae"

    #Prunning
treeDataComb <- drop.tip(treeSp,drop.tip(treeSp,as.character(DataComb4tree$Taxon))$tip)
#OLD: 14 species are badly spelled???!: DataComb4tree$Taxon[is.na(match(as.character(DataComb4tree$Taxon),treeSp$tip.label))]

  # MeanC  
MeanC <- DataComb4tree$MeanC
names(MeanC)<-DataComb4tree$Taxon
name.check(treeDataComb, MeanC)
#dev.off()
#pdf(file ='Rplot-Kmap.pdf',useDingbats=FALSE,bg="white",width= 4.99,height =6.50)
par(oma=c(0,0,1,0),cex=0.5)
MeanCmap <- contMap(treeDataComb,MeanC,plot=FALSE)
plot(setMap(MeanCmap,invert=TRUE,show.tip.label=F),lwd=3.5,outline=T
     #,type="fan"
     )
title("Mean global compactness", outer=T,cex=12)
#dev.off()

#library("strap") #for scale of timetree (geoscalePhylo function)
#treeDataComb$root.time<-max(nodeHeights(treeDataComb))
#geoscalePhylo(treeDataComb, cex.ts=0.6, cex.tip=0.2,cex.age=0.6, direction="upwards",width=2,label.offset=0.1)

#### 6. Table ####
  #Just MeanC
TabRes<-data.frame()
TabRes["Global compactness","Mean"] <- mean(DataComb$MeanC)
TabRes["Global compactness","SD"] <- sd(DataComb$MeanC)
TabRes["Global compactness","Min"] <- min(DataComb$MeanC)
TabRes["Global compactness","Max"] <- max(DataComb$MeanC)

TabRes["X-sectional area","Mean"] <- mean(DataComb$MeanArea)
TabRes["X-sectional area","SD"] <- sd(DataComb$MeanArea)
TabRes["X-sectional area","Min"] <- min(DataComb$MeanArea)
TabRes["X-sectional area","Max"] <- max(DataComb$MeanArea)

write.table(round(TabRes,digits = 2),file="~/Documents/utilitaire/R/MammalBoneStruc/WholeVertCompRes/TabRes.xls",col.names=NA,sep="	")


  #Comparison between MeanC (local outliers) and MeanC_Loess [MGC=Mean global compactness; MXA=mean cross-sectional area]
TabResC<-data.frame()
TabResC["All slices include","MGC"] <-  paste(round(mean(DataComb$MeanCraw),2)," (", round(sd(DataComb$MeanCraw),2),")",sep = "")
TabResC["All slices include","MXA"] <-  paste(round(mean(DataComb$MeanAreaRaw),2)," (", round(sd(DataComb$MeanAreaRaw),2),")",sep = "")
TabResC["'Accurate' data acquisition",] <- NA
TabResC["No correction1","MGC"] <- paste(round(mean(DataComb$MeanCNoP),2)," (", round(sd(DataComb$MeanCNoP),2),")",sep = "")
TabResC["No correction1","MXA"] <- paste(round(mean(DataComb$MeanAreaNoP),2)," (", round(sd(DataComb$MeanAreaNoP),2),")",sep = "")
TabResC["Local outliers1","MGC"] <- paste(round(mean(DataComb$MeanC),2)," (", round(sd(DataComb$MeanC),2),")",sep = "")
TabResC["Local outliers1","MXA"] <- paste(round(mean(DataComb$MeanArea),2)," (", round(sd(DataComb$MeanArea),2),")",sep = "")
TabResC["Loess1","MGC"] <- paste(round(mean(DataComb$MeanC_Loess),2)," (", round(sd(DataComb$MeanC_Loess),2),")",sep = "")
TabResC["Loess1","MXA"] <- paste(round(mean(DataComb$MeanA_Loess),2)," (", round(sd(DataComb$MeanA_Loess),2),")",sep = "")
TabResC["'Rough' data acquisition",] <- NA
TabResC["No correction2","MGC"] <- paste(round(mean(DataComb$MeanCNoPRough),2)," (", round(sd(DataComb$MeanCNoPRough),2),")",sep = "")
TabResC["No correction2","MXA"] <- paste(round(mean(DataComb$MeanAreaNoPRough),2)," (", round(sd(DataComb$MeanAreaNoPRough),2),")",sep = "")
TabResC["Local outliers2","MGC"] <- paste(round(mean(DataComb$MeanCRough),2)," (", round(sd(DataComb$MeanCRough),2),")",sep = "")
TabResC["Local outliers2","MXA"] <- paste(round(mean(DataComb$MeanAreaRough),2)," (", round(sd(DataComb$MeanAreaRough),2),")",sep = "")
TabResC["Loess2","MGC"] <- paste(round(mean(DataComb$MeanC_LoessRough),2)," (", round(sd(DataComb$MeanC_LoessRough),2),")",sep = "")
TabResC["Loess2","MXA"] <- paste(round(mean(DataComb$MeanA_LoessRough),2)," (", round(sd(DataComb$MeanA_LoessRough),2),")",sep = "")
TabResC
write.table(TabResC,file="~/Documents/utilitaire/R/MammalBoneStruc/WholeVertCompRes/TabResComp.xls",col.names=NA,sep="	")

