library(xcms)
library(MSnbase)
library(dplyr)

N <- 5 #1to5 represents 5 software
masstol <- 0.01
RTtol <- 30
exportfilename <- "OpenMSdetectedNOTNew.csv"

directory <- "F:/Jian_Guo/SoftwareComparison_20211103/singleDatafilecompare20220207/New10datasetResults20220204/BrukerUrineRPorigianlDDA/ReasonforDetectability20220228"
setwd(directory)
Alignedtable <- read.csv("Aligned.csv", stringsAsFactors = F)
ms1data <- readMSData(files = "UrineOriginal1.mzXML", mode = "onDisk", msLevel. = 1)
mzData <- mz(ms1data)
intData <- intensity(ms1data)
rtime <- rtime(ms1data)
################################################################################################################
Detectable <- Alignedtable[which(Alignedtable[,3*N]== 0), c((3*N-2):(3*N-1))]
peakslope <- c()
peaksharpness <- c()
peakHeight <- c()
peakSNratio <- c()
peakscannumber <- c()
peakwidth <- c()
peakmassaccuracy <- c()
peakmassaccuracyDa <- c()
for(i in 1:nrow(Detectable)) {
  currmzRange <- c(Detectable[i,1]-masstol, Detectable[i,1]+masstol)
  tmpMZdata <- mzData
  tmpINTdata <- intData
  for(j in 1:length(mzData)){
    index <- which(tmpMZdata[[j]] >= currmzRange[1] & tmpMZdata[[j]] < currmzRange[2])
    tmpMZdata[[j]] <- tmpMZdata[[j]][index]
    tmpINTdata[[j]] <- tmpINTdata[[j]][index]
  }
  # Extract the intensity vectors from each m/z bin 
  eicINT <- c()
  eicRT <- c()
  for(k in 1:length(mzData)){
    if(length(tmpINTdata[[k]]) > 0){
      eicINT[k] <- mean(tmpINTdata[[k]])
    }else{
      eicINT[k] <- 0
    }
    eicRT[k] <- rtime[k]
  }
  if(sum(eicINT != 0) == 0) next()
  # Sort the intensity vectors from each m/z bin, estimate the noise cut off and average
  eicNon0 <- sort(eicINT[eicINT > 0])
  if(length(eicNon0) > 10){
    for(x in seq(10,length(eicNon0), 10)){
      sd <- sd(eicNon0[1:x])
      blk <- sum(eicNon0[1:x])/x
      thres <- blk + 3*sd
      if(x+1 <= length(eicNon0)){
        if(eicNon0[x+1] >= thres) break()
      }
    }
    cutOFF <- eicNon0[x]
  }else{
    cutOFF <- max(eicNon0)
    sd <- 0
    blk <- 0
  }
  aboveTHindex <- which(eicINT > cutOFF & eicRT <= Detectable[i,2]+RTtol & eicRT >= Detectable[i,2]-RTtol)
  if(length(aboveTHindex) == 0) next()
  peakInd <- aboveTHindex[which(eicINT[aboveTHindex] == max(eicINT[aboveTHindex]))[1]]
  highestINT <- which(tmpINTdata[[peakInd]] == max(tmpINTdata[[peakInd]]))[1]
  refMZvec <- tmpMZdata[[peakInd]][highestINT]
  currSamePeakMass <- c()
  currSamePeakMass <- c(currSamePeakMass, refMZvec)
  LeftPeakslope <- c()
  RightPeakslope <- c()
  Denominatorslope <- c()
  LeftPeaksharp <- c()
  RightPeaksharp <- c()
  leftInd <- peakInd-1
  rightInd <- peakInd+1
  if(leftInd >= min(aboveTHindex)){
    while (length(tmpMZdata[[leftInd]]) > 0 & mean(tmpINTdata[[leftInd]]) >= cutOFF) {
      if (length(tmpMZdata[[leftInd]]) == 1){
        currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[leftInd]])
        if(eicINT[leftInd] > eicINT[leftInd+1] & length(currSamePeakMass) > 5){
          Q1 <- as.numeric(summary(currSamePeakMass)[2])
          Q3 <- as.numeric(summary(currSamePeakMass)[5])
          LB <- Q1 - 1.5 *(Q3 - Q1)
          RB <- Q3 + 1.5 *(Q3 - Q1)
          if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
        }
      } else {
        abvector <- abs(tmpMZdata[[leftInd]] - refMZvec)
        NearInd <- which(abvector == min(abvector))[1]
        currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[leftInd]][NearInd])
        if(eicINT[leftInd] > eicINT[leftInd+1] & length(currSamePeakMass) > 5){
          Q1 <- as.numeric(summary(currSamePeakMass)[2])
          Q3 <- as.numeric(summary(currSamePeakMass)[5])
          LB <- Q1 - 1.5 *(Q3 - Q1)
          RB <- Q3 + 1.5 *(Q3 - Q1)
          if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
        }
      }
      if(eicINT[peakInd] == 0){
        LeftPeaksharp <- c(LeftPeaksharp, 0) 
      }
      if(leftInd > 0 & leftInd < length(tmpINTdata) & eicINT[peakInd] != 0) {
        LeftPeaksharp <- c(LeftPeaksharp, abs(eicINT[peakInd] - eicINT[leftInd])/(abs(peakInd - leftInd)*sqrt(eicINT[peakInd])))
      }
      
      if(leftInd > 2 & leftInd < length(tmpINTdata) - 2) {
        Firstderivative <- (-2 * eicINT[leftInd-2] - eicINT[leftInd-1] + eicINT[leftInd + 1] + 2 * eicINT[leftInd + 2])/10
      } else {
        Firstderivative <- 0
      }
      if (Firstderivative >= 0) {
        LeftPeakslope <- c(LeftPeakslope, Firstderivative)
      } 
      Denominatorslope <- c(Denominatorslope, abs(Firstderivative))
      
      leftInd <- leftInd-1
      if(leftInd <= 0) break()
    }
  }
  if(rightInd <= max(aboveTHindex)){
    while (length(tmpMZdata[[rightInd]]) > 0 & mean(tmpINTdata[[rightInd]]) >= cutOFF) {
      if (length(tmpMZdata[[rightInd]]) == 1){
        currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[rightInd]])
        if(eicINT[rightInd] > eicINT[rightInd-1] & length(currSamePeakMass) > 5){
          Q1 <- as.numeric(summary(currSamePeakMass)[2])
          Q3 <- as.numeric(summary(currSamePeakMass)[5])
          LB <- Q1 - 1.5 *(Q3 - Q1)
          RB <- Q3 + 1.5 *(Q3 - Q1)
          if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
        }
      } else {
        abvector <- abs(tmpMZdata[[rightInd]] - refMZvec)
        NearInd <- which(abvector == min(abvector))[1]
        currSamePeakMass <- c(currSamePeakMass, tmpMZdata[[rightInd]][NearInd])
        if(eicINT[rightInd] > eicINT[rightInd-1] & length(currSamePeakMass) > 5){
          Q1 <- as.numeric(summary(currSamePeakMass)[2])
          Q3 <- as.numeric(summary(currSamePeakMass)[5])
          LB <- Q1 - 1.5 *(Q3 - Q1)
          RB <- Q3 + 1.5 *(Q3 - Q1)
          if (currSamePeakMass[length(currSamePeakMass)] < LB || currSamePeakMass[length(currSamePeakMass)] > RB) break()
        }
      }
      
      if(eicINT[peakInd] == 0) {
        RightPeaksharp <- c(RightPeaksharp, 0)
      }
      if(rightInd > 0 & rightInd < length(tmpINTdata) & eicINT[peakInd] != 0) {
        RightPeaksharp <- c(RightPeaksharp, abs(eicINT[peakInd] - eicINT[rightInd])/(abs(peakInd - rightInd)*sqrt(eicINT[peakInd])))
      } 
      if(rightInd > 2 & rightInd < length(tmpINTdata) - 2) {
        Firstderivative <- (-2 * eicINT[rightInd-2] - eicINT[rightInd-1] + eicINT[rightInd + 1] + 2 * eicINT[rightInd + 2])/10
      } else {
        Firstderivative <- 0
      }
      if (Firstderivative <= 0) {
        RightPeakslope <- c(RightPeakslope, abs(Firstderivative))
      } 
      Denominatorslope <- c(Denominatorslope, abs(Firstderivative))
      
      rightInd <- rightInd+1
      if(rightInd > length(tmpMZdata)) break()
    }
  }
  if(sum(is.na(currSamePeakMass)) > 0) next()
  if(length(currSamePeakMass) > 1 && length(currSamePeakMass) < 200){
    peakmassaccuracy <- c(peakmassaccuracy, (2*sd(currSamePeakMass))/refMZvec * 1e6)
    peakmassaccuracyDa <- c(peakmassaccuracyDa, 2*sd(currSamePeakMass))
    peakwidth <- c(peakwidth, rtime[[rightInd - 1]] - rtime[[leftInd + 1]])
    peakscannumber <- c(peakscannumber, rightInd - leftInd - 1)
    if (cutOFF > 0){
      peakSNratio <- c(peakSNratio, abs(eicINT[peakInd]-blk)/sd)
    } else {
      peakSNratio <- c(peakSNratio, eicINT[peakInd]) 
    }
    peakHeight <- c(peakHeight, eicINT[peakInd])
    slope <- (sum(LeftPeakslope) + sum(RightPeakslope))/sum(Denominatorslope)
    peakslope <- c(peakslope, slope)
    if(length(LeftPeaksharp) > 0 & length(RightPeaksharp) > 0){
      sharpness <- (max(LeftPeaksharp) + max(RightPeaksharp))/2
    }
    if(length(LeftPeaksharp) > 0 & length(RightPeaksharp) == 0){
      sharpness <- max(LeftPeaksharp)/2
    }
    if(length(LeftPeaksharp) == 0 & length(RightPeaksharp) > 0){
      sharpness <- max(RightPeaksharp)/2
    }
    if(length(LeftPeaksharp) == 0 & length(RightPeaksharp) == 0){
      sharpness <- 0
    }
    peaksharpness <- c(peaksharpness, sharpness)
  }
}

denspeakslope <- density(peakslope, from = 0, to = 1, bw = 0.001)
denspeaksharpness <- density(peaksharpness, from = 0)
denspeakHeight <- density(peakHeight, from = 0)
denspeakSNratio <- density(peakSNratio, from = 0)
denspeakscannumber <- density(peakscannumber, from = 0)
denspeakwidth <- density(peakwidth, from = 0)
denspeakmassaccuracy <- density(peakmassaccuracy, from = 0)
denspeakmassaccuracyDa <- density(peakmassaccuracyDa, from = 0)
##################################################################################################


Densitytable <- data.frame(matrix(nrow = 512, ncol = 16))
colnames(Densitytable) <- c("slopeX","slopeY","sharpX","sharpY","HeightX","HeightY","SNratioX","SNratioY",
                            "scanX","scanY","widthX","widthY","ppmX","ppmY", "DaX", "DaY")
Densitytable[,1] <- denspeakslope[["x"]]
Densitytable[,2] <- denspeakslope[["y"]]
Densitytable[,3] <- denspeaksharpness[["x"]]
Densitytable[,4] <- denspeaksharpness[["y"]]
Densitytable[,5] <- denspeakHeight[["x"]]
Densitytable[,6] <- denspeakHeight[["y"]]
Densitytable[,7] <- denspeakSNratio[["x"]]
Densitytable[,8] <- denspeakSNratio[["y"]]
Densitytable[,9] <- denspeakscannumber[["x"]]
Densitytable[,10] <- denspeakscannumber[["y"]]
Densitytable[,11] <- denspeakwidth[["x"]]
Densitytable[,12] <- denspeakwidth[["y"]]
Densitytable[,13] <- denspeakmassaccuracy[["x"]]
Densitytable[,14] <- denspeakmassaccuracy[["y"]]
Densitytable[,15] <- denspeakmassaccuracyDa[["x"]]
Densitytable[,16] <- denspeakmassaccuracyDa[["y"]]
setwd(directory)
write.csv(Densitytable, file = exportfilename, row.names = FALSE)

