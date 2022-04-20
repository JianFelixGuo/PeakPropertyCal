###This is a script to identify overlaping MS2 spectras##
#Jian Guo
#2021-05-10
###################################################################
# Parameter setting
library(dplyr)
directory <- "F:/Jian_Guo/SoftwareComparison_20211103/singleDatafilecompare20220207/New10datasetResults20220204/0REDO_20220404/FalseNegRate/AgilentExposomeDDA/FalseNegRate"
# file1name <- "fold4level3.csv"
# file2name <-"fold2level1&2.csv"
mass.tol <- 0.02
rt.tol <- 20
#sample file number
n <- 5
###################################################################
# Feature alignment
setwd(directory)
filename <- list.files(pattern = ".csv")
file1 <- read.csv(filename[1], stringsAsFactors = F)
file2 <- read.csv(filename[2], stringsAsFactors = F)
file1T <- data.frame(matrix(nrow = 0, ncol = 3))
file1T <- file1[, c("mz", "rt", "intensity")]
colnames(file1T)[3] <- "1"
file2T <- data.frame(matrix(nrow = 0, ncol = 3))
file2T <- file2[, c("mz", "rt", "intensity")]
colnames(file2T)[3] <- "2"
output1 <- data.frame(matrix(nrow = 0, ncol = 2*ncol(file1T)))
for(i in 1:nrow(file1T)){
    mass.lower.limit <- as.numeric(file1T$mz[i]) - mass.tol
    mass.upper.limit <- as.numeric(file1T$mz[i]) + mass.tol
    rt.lower.limit <- as.numeric(file1T$rt[i]) - rt.tol
    rt.upper.limit <- as.numeric(file1T$rt[i]) + rt.tol
    temp <- file2T[(as.numeric(file2T$mz) >= as.numeric(mass.lower.limit) &
                     as.numeric(file2T$mz) <= as.numeric(mass.upper.limit) &
                     as.numeric(file2T$rt) >= as.numeric(rt.lower.limit) &
                     as.numeric(file2T$rt) <= as.numeric(rt.upper.limit)),]
  temp <- temp[complete.cases(temp),]
  if(nrow(temp) > 0) {
    tmprow <- cbind(file1T[i,], temp[1,])
    colnames(tmprow) <- colnames(output1)
    output1 <- rbind(output1, tmprow)
  }
}
colnames(output1) <- c(colnames(file1T), colnames(file2T))
colnames(output1)[3:6] <- c("1", "mz 2", "rt 2", "2")
unique1_1 <- anti_join(file1T, as.data.frame(output1[1:ncol(file1T)]))
unique1_2 <- anti_join(file2T, as.data.frame(output1[(ncol(file1T)+1):ncol(output1)]))
if (nrow(unique1_1) > 0 & nrow(unique1_2) > 0){
unique1_1T <- data.frame(matrix(nrow = 0, ncol = 2 * ncol(unique1_1)))
unique1_1T <- cbind(unique1_1, unique1_1[,c(1, 2)], 0)
colnames(unique1_1T)[3:6] <- c("1", "mz 2", "rt 2", "2")
unique1_2T <- data.frame(matrix(nrow = 0, ncol = 2 * ncol(unique1_2)))
unique1_2T <- cbind(unique1_2[,c(1, 2)], 0, unique1_2)
colnames(unique1_2T)[3:6] <- c("1", "mz 2", "rt 2", "2")
AlignTable1 <- data.frame(matrix(nrow = 0, ncol = ncol(output1)))
AlignTable1 <- rbind(output1, unique1_1T, unique1_2T)
}
if (nrow(unique1_1) > 0 & nrow(unique1_2) == 0) {
  unique1_1T <- data.frame(matrix(nrow = 0, ncol = 2 * ncol(unique1_1)))
  unique1_1T <- cbind(unique1_1, unique1_1[,c(1, 2)], 0)
  colnames(unique1_1T)[3:6] <- c("1", "mz 2", "rt 2", "2")
  AlignTable1 <- data.frame(matrix(nrow = 0, ncol = ncol(output1)))
  AlignTable1 <- rbind(output1, unique1_1T)
}
if (nrow(unique1_1) == 0 & nrow(unique1_2) > 0) {
  unique1_2T <- data.frame(matrix(nrow = 0, ncol = 2 * ncol(unique1_2)))
  unique1_2T <- cbind(unique1_2[,c(1, 2)], 0, unique1_2)
  colnames(unique1_2T)[3:6] <- c("1", "mz 2", "rt 2", "2")
  AlignTable1 <- data.frame(matrix(nrow = 0, ncol = ncol(output)))
  AlignTable1 <- rbind(output, unique_2T)
}
###################################################################
# Feature alignment loop

for(k in 3:n){
  file <- read.csv(filename[k], stringsAsFactors = F)
  fileT <- data.frame(matrix(nrow = 0, ncol = 3))
  fileT <- file[, c("mz", "rt", "intensity")]
  colnames(fileT) <- c(paste("mz", k), paste("rt", k), k)
  output <- data.frame(matrix(nrow = 0, ncol = (ncol(AlignTable1) + 3)))
  for(j in 1:nrow(AlignTable1)){
    mass.lower.limit <- as.numeric(AlignTable1$mz[j]) - mass.tol
    mass.upper.limit <- as.numeric(AlignTable1$mz[j]) + mass.tol
    rt.lower.limit <- as.numeric(AlignTable1$rt[j]) - rt.tol
    rt.upper.limit <- as.numeric(AlignTable1$rt[j]) + rt.tol
    temp <- fileT[(as.numeric(fileT$mz) >= as.numeric(mass.lower.limit) &
                      as.numeric(fileT$mz) <= as.numeric(mass.upper.limit) &
                      as.numeric(fileT$rt) >= as.numeric(rt.lower.limit) &
                      as.numeric(fileT$rt) <= as.numeric(rt.upper.limit)),]
    temp <- temp[complete.cases(temp),]
    if(nrow(temp) > 0) {
      tmprow <- cbind(AlignTable1[j,], temp[1,])
      colnames(tmprow) <- colnames(output)
      output <- rbind(output, tmprow)
    }
  }
  colnames(output) <- c(colnames(AlignTable1), colnames(fileT))
  unique_1 <- anti_join(AlignTable1, as.data.frame(output[1:ncol(AlignTable1)]))
  unique_2 <- anti_join(fileT, as.data.frame(output[(ncol(AlignTable1)+1):ncol(output)]))
  if (nrow(unique_1) > 0 & nrow(unique_2) > 0) {
  unique_1T <- data.frame(matrix(nrow = 0, ncol = (ncol(unique_1) + 3)))
  unique_1T <- cbind(unique_1, unique_1[,c(1, 2)], 0)
  colnames(unique_1T) <- colnames(output)
  unique_2T <- data.frame(matrix(nrow = 0, ncol = (ncol(unique_2) * k)))
  m <- cbind(unique_2[,c(1, 2)], 0)
  unique_2T <- cbind(do.call(cbind, replicate(k-1, m, simplify=FALSE)), unique_2)
  colnames(unique_2T) <- colnames(output)
  AlignTable1 <- data.frame(matrix(nrow = 0, ncol = ncol(output)))
  AlignTable1 <- rbind(output, unique_1T, unique_2T)
  }
  if (nrow(unique_1) > 0 & nrow(unique_2) == 0) {
    unique_1T <- data.frame(matrix(nrow = 0, ncol = (ncol(unique_1) + 3)))
    unique_1T <- cbind(unique_1, unique_1[,c(1, 2)], 0)
    colnames(unique_1T) <- colnames(output)
    AlignTable1 <- data.frame(matrix(nrow = 0, ncol = ncol(output)))
    AlignTable1 <- rbind(output, unique_1T)
  }
  if (nrow(unique_1) == 0 & nrow(unique_2) > 0) {
    unique_2T <- data.frame(matrix(nrow = 0, ncol = (ncol(unique_2) * k)))
    m <- cbind(unique_2[,c(1, 2)], 0)
    unique_2T <- cbind(do.call(cbind, replicate(k-1, m, simplify=FALSE)), unique_2)
    colnames(unique_2T) <- colnames(output)
    AlignTable1 <- data.frame(matrix(nrow = 0, ncol = ncol(output)))
    AlignTable1 <- rbind(output, unique_2T)
  }
}
setwd(directory)
write.csv(AlignTable1, file="Aligned.csv", row.names = FALSE)

colSums(AlignTable1 != 0)
dereplicated <- data.frame(matrix(ncol = ncol(AlignTable1), nrow = 0))
colnames(dereplicated) <- colnames(AlignTable1)
for(n in 1:nrow(AlignTable1)) {
  mass.lower.limit <- AlignTable1$mz[n] - mass.tol
  mass.upper.limit <- AlignTable1$mz[n] + mass.tol
  rt.lower.limit <- AlignTable1$rt[n] - rt.tol
  rt.upper.limit <- AlignTable1$rt[n] + rt.tol
  temp <- dereplicated[dereplicated$mz >= mass.lower.limit & dereplicated$mz <= mass.upper.limit,]
  temp <- temp[temp$rt >= rt.lower.limit & temp$rt <= rt.upper.limit,]
  if(nrow(temp) == 0) {
    dereplicated[nrow(dereplicated) + 1,] = AlignTable1[n,]
  }
}
setwd(directory)
write.csv(dereplicated, file="dereplicatedAligned.csv", row.names = FALSE)



