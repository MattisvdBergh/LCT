# source functions for LCT (available through https://github.com/MattisvdBergh/LCT)
source("C:/Users/Huub/Documents/Git/LCT/Growth.R")

# location of the Latent Gold 5.1 program
LG = "C:/Users/Huub/Latent_Gold_5.1/LatentGOLD5.1/lg51.exe"

# setwd to the analysis directory (line below works only in Rstudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Make the spss data into .dat file and add weights (only the first time)
library(foreign)
dataC = read.spss(file = "C:\\Users\\Huub\\surfdrive\\PhD\\LCTpackage\\data\\DFG.sav",
                 use.value.labels = FALSE, to.data.frame = TRUE)
dataTotalC = cbind(dataC, weight = 1)
dataTotalC[,"well"] = as.factor(dataTotalC[,"well"])

write.table(dataTotalC,
            file = "C:\\Users\\Huub\\surfdrive\\PhD\\LCTpackage\\data\\Crayen.dat",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

dataDir = "C:\\Users\\Huub\\surfdrive\\PhD\\LCTpackage\\data\\Crayen.dat"
dataDir = dataTotalC


dataDir = paste0(getwd(), "/Crayen.dat")
LGS = paste0(getwd(), "/SyntaxCrayen.LGS")

Results.C2 = LCGT(dataDir = dataDir,
                  LG = LG,
                  dependent = "well",
                  independent = c("time_cont", "time_cont2", "time_cont3"),
                  caseid = "UserID",
                  levelsDependent = 3,
                  resultsName = "testLCGT")

Results.C3 = LCGT(dataDir = dataDir,
                  LG = LG,
                  dependent = "well",
                  independent = c("time_cont", "time_cont2", "time_cont3"),
                  caseid = "UserID",
                  levelsDependent = 3,
                  resultsName = "testLCGT3",
                  maxClassSplit1 = 3)
