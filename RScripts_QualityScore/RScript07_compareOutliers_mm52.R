rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script extends RScript03 and produces a table with more cluster
## outputs, to better classify "bad" molecules
########################################################################

########################################################################
## Load header files and source functions
########################################################################
RegistrationPath <- '~/R_Packages/Registration/R/'
Files <- list.files(path = RegistrationPath, pattern = ".R")
for(Script in Files) source(paste(RegistrationPath, Script, sep = ''))

PackagesLoaded <- loadPackages()
Packages <- PackagesLoaded$Packagesx
Packages_Par <- PackagesLoaded$Packages_Par

RScriptPath <- '~/Project_QualityScore/RScripts_QualityScore/'
source(paste(RScriptPath, 'fn_Library_CurveReg.R', sep = ''))
source(paste(RScriptPath, 'fn_Library_GC_Content.R', sep = ''))
source(paste(RScriptPath, 'fn_Library_Mflorum.R', sep = ''))
source(paste(RScriptPath, 'fn_Library_mm52.R', sep = ''))
source(paste(RScriptPath, 'fn_Library_Simulation.R', sep = ''))

source(paste(RScriptPath, 'fn_Library_QualityScore.R', sep = ''))

RPlotPath <- '~/Project_QualityScore/Plots/'
ProjectPath <- '~/Project_QualityScore/'
RDataPath <- '~/Project_QualityScore/RData/'
DataPath <- '~/Project_QualityScore/Data/'
DataPath.mf <- '/z/Proj/newtongroup/snandi/MF_cap348/'
DataPath.mm52 <- '/z/Proj/newtongroup/snandi/mm52-all7341/'

########################################################################
ConversionFactor <- 206
BasePairInterval <- ConversionFactor
BackbonePixels <- 1

Chr <- 'chr13'

DataPath.mm52_Intensities <- paste(DataPath.mm52, 'intensities_inca34_', BackbonePixels, 'pixel/', sep = '')
DataPath.mm52_Quality <- paste(DataPath.mm52, 'Project_QualityScore/', sep = '')

Filename <- paste0(DataPath.mm52_Quality, Chr, '/Chr13_Fragments.txt')
Frags <- read.table(file = Filename, header = FALSE, sep = ',')
Frags <- Frags$V1

AllOutliers <- c()
for(FragIndex in Frags){
  Folderpath_Quality <- paste( DataPath.mm52_Quality, Chr, '/refFrag_', FragIndex, '/', sep = '' )
  Filename <- paste0( Folderpath_Quality, 'OutlierComparisons.RData' )
  load(Filename)
  OutlierComparisons$FragIndex <- FragIndex
  AllOutliers <- rbind(AllOutliers, OutlierComparisons)
}

OutlierNames <- c(
  "Outliers_FM_Trim", 
  "Outliers_FM", 
  "Outliers_Mode_Trim", 
  "Outliers_Mode", 
  "Outliers_RTukey_Trim", 
  "Outliers_RTukey", 
  "Outliers_RProj_Trim", 
  "Outliers_RProj"
)
OutlierNames <- c(OutlierNames, 'Outliers_CumScore', 'Outliers_Any')

colnames(AllOutliers) <- gsub(pattern = 'Outliers_', replacement = '', x = colnames(AllOutliers))
OutlierNames <- gsub(pattern = 'Outliers_', replacement = '', x = OutlierNames)

#AllOutliers$CumScore <- (AllOutliers$CumScore > 0)
#View(AllOutliers[,c('FragIndex', 'Discard', 'Any', 'CumScore')])

AllOutliers$MoleculeID <- as.factor(AllOutliers$MoleculeID)
library(gdata)
AllOutliers$Curve <- as.numeric(AllOutliers$MoleculeID)

AllOutliers$Discard_Num <- (AllOutliers$Discard == 'Discard')
str(AllOutliers)

aggregate(Discard_Num ~ FragIndex, data = AllOutliers, FUN = sum)

subset(AllOutliers, FragIndex == 3326)[,c('Curve', 'Discard_Num', OutlierNames)]

## Type II error
n_D <- nrow(subset(AllOutliers, Discard == 'Discard'))
n_D/nrow(AllOutliers)
n_A <- colSums(subset(AllOutliers, Discard == 'Discard')[,c('Discard_Model', OutlierNames)])

Type2Errors <- c()
for(Frag in unique(AllOutliers$FragIndex)){
  Subset <- subset(AllOutliers, FragIndex == Frag & Discard == 'Discard')
  n_D <- nrow(Subset)
  if(n_D > 0)
  n_A <- colSums(Subset[,OutlierNames])/n_D
  Type2Errors <- rbind(Type2Errors, cbind(FragIndex = Frag, n_D = n_D, t(cbind(n_A))))
}

FragIndex <- 3326
