rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

########################################################################
## This script reads in the pixel intensity values produced by Chengyue's
## python script and plots histograms. There will be three levels of pixel
## intensity values, one pixel, two pixels and three pixels around the 
## backbone of the molecules, after leaving out +/-2 pixels all around 
########################################################################

########################################################################
## Load header files and source functions
########################################################################
RegistrationPath <- '~/R_Packages/Registration/R/'
Files <- list.files(path = RegistrationPath, pattern = ".R")
for(Script in Files) source(paste(RegistrationPath, Script, sep = ''))

PackagesLoaded <- loadPackages()
Packages <- PackagesLoaded$Packages
Packages_Par <- PackagesLoaded$Packages_Par

RScriptPath <- '~/Project_QualityScore/RScripts_QualityScore/'
source(paste(RScriptPath, 'fn_Library_CurveReg.R', sep = ''))
source(paste(RScriptPath, 'fn_Library_GC_Content.R', sep = ''))
source(paste(RScriptPath, 'fn_Library_Mflorum.R', sep = ''))
source(paste(RScriptPath, 'fn_Library_QualityScore.R', sep = ''))

RPlotPath <- '~/Project_QualityScore/Plots/'
ProjectPath <- '~/Project_QualityScore/'
RDataPath <- '~/Project_QualityScore/RData/'
DataPath <- '~/Project_QualityScore/Data/'
DataPath.mf <- '/z/Proj/newtongroup/snandi/MF_cap348/'

########################################################################
## Defining some constants and important variables
########################################################################
ConversionFactor <- 209
BasePairInterval <- ConversionFactor
BackbonePixels <- 1
DataPath.mf_Intensities <- paste(DataPath.mf, 'intensities_inca34_', BackbonePixels, 'pixel/', sep = '')
DataPath.mf_Quality <- paste(DataPath.mf, 'Project_QualityScore/', sep = '')

bp.loc <- fn_load_MF_bploc(ConversionFactor = ConversionFactor)

Filename.Alchunk <- paste(DataPath.mf_Intensities, 'MF_cap348_inca34_cf209_minSize50_minFrag5_alignmentChunks.RData', sep = '')
load(Filename.Alchunk)

########################################################################
## Get MoleculeIDs for a fragIndex
########################################################################
FragIndex <- 9

## Get only those molecules that have punctates both, at the beginning and end of the interval
AlChunk.Frag <- subset(AlChunk, refStartIndex ==  FragIndex & refEndIndex ==  (FragIndex + 1))
AlChunk.Frag$molID <- as.factor(AlChunk.Frag$molID)
str(AlChunk.Frag)
MoleculeID.Table <- table(AlChunk.Frag$molID)

## Discard moleculeIDs that have more than one fragment aligned to the same reference fragment
MoleculeID.Table[MoleculeID.Table > 1]
MoleculeIDs_MultipleFrag <- names(MoleculeID.Table[MoleculeID.Table > 1])
MoleculeID.Table <- MoleculeID.Table[MoleculeID.Table ==  1]
MoleculeIDs <- names(MoleculeID.Table)

########################################################################
## Read in data for a groupNum, frameNum & MoleculeID
########################################################################
# groupNum <- '2433096'
# frameNum <- 25
# MoleculeID <-  99
Count <- 0
CountZero <- 0
CountPNG <- 0
MoleculesZero <- c()

i <- 1
AllPixelData <- c()

for(i in 1:length(MoleculeIDs)){
#for(i in 1:5){
  MoleculeID <- MoleculeIDs[i]
  
  groupNum <- substr(MoleculeID, start = 1, stop = 7)
  MoleculeNum <- as.numeric(substr(MoleculeID, start = 13, stop = 19)) %% (ConversionFactor * 10000)
  
  Folderpath_Quality <- paste(DataPath.mf_Quality, 'refFrag_', FragIndex, '/group1-', groupNum, 
                              '-inca34-outputs/', sep = '')
  Files <- try(list.files(path = Folderpath_Quality, pattern = paste('molecule', MoleculeNum, sep = '')))
  MoleculeFiles <- try(Files[grep(pattern = '.txt', x = Files)])
  PngFiles <- Files[grep(pattern = '.png', x = Files)]

  if(length(MoleculeFiles) ==  1){
    if(length(PngFiles) > 0){
      pngFile <- try(paste0(Folderpath_Quality, PngFiles))
      pngImage <- try(readPNG(source = pngFile, native = TRUE))
      #plot(pngImage)
      CountPNG <- CountPNG + 1
    }
    Count <- Count + 1
#     print(groupNum)
#     print(MoleculeFiles)
    Filename.txt <- paste(Folderpath_Quality, MoleculeFiles, sep = '')
    Data <- read.table(Filename.txt, sep = ' ', header = T, stringsAsFactors = F)
    Data <- Data[,1:3]
    Xlim <- range(Data[Data>0])
    
    Pixel1 <- subset(Data, intensity1 > 0)[,'intensity1']
    Pixel2 <- subset(Data, intensity2 > 0)[,'intensity2']
    Pixel3 <- subset(Data, intensity3 > 0)[,'intensity3']

    Pixel1_Norm <- Pixel1/median(c(Pixel1, Pixel2, Pixel3))
    Pixel2_Norm <- Pixel2/median(c(Pixel1, Pixel2, Pixel3))
    Pixel3_Norm <- Pixel3/median(c(Pixel1, Pixel2, Pixel3))
    
    PixelData <- as.data.frame(cbind(MoleculeID = MoleculeID, Pixel1 = Pixel1, Pixel1_Norm = Pixel1_Norm, 
                                     Pixel2 = Pixel2, Pixel2_Norm = Pixel2_Norm,
                                     Pixel3 = Pixel3, Pixel3_Norm = Pixel3_Norm))
    PixelData <- within(data = PixelData,{
      Pixel1 <- as.numeric(as.vector(Pixel1))
      Pixel1_Norm <- as.numeric(as.vector(Pixel1_Norm))
      Pixel2 <- as.numeric(as.vector(Pixel2))
      Pixel2_Norm <- as.numeric(as.vector(Pixel2_Norm))
      Pixel3 <- as.numeric(as.vector(Pixel3))
      Pixel3_Norm <- as.numeric(as.vector(Pixel3_Norm))
    })
    #AllPixelData <- rbind(AllPixelData, PixelData)
    Filename.ClustPlot <- paste(Folderpath_Quality, 'group', groupNum, '_molecule', MoleculeNum, 
                                '_GaussianCluster.pdf', sep='')
    pdf(file = Filename.ClustPlot)
    fn_ClusterAndPlot(PixelData = PixelData, Molecule = MoleculeID)
    dev.off()
  } 
  if(length(MoleculeFiles) ==  0){
    CountZero <- CountZero + 1
    MoleculesZero <- c(MoleculesZero, MoleculeID)    
  } 
}

Count                         ## No. of molecules with no text files in the folder
CountZero                     ## No. of molecules which are completely in 1 frame
CountPNG                      ## No. of molecules with png files
CountMultipleFrames <- nrow(AlChunk.Frag) - Count
CountMultipleFrames/Count     ## % of molecules which span across multiple frames

#PixelData_Split <- split(x = AllPixelData, f = AllPixelData$MoleculeID)
# Molecule = MoleculeIDs[25]
# 
# PixelData <- subset(AllPixelData, MoleculeID == Molecule)
# 
# PixelData_Long <- melt(data = PixelData[,c('MoleculeID', 'Pixel1', 'Pixel2', 'Pixel3')], 
#                        id.vars = 'MoleculeID', 
#                        measure.vars = c('Pixel1', 'Pixel2', 'Pixel3'))
# str(PixelData_Long)
# Maintitle <- paste('Reference Fragment', FragIndex, 'Molecule', Molecule)
# Hist1 <- ggplot(data = PixelData_Long, aes(x = log(value), col = variable, fill = variable)) + 
#   geom_histogram(binwidth = 0.01) + 
#   facet_wrap(~variable, ncol = 1) + 
#   geom_density(kernel = 'epanechnikov', col = 'gray60', lwd = 1) + 
#   ggtitle(label = Maintitle) + xlab(label = 'log of pixel intensities')
# 
# PixelData_Norm_Long <- melt(data = PixelData[,c('MoleculeID', 'Pixel1_Norm', 'Pixel2_Norm', 'Pixel3_Norm')], 
#                        id.vars = 'MoleculeID', 
#                        measure.vars = c('Pixel1_Norm', 'Pixel2_Norm', 'Pixel3_Norm'))
# 
# Model3 <- Mclust(data = log(PixelData[,'Pixel3_Norm']), G=1:4)
# summary(Model3)
# NClust3 <- Model3$G
# Model2 <- Mclust(data = log(PixelData[,'Pixel2_Norm']), G=1:4)
# summary(Model2)
# NClust2 <- Model2$G
# Model1 <- Mclust(data = log(PixelData[,'Pixel1_Norm']), G=1:4)
# summary(Model1)
# NClust1 <- Model1$G
# 
# PixelData$Cluster1 <- as.factor(Model1$classification)
# PixelData$Cluster2 <- as.factor(Model2$classification)
# PixelData$Cluster3 <- as.factor(Model3$classification)
# 
# ClusterMeans3 <- unlist(lapply(X=split(x = PixelData[,'Pixel3_Norm'], f = PixelData$Cluster3), FUN=mean))
# MaxDiff3 <- max(ClusterMeans3) - min(ClusterMeans3)
# ClusterMeans2 <- unlist(lapply(X=split(x = PixelData[,'Pixel2_Norm'], f = PixelData$Cluster2), FUN=mean))
# MaxDiff2 <- max(ClusterMeans2) - min(ClusterMeans2)
# ClusterMeans1 <- unlist(lapply(X=split(x = PixelData[,'Pixel1_Norm'], f = PixelData$Cluster1), FUN=mean))
# MaxDiff1 <- max(ClusterMeans1) - min(ClusterMeans1)
# 
# Discard <- ifelse(test = (MaxDiff3 > 0.5) && (MaxDiff2 > 0.3) , yes = 'Discard', no = 'Keep')
# MaxDiff3
# 
# SurfaceNoiseScore <- exp(-MaxDiff3)*exp(-MaxDiff2)*exp(-MaxDiff1)
# Maintitle <- paste('Reference Fragment', FragIndex, 'Molecule', Molecule, '\n', 
#                    'Surface Noise Score', round(SurfaceNoiseScore, 4), 
#                    ', Num of Clusters', NClust3, ',', Discard)
# 
# Hist1_Norm <- ggplot(data = PixelData_Norm_Long, aes(x = value, col = variable, fill = variable)) + 
#   geom_histogram(binwidth = 0.01) + 
#   facet_wrap(~variable, ncol = 1) +
#   geom_density(kernel = 'epanechnikov', col = 'gray50', lwd = 1) + 
#   ggtitle(label = Maintitle) + xlab(label = 'Normalized pixel intensities')
# 
# Hist1_Norm