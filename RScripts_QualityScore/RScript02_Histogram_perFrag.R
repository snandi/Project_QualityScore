rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
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
Files <- list.files(path=RegistrationPath, pattern=".R")
for(Script in Files) source(paste(RegistrationPath, Script, sep=''))

PackagesLoaded <- loadPackages()
Packages <- PackagesLoaded$Packages
Packages_Par <- PackagesLoaded$Packages_Par

RScriptPath <- '~/Project_QualityScore/RScripts_QualityScore/'
source(paste(RScriptPath, 'fn_Library_CurveReg.R', sep=''))
source(paste(RScriptPath, 'fn_Library_GC_Content.R', sep=''))
source(paste(RScriptPath, 'fn_Library_Mflorum.R', sep=''))
source(paste(RScriptPath, 'fn_Library_QualityScore.R', sep=''))

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
DataPath.mf_Intensities <- paste(DataPath.mf, 'intensities_inca34_', BackbonePixels, 'pixel/', sep='')
DataPath.mf_Quality <- paste(DataPath.mf, 'Project_QualityScore/', sep='')

bp.loc <- fn_load_MF_bploc(ConversionFactor=ConversionFactor)

Filename.Alchunk <- paste(DataPath.mf_Intensities, 'MF_cap348_inca34_cf209_minSize50_minFrag5_alignmentChunks.RData', sep='')
load(Filename.Alchunk)

########################################################################
## Get MoleculeIDs for a fragIndex
########################################################################
FragIndex <- 9

## Get only those molecules that have punctates both, at the beginning and end of the interval
AlChunk.Frag <- subset(AlChunk, refStartIndex == FragIndex & refEndIndex == (FragIndex + 1))
AlChunk.Frag$molID <- as.factor(AlChunk.Frag$molID)
str(AlChunk.Frag)
MoleculeID.Table <- table(AlChunk.Frag$molID)

## Discard moleculeIDs that have more than one fragment aligned to the same reference fragment
MoleculeID.Table[MoleculeID.Table > 1]
MoleculeIDs_MultipleFrag <- names(MoleculeID.Table[MoleculeID.Table > 1])
MoleculeID.Table <- MoleculeID.Table[MoleculeID.Table == 1]
MoleculeIDs <- names(MoleculeID.Table)

########################################################################
## Read in data for a groupNum, frameNum & MoleculeID
########################################################################
# groupNum <- '2433096'
# frameNum <- 25
# MoleculeID <-  99
Count <- 0
CountZero <- 0
MoleculesZero <- c()

AllPixelData <- c()
for(i in 1:length(MoleculeIDs)){
#for(i in 1:5){
    MoleculeID <- MoleculeIDs[i]
  
  groupNum <- substr(MoleculeID, start = 1, stop = 7)
  MoleculeNum <- as.numeric(substr(MoleculeID, start = 13, stop = 19)) %% (ConversionFactor * 10000)
  
  Folderpath_Quality <- paste(DataPath.mf_Quality, 'refFrag_', FragIndex, '/group1-', groupNum, 
                              '-inca34-outputs/', sep = '')
  MoleculeFiles <- try(list.files(path = Folderpath_Quality, pattern = paste('molecule', MoleculeNum, sep='')))
  if(length(MoleculeFiles) == 1){
    Count <- Count + 1
#     print(groupNum)
#     print(MoleculeFiles)
    Filename <- paste(Folderpath_Quality, MoleculeFiles, sep = '')
    Data <- read.table(Filename, sep=' ', header=T, stringsAsFactors=F)
    Data <- Data[,1:3]
    Xlim <- range(Data[Data>0])
    
    Pixel1 <- subset(Data, intensity1 > 0)[,'intensity1']
    Pixel1_Norm <- Pixel1/median(Pixel1)
    
    PixelData <- as.data.frame(cbind(Pixel1=Pixel1, Pixel1_Norm=Pixel1_Norm, MoleculeID=MoleculeID))
    PixelData <- within(data=PixelData,{
      Pixel1 <- as.numeric(as.vector(Pixel1))
      Pixel1_Norm <- as.numeric(as.vector(Pixel1_Norm))
    })
    AllPixelData <- rbind(AllPixelData, PixelData)
  } 
  if(length(MoleculeFiles) == 0){
    CountZero <- CountZero + 1
    MoleculesZero <- c(MoleculesZero, MoleculeID)    
  } 
}

Count
CountZero
CountMultipleFrames <- nrow(AlChunk.Frag) - Count
CountMultipleFrames/Count 

L <- split(x=AllPixelData, f=AllPixelData$MoleculeID)
fn_returnGammaPar <- function(DataToFit, Colname='Pixel1_Norm'){
  DataVectorToFit <- DataToFit[,Colname]
  MoleculeID <- as.vector(DataToFit$MoleculeID)[1]
  DistFit <- fitdistr(x=DataVectorToFit, densfun='gamma')
  Shape <- DistFit$estimate[['shape']]
  Rate <- DistFit$estimate[['rate']]

  KSTest <- ks.test(x=DataVectorToFit, y='pgamma', rate=Rate, shape=Shape)
  pValue <- KSTest$p.value
  Max <- max(DataVectorToFit)
  q95 <- quantile(x=DataVectorToFit, probs=0.95)
  #return(list(MoleculeID=MoleculeID, Shape=Shape, Rate=Rate))
  return(c(MoleculeID, Shape, Rate, pValue, Max, q95))
}

GammaParameters <- as.data.frame(do.call(what=rbind, lapply(X=L, FUN=fn_returnGammaPar, Colname='Pixel1_Norm')), 
                                 stringsAsFactors=FALSE)
colnames(GammaParameters) <- c('MoleculeID', 'Shape', 'Rate', 'pValue', 'Max', 'q95')
GammaParameters <- within(data=GammaParameters,{
  MoleculeID <- factor(MoleculeID)
  Shape <- as.numeric(Shape)
  Rate <- as.numeric(Rate)
  pValue <- round(as.numeric(pValue), 12)
  Max <- round(as.numeric(Max), 2)
  q95 <- round(as.numeric(q95), 2)
})
str(GammaParameters)

# ggplot() + geom_point(aes(x = Shape, y = Rate), data = GammaParameters)
# 
# ggplot() + geom_histogram(aes(x = Shape), data = GammaParameters)
# 
# ggplot() + geom_point(aes(x = Max, y = pValue), data = GammaParameters)
# 
# ggplot() + geom_histogram(aes(x = q95), data = GammaParameters) + 
#   geom_vline(xintercept=quantile(x=GammaParameters$q95, probs=0.95), col='royalblue1')
# 
# ggplot() + geom_point(aes(x = q95, y = pValue), data = GammaParameters, size=3) + 
#   geom_vline(xintercept=quantile(x=GammaParameters$q95, probs=0.95), col='royalblue1')

MaxPixel1_Norm <- max(subset(GammaParameters, pValue >= 0.1)[,'Max'])
MaxPixel1_Norm

Cuttoff <- min(quantile(x=GammaParameters$q95, probs=0.95), MaxPixel1_Norm)

nrow(subset(GammaParameters, Max<=Cuttoff))

PP <- L[[10]]$Pixel1_Norm
DistFit <- fitdistr(x=PP, densfun='gamma')
pdf.ke <- pdfCluster::kepdf(PP)

Molecule <- as.vector(L[[10]]$MoleculeID[1])
Shape <- subset(GammaParameters, MoleculeID==Molecule)[,'Shape']
Rate <- subset(GammaParameters, MoleculeID==Molecule)[,'Rate']

Discard <- ifelse(test=subset(GammaParameters, MoleculeID==Molecule)[,'Max'] < Cuttoff, yes='Keep', no='Discard')

Maintitle <- paste('Reference Fragment', FragIndex, 'Molecule', Molecule, Discard)
Hist1_Dens1 <- ggplot(data = subset(AllPixelData, MoleculeID==Molecule), aes(x = Pixel1_Norm)) + 
  geom_histogram(fill='gray60') + 
  geom_density(kernel = 'epanechnikov', col = 'gray20', lwd=1) + 
  stat_function(fun = dgamma, args=c(shape=Shape, rate=Rate), col='red', size=1) +
  geom_line(aes(x = pdf.ke@x, y = pdf.ke@estimate), col = 'royalblue1', size = 1) +
  ggtitle(label=Maintitle) + 


Hist1_Dens1
MoleculeIDs.Final <- unique(as.vector(AllPixelData$MoleculeID))
Molecule <- MoleculeIDs.Final[23]

fn_plotDensities <- function(AllPixelData, Molecule, GammaParameters, Cuttoff){
 PixelData <- as.data.frame(subset(AllPixelData, MoleculeID == Molecule))
 PP <- PixelData$Pixel1_Norm
 DistFit <- fitdistr(x=PP, densfun='gamma')
 pdf.ke <- pdfCluster::kepdf(PP)
 Shape <- subset(GammaParameters, MoleculeID == Molecule)[,'Shape']
 Rate <- subset(GammaParameters, MoleculeID == Molecule)[,'Rate']
 Discard <- ifelse(test=subset(GammaParameters, MoleculeID==Molecule)[,'Max'] < Cuttoff, yes='Keep', no='Discard')
 Maintitle <- paste('Reference Fragment', FragIndex, 'Molecule', Molecule, Discard)
 Hist_Dens <- ggplot(data = PixelData, aes(x = Pixel1_Norm)) + 
   geom_histogram(fill='gray60') + 
   geom_density(kernel = 'epanechnikov', col = 'gray20', lwd=1) + 
   stat_function(fun = dgamma, args=c(shape=Shape, rate=Rate), col='red', size=1) +
   geom_line(aes(x = pdf.ke@x, y = pdf.ke@estimate), col = 'royalblue1', size = 1) +
   ggtitle(label=Maintitle) + ylab(label='') + xlab(label='Normalized Pixel intensities')
 return(Hist_Dens)
}

Filename.plot <- paste(RPlotPath, 'refFrag', FragIndex, '_DensityPlots.pdf', sep='')
pdf(file=Filename.plot, onefile=TRUE)
for(Molecule in MoleculeIDs.Final){
#   Plot <- fn_plotDensities(AllPixelData=AllPixelData, Molecule=Molecule, 
#                            GammaParameters=GammaParameters, Cuttoff=Cuttoff)
#   try(print(Plot))
  PixelData <- as.data.frame(subset(AllPixelData, MoleculeID == Molecule))
  PP <- PixelData$Pixel1_Norm
  DistFit <- fitdistr(x=PP, densfun='gamma')
  pdf.ke <- pdfCluster::kepdf(PP)
  Shape <- subset(GammaParameters, MoleculeID == Molecule)[,'Shape']
  Rate <- subset(GammaParameters, MoleculeID == Molecule)[,'Rate']
  Discard <- ifelse(test=subset(GammaParameters, MoleculeID==Molecule)[,'Max'] < Cuttoff, yes='Keep', no='Discard')
  Maintitle <- paste('Reference Fragment', FragIndex, 'Molecule', Molecule, Discard)
  Hist_Dens <- ggplot(data = PixelData, aes(x = Pixel1_Norm)) + 
    geom_histogram(fill='gray60') + 
    geom_density(kernel = 'epanechnikov', col = 'gray20', lwd=1) + 
    stat_function(fun = dgamma, args=c(shape=Shape, rate=Rate), col='red', size=1) +
    geom_line(aes(x = pdf.ke@x, y = pdf.ke@estimate), col = 'royalblue1', size = 1) +
    ggtitle(label=Maintitle) + ylab(label='') + xlab(label='Normalized Pixel intensities')
  try(print(Hist_Dens))
}
dev.off()
