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
FragIndex <- 30

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

#for(i in 1:length(MoleculeIDs)){
for(i in 1:5){
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
  } 
  if(length(MoleculeFiles) == 0){
    CountZero <- CountZero + 1
    MoleculesZero <- c(MoleculesZero, MoleculeID)    
  } 
}

Count
CountZero

########################################################################
## Draw histogram
########################################################################
## 1 pixel
Pixel1 <- subset(Data, range1 > 0)[,'range1']
MainTitle <- paste('Group', groupNum, 'Molecule', MoleculeID)
Xlabel <- '1 pixel dilation'

Pixel1_Norm <- Pixel1/median(Pixel1)

Hist1 <- qplot() + geom_histogram(aes(x = Pixel1)) + xlim(Xlim) + 
	ggtitle(label=MainTitle) + xlab(label=Xlabel) + 
  theme(plot.title = element_text(size = 10, colour = "gray20"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)
  )

Hist1_Norm <- qplot() + geom_histogram(aes(x = Pixel1_Norm)) + 
  ggtitle(label=paste(MainTitle, 'Normalized')) + xlab(label=Xlabel) + 
  theme(plot.title = element_text(size = 10, colour = "gray20"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)
  )

Density1 <- qplot() + geom_density(aes(x = Pixel1), kernel = 'epanechnikov', 
                       fill = 'gray20', col = 'gray20') + xlim(Xlim) + 
	ggtitle(label=MainTitle) + xlab(label=Xlabel) + 
  theme(plot.title = element_text(size = 10, colour = "gray20"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)
  )

Density1_norm <- qplot() + geom_density(aes(x = Pixel1_Norm), kernel = 'epanechnikov', 
                                        fill = 'gray20', col = 'gray20') + 
  ggtitle(label=MainTitle) + xlab(label=Xlabel) + 
  theme(plot.title = element_text(size = 10, colour = "gray20"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)
  )

DistFit <- fitdistr(x=Pixel1_Norm, densfun='gamma')

fitdistr(x=Pixel1, densfun='gamma')

xgamma <- rgamma(1000, shape = 5, rate = 0.1)
fitdistr(xgamma, "gamma")
qplot() + geom_density(aes(x = xgamma), kernel = 'epanechnikov')
