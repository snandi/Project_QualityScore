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
## Defining some constants and important variables
########################################################################
ConversionFactor <- 206
BasePairInterval <- ConversionFactor
BackbonePixels <- 1
DataPath.mm52_Intensities <- paste(DataPath.mm52, 'intensities_inca34_', BackbonePixels, 'pixel/', sep = '')
DataPath.mm52_Quality <- paste(DataPath.mm52, 'Project_QualityScore/', sep = '')

########################################################################
## Enter Fragment Details                                             
########################################################################
Chr <- 'chr13'
ChrNum <- gsub(pattern = 'chr', replacement = '', x = Chr)
FragIndex <- 7491

## Load GC Signal file
GC_Signal <- fn_load_GC_Signal(
  Chr = Chr, 
  DataPath.mm52 = DataPath.mm52_Intensities
)

GC_Signal <- subset(GC_Signal, Length_Pixels > 19)
FragIndices <- GC_Signal$FragIndex

bp.loc <- fn_load_bploc_mm52(
  Folder = '/ua/snandi/human_nMaps/GC_Content/', 
  Filename = paste0('mm52_all7431.goldOnly.bploc_', Chr), 
  ConversionFactor = ConversionFactor
)
########################################################################

IntensityData_inRange <- try(fn_saveTruncData_mm52(
  Chr = Chr, 
  FragIndex = FragIndex, 
  DataPath = DataPath.mm52_Intensities, 
  Truncate = FALSE, 
  TruncateLength = 0, 
  StretchPercentAllowed = 50, 
  Save = FALSE, 
  bp.loc = bp.loc
))

Filename <- paste0(DataPath.mm52, 'alignmentChunks/alignmentChunks.withLength.all7134Groups.goldOnly_', Chr)
Colnames <- c('refChr', 'refStartIndex', 'refEndIndex', 'molID', 'molStartIndex', 
              'molEndIndex', 'refStartCoord', 'refEndCoord', 'molStartCoord', 
              'molEndCoord', 'orient', 'lengthRatio')

Colnames.Steve <- c('refID', 'refStartIndex', 'refEndIndex', 'opID', 
                    'opStartIndex', 'opEndIndex', 'refStartCoord', 'refEndCoord', 'opStartCoord', 'opEndCoord', 
                    'orient', 'lengthRatio')

AlChunk <- fn_readAlignmentChunks(Filename = Filename, Header = FALSE, Colnames = Colnames, 
                                  Colnames.Steve = Colnames.Steve)

########################################################################
## Get MoleculeIDs for a fragIndex
########################################################################
print(FragIndex)
## Get only those molecules that have punctates both, at the beginning and end of the interval
AlChunk.Frag <- subset(AlChunk, refStartIndex ==   FragIndex & refEndIndex ==   (FragIndex + 1))
AlChunk.Frag$molID <- as.factor(AlChunk.Frag$molID)
#   str(AlChunk.Frag)
MoleculeID.Table <- table(AlChunk.Frag$molID)

## Discard moleculeIDs that have more than one fragment aligned to the same reference fragment
MoleculeID.Table[MoleculeID.Table > 1]
MoleculeIDs_MultipleFrag <- names(MoleculeID.Table[MoleculeID.Table > 1])
MoleculeID.Table <- MoleculeID.Table[MoleculeID.Table ==   1]
MoleculeIDs <- names(MoleculeID.Table)

MoleculeIDs_fromIntensity <- unique(as.vector(IntensityData_inRange$MoleculeID))

## Check if the same molecules are loaded, from Intensity profile and from Alchunk
sum(MoleculeIDs != MoleculeIDs_fromIntensity)  ## If this sum is non-zero, there is a problem

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
Filepath.Below75 <- paste0(DataPath.mm52_Quality, Chr, '/refFrag_', FragIndex, '/refFrag', FragIndex, '_Below75.txt')
system(command = paste('rm -f', Filepath.Below75))

i <- 8
ClusterMetrics <- c()
FilenamesBelow75 <- c()

#File_Straight <- fn_readStraightScoreFile(DataPath.mm52_Quality, FragIndex)

for(i in 1:length(MoleculeIDs)){
  #for(i in 1:50){
  MoleculeID <- MoleculeIDs[i]
  #StraightScore <- round(File_Straight$score[File_Straight$MoleculeID == MoleculeID], 4)
  StraightScore <- 1
  groupNum <- substr(MoleculeID, start = 1, stop = 7)
  MoleculeNum <- as.numeric(substr(MoleculeID, start = 13, stop = 19)) %% (ConversionFactor * 10000)
  
  Folderpath_Quality <- paste(DataPath.mm52_Quality, Chr, '/refFrag_', FragIndex, '/group1-', groupNum, 
                              '-inca34-outputs/', sep = '')
  Files <- try(list.files(path = Folderpath_Quality, pattern = paste('molecule', MoleculeNum, sep = '')))
  MoleculeFiles <- try(Files[grep(pattern = '.txt', x = Files)])
  PngFiles <- Files[grep(pattern = '.png', x = Files)]
  
  if(length(MoleculeFiles) ==   1){
    Filename.txt <- paste(Folderpath_Quality, MoleculeFiles, sep = '')
    Data <- read.table(Filename.txt, sep = ' ', header = T, stringsAsFactors = F)
    Data <- Data[,c('intensity1', 'intensity2', 'intensity3')]
  } else if(length(MoleculeFiles) > 1){
    print(groupNum)
    print(MoleculeFiles)
    Data <- c()
    for(MoleculeFile in MoleculeFiles){
      Filename.txt <- paste(Folderpath_Quality, MoleculeFile, sep = '')
      Data1 <- read.table(Filename.txt, sep = ' ', header = T, stringsAsFactors = F)
      Data1 <- Data1[,c('intensity1', 'intensity2', 'intensity3')]
      Data <- rbind(Data, Data1)
    }
  }
  pngImages <- vector("list", length(PngFiles))
  i_png <- 1
  
  if(length(PngFiles) > 0){
    for(PngFile in PngFiles){
      pngFile <- try(paste0(Folderpath_Quality, PngFile))
      pngImages[[i_png]] <- try(readPNG(source = pngFile, native = TRUE))
      #plot(pngImages[[i_png]])
      CountPNG <- CountPNG + 1
      i_png <- i_png + 1
    }
  }
  if(length(MoleculeFiles) > 0){
    Count <- Count + 1
    Xlim <- range(Data[ Data > 0 ])
    
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
    
    Filename.ClustPlot <- paste0('group', groupNum, '_molecule', MoleculeNum, '_GaussianCluster.pdf')
    Filepath.ClustPlot <- paste0(Folderpath_Quality, Filename.ClustPlot)
    
    pdf(file = Filepath.ClustPlot)
    
    ClusterOutput <- fn_ClusterPlotOutput(
      PixelData     = PixelData, 
      Molecule      = MoleculeID, 
      FragIndex     = FragIndex, 
      StraightScore = StraightScore
    )
    ClusterMetrics <- as.data.frame(rbind(ClusterMetrics, unlist(ClusterOutput[['ClusterMetrics']])), 
                                    stringsAsFactors = F)
    
    if(ClusterOutput[['ClusterMetrics']][['SurfaceNoiseScore']] <=  0.75){
      cat(paste0('group1-', groupNum, '-inca34-outputs/', Filename.ClustPlot), 
          file = Filepath.Below75, sep = '\n', append = T)  
    }
    print(ClusterOutput[['HistPlot']])
    
    library(grid)
    for(i_png in 1:length(pngImages)){
      pngImage <- pngImages[[i_png]]
      print(i_png)
      if(!(class(pngImage) ==  'try-error')){
        plot.new()
        try(grid.raster(pngImage, x = 0.5, y = 0.5, width = 1, height = 1))
      } else{
        print(paste(pngFile, 'Not valid'))
      }
    }
    dev.off()
  } 
  
  if(length(MoleculeFiles) ==   0){
    CountZero <- CountZero + 1
    MoleculesZero <- c(MoleculesZero, MoleculeID)    
  } 
}

ClusterMetrics.DF <- fn_formatClusterMetrics(ClusterMetrics)
fn_saveClusterMetrics_mm52(ClusterMetrics, DataPath.mm52_Quality, Chr, FragIndex)

IntensityData_inRange <- merge(
  x = IntensityData_inRange, 
  y = ClusterMetrics.DF[,c('MoleculeID', 'Discard')],
  by = 'MoleculeID'
)

Plot1 <- fn_plotCurves_diffLength_wOutliers(
  Data = IntensityData_inRange[,c('PixelNum', 'MoleculeID', 'Intensity', 'Discard')],
  Xlab = 'PixelNum',
  Ylab = 'Intensity',
  MainTitlePhrase = 'curves, with Image processed Outliers'
)

fn_loadFragmentFiles(
  Chr = Chr, 
  FragIndex = FragIndex,
  saveIn = .GlobalEnv,
  DataPath.mm52 = DataPath.mm52
)

########################################################################
## Drop nMaps that are beyond permissible stretch 
########################################################################
Pixels <- 1:nrow(Intensity_Eval)
StretchAllowed <- 0.15
Stretch <- aggregate(PixelNum ~ MoleculeID, data = IntensityData_inRange, FUN = max)
Stretch$PixelNum <- Stretch$PixelNum - 10 * (attributes(IntensityData_inRange)$TruncateLength == 0)
PixelMin <- (1 - StretchAllowed)*length(Pixels)
PixelMax <- (1 + StretchAllowed)*length(Pixels)
Stretch$StretchDrop <- (Stretch$PixelNum > PixelMax) | (Stretch$PixelNum < PixelMin)

Drop_forStretch <- as.vector(subset(Stretch, StretchDrop == TRUE)$MoleculeID)

if(sum(Stretch$StretchDrop > 0)) {
  Intensity_Eval <- Intensity_Eval[ , (colnames(Intensity_Eval) %w/o% Drop_forStretch) ]
  Intensity_Eval.D1 <- Intensity_Eval.D1[ , (colnames(Intensity_Eval.D1) %w/o% Drop_forStretch) ]
}

ClusterMetrics.DF <- merge(
  x = ClusterMetrics.DF,
  y = Stretch,
  by = 'MoleculeID'
  )

########################################################################
## Check for functional outliers
########################################################################
Outliers <- fn_getFunctionalOutliers(
  Curves = Intensity_Eval, 
  Xaxis = Pixels, 
  Names = list(main = paste(Chr, 'Frag', FragIndex), xlab = 'PixelPosition', ylab = 'Intensity'),
  N_Bootstrap = 500,
  trim = 0.10
)
Outliers$outliers

