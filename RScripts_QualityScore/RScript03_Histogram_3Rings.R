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
Packages <- PackagesLoaded$Packagesx
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
FragIndex <- 28

for(FragIndex in 28:37){
  print(FragIndex)
  ## Get only those molecules that have punctates both, at the beginning and end of the interval
  AlChunk.Frag <- subset(AlChunk, refStartIndex ==  FragIndex & refEndIndex ==  (FragIndex + 1))
  AlChunk.Frag$molID <- as.factor(AlChunk.Frag$molID)
  #   str(AlChunk.Frag)
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
      if(!(class(pngImage) == 'try-error')){
        plot.new()
        try(grid.raster(pngImage, x=0.5, y=0.5, width=1, height=1))
      } else{
        print(paste(pngFile, 'Not valid'))
      }
      dev.off()
    } 
    if(length(MoleculeFiles) ==  0){
      CountZero <- CountZero + 1
      MoleculesZero <- c(MoleculesZero, MoleculeID)    
    } 
  }
}
# Count                         ## No. of molecules with no text files in the folder
# CountZero                     ## No. of molecules which are completely in 1 frame
# CountPNG                      ## No. of molecules with png files
# CountMultipleFrames <- nrow(AlChunk.Frag) - Count
# CountMultipleFrames/Count     ## % of molecules which span across multiple frames

