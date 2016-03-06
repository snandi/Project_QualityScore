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

# Args <- (commandArgs(TRUE))
# for(i in 1:length(Args)){
#   eval(parse(text = Args[[i]]))
# }
########################################################################

fn_detectOutliers_mm52 <- function(
  Chr, 
  FragIndex, 
  DataPath.mm52
){
  
  ########################################################################
  ## Defining some constants and important variables
  ########################################################################
  ConversionFactor <- 206
  BasePairInterval <- ConversionFactor
  BackbonePixels <- 1
  
  ########################################################################
  ## Enter Fragment Details                                             
  ########################################################################
  ChrNum <- gsub(pattern = 'chr', replacement = '', x = Chr)
  #FragIndex <- 7465
  
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
  
  ClusterMetrics <- fn_loadClusterMetrics_mm52(ClusterMetrics, DataPath.mm52_Quality, Chr, FragIndex)
  
  IntensityData_inRange <- merge(
    x = IntensityData_inRange, 
    y = ClusterMetrics[,c('MoleculeID', 'Discard', 'Discard_Model')],
    by = 'MoleculeID'
  )
  
  Plot_Discard <- fn_plotCurves_diffLength_wOutliers(
    Data = IntensityData_inRange[,c('PixelNum', 'MoleculeID', 'Intensity', 'Discard')],
    Xlab = 'PixelNum',
    Ylab = 'Intensity',
    MainTitlePhrase = 'Original curves, with Image processed Outliers'
  )
  Plot_Discard <- Plot_Discard + scale_color_manual(values=c("red", "gray60"))
  
  Plot_Discard_Model <- fn_plotCurves_diffLength_wOutliers(
    Data = IntensityData_inRange[,c('PixelNum', 'MoleculeID', 'Intensity', 'Discard_Model')],
    Xlab = 'PixelNum',
    Ylab = 'Intensity',
    MainTitlePhrase = 'Original curves, with Image processed Outliers, Model based'
  )
  Plot_Discard_Model <- Plot_Discard_Model + scale_color_manual(values=c("gray60", "red"))
  
  fn_loadFragmentFiles(
    Chr = Chr, 
    FragIndex = FragIndex,
    #saveIn = .GlobalEnv,
    saveIn = as.environment(-1),
    DataPath.mm52 = DataPath.mm52
  )
  
  ########################################################################
  ## Drop nMaps that are beyond permissible stretch 
  ########################################################################
  Pixels <- 1:nrow(Intensity_Eval)
  StretchAllowed <- 0.20
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
  
  ClusterMetrics <- merge(
    x = ClusterMetrics,
    y = Stretch,
    by = 'MoleculeID'
  )
  
  ########################################################################
  ## Check for functional outliers using 
  ########################################################################
  ## FM depth, Trim Yes
  Outliers_FM_Trim <- fn_getOutliers (
    Curves = Intensity_Eval, 
    Xaxis = Pixels, 
    Names = list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity'),
    DepthType = 'FM',
    N_Bootstrap = 500,
    Trim = 'Yes',
    TrimPct = 0.15
  )
  Outliers_FM_Trim$outliers
  ClusterMetrics$Outliers_FM_Trim <- ClusterMetrics$MoleculeID %in% Outliers_FM_Trim$outliers
  
  ## FM depth, Trim No
  Outliers_FM <- fn_getOutliers (
    Curves = Intensity_Eval, 
    Xaxis = Pixels, 
    Names = list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity'),
    DepthType = 'FM',
    N_Bootstrap = 500,
    Trim = 'No',
    TrimPct = 0
  )
  Outliers_FM$outliers
  ClusterMetrics$Outliers_FM <- ClusterMetrics$MoleculeID %in% Outliers_FM$outliers
  
  ## Modal depth, Trim Yes
  Outliers_Mode_Trim <- fn_getOutliers (
    Curves = Intensity_Eval, 
    Xaxis = Pixels, 
    Names = list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity'),
    DepthType = 'Mode',
    N_Bootstrap = 500,
    Trim = 'Yes',
    TrimPct = 0.15
  )
  Outliers_Mode_Trim$outliers
  ClusterMetrics$Outliers_Mode_Trim <- ClusterMetrics$MoleculeID %in% Outliers_Mode_Trim$outliers
  
  ## Modal depth, Trim No
  Outliers_Mode <- fn_getOutliers (
    Curves = Intensity_Eval, 
    Xaxis = Pixels, 
    Names = list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity'),
    DepthType = 'Mode',
    N_Bootstrap = 500,
    Trim = 'No',
    TrimPct = 0
  )
  Outliers_Mode$outliers
  ClusterMetrics$Outliers_Mode <- ClusterMetrics$MoleculeID %in% Outliers_Mode$outliers
  
  ## RTukey depth, Trim Yes
  Outliers_RTukey_Trim <- fn_getOutliers (
    Curves = Intensity_Eval, 
    Xaxis = Pixels, 
    Names = list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity'),
    DepthType = 'RTukey',
    N_Bootstrap = 500,
    Trim = 'Yes',
    TrimPct = 0.15
  )
  Outliers_RTukey_Trim$outliers
  ClusterMetrics$Outliers_RTukey_Trim <- ClusterMetrics$MoleculeID %in% Outliers_RTukey_Trim$outliers
  
  ## RTukey depth, Trim No
  Outliers_RTukey <- fn_getOutliers (
    Curves = Intensity_Eval, 
    Xaxis = Pixels, 
    Names = list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity'),
    DepthType = 'RTukey',
    N_Bootstrap = 500,
    Trim = 'No',
    TrimPct = 0
  )
  Outliers_RTukey$outliers
  ClusterMetrics$Outliers_RTukey <- ClusterMetrics$MoleculeID %in% Outliers_RTukey$outliers
  
  ## RProj depth, Trim Yes
  Outliers_RProj_Trim <- fn_getOutliers (
    Curves = Intensity_Eval, 
    Xaxis = Pixels, 
    Names = list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity'),
    DepthType = 'RProj',
    N_Bootstrap = 500,
    Trim = 'Yes',
    TrimPct = 0.15
  )
  Outliers_RProj_Trim$outliers
  ClusterMetrics$Outliers_RProj_Trim <- ClusterMetrics$MoleculeID %in% Outliers_RProj_Trim$outliers
  
  ## RProj depth, Trim No
  Outliers_RProj <- fn_getOutliers (
    Curves = Intensity_Eval, 
    Xaxis = Pixels, 
    Names = list(main = 'Main', xlab = 'PixelPosition', ylab = 'Intensity'),
    DepthType = 'RProj',
    N_Bootstrap = 500,
    Trim = 'No',
    TrimPct = 0
  )
  Outliers_RProj$outliers
  ClusterMetrics$Outliers_RProj <- ClusterMetrics$MoleculeID %in% Outliers_RProj$outliers
  
  ########################################################################
  ## Cumulative Outlier score
  ########################################################################
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
  ClusterMetrics$Outliers_CumScore <- rowSums(x = ClusterMetrics[,OutlierNames])
  
  ClusterMetrics$Outliers_Any <- ( ClusterMetrics$Outliers_CumScore > 0 )
  
  OutlierNames <- c(OutlierNames, 'Outliers_CumScore', 'Outliers_Any')
  #View(ClusterMetrics[,c('MoleculeID', 'Discard', 'Discard_Model', 'Outliers_CumScore', 'Outliers_Any')])
  
  #PlotSmooth
  
  PlotSmooth_Discard <- fn_plotMultCurves_inGroups(
    Data = Intensity_Eval, 
    Groups = ClusterMetrics[,c('MoleculeID', 'Discard')],
    GroupColors = c("red", "gray50"),
    XVar = Pixels, 
    MainTitlePhrase = 'Smoothed curves, with Image processed Outliers')
  PlotSmooth_Discard <- PlotSmooth_Discard + theme_gray() +
    theme(legend.position = '')
  
  PlotSmooth_Outliers_Any <- fn_plotMultCurves_inGroups(
    Data = Intensity_Eval, 
    Groups = ClusterMetrics[,c('MoleculeID', 'Outliers_Any')],
    GroupColors = c("gray50", "red"),
    XVar = Pixels, 
    MainTitlePhrase = 'Smoothed curves, with Depth-based Outliers')
  PlotSmooth_Outliers_Any <- PlotSmooth_Outliers_Any + theme_gray() +
    theme(legend.position = '')
  
  # grid.arrange(PlotSmooth_Discard, PlotSmooth_Outliers_Any, ncol = 1)
  
  OutlierComparisons <- ClusterMetrics[,c('MoleculeID', 'Discard', 'Discard_Model', OutlierNames)]
  
  fn_saveOutlierComparison_mm52( OutlierComparisons, DataPath.mm52_Quality, Chr, FragIndex )
  
  Folderpath_Quality <- paste( DataPath.mm52_Quality, Chr, '/refFrag_', FragIndex, '/', sep = '' )
  Filename <- paste0( Folderpath_Quality, 'OutlierComparisons_', FragIndex, '.pdf' )
  pdf(file = Filename, onefile = TRUE)
  try(plot(Plot_Discard))
  try(plot(Plot_Discard_Model))
  try(plot(PlotSmooth_Discard))
  try(plot(PlotSmooth_Outliers_Any))
  try(dev.off())
  try(dev.off())
}

# fn_detectOutliers_mm52(
#   Chr = Chr, 
#   FragIndex = FragIndex, 
#   DataPath.mm52 = DataPath.mm52
# )



########################################################################
## For parallel execution
########################################################################
Filename <- paste0(DataPath.mm52_Quality, Chr, '/Chr13_Fragments.txt')
Frags <- read.table(file = Filename, header = FALSE, sep = ',')
Frags <- Frags$V1
NCores <- 22
#cl <- makeCluster(NCores)
cl <- makePSOCKcluster(NCores)
doParallel::registerDoParallel(cl)
foreach(FragIndex = Frags, .inorder=FALSE, .packages=Packages_Par) %dopar% fn_detectOutliers_mm52(Chr = Chr, FragIndex, DataPath.mm52)
stopCluster(cl)
gc()
