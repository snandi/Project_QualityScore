rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is copied from 
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
source(paste(RScriptPath, 'fn_Library_CurveReg.R', sep = ''))
source(paste(RScriptPath, 'fn_Library_GC_Content.R', sep = ''))
source(paste(RScriptPath, 'fn_Library_Mflorum.R', sep = ''))
source(paste(RScriptPath, 'fn_Library_QualityScore.R', sep = ''))

RPlotPath <- '~/Project_QualityScore/Plots/'
ProjectPath <- '~/Project_QualityScore/'
RDataPath <- '~/Project_QualityScore/RData/'
DataPath <- '~/Project_QualityScore/Data/'

DataPath.ngroup <- '/z/Proj/newtongroup/snandi/'
DataPath.ngroup.mf <- '/z/Proj/newtongroup/snandi/MF_cap348/'

########################################################################
## Outline of this script
########################################################################
## Step 1: Load Subinterval-1 pre-process the Nmaps and Register iteratively
## Step 2: Load Subinterval-2 pre-process the Nmaps and Register iteratively
## Step 3: Estimate T-Statistic of the registered Nmaps
## Step 4: Permute steps 1 to 3 a thousand times

########################################################################
## Defining some constants and important variables
########################################################################
ConversionFactor <- 209
BasePairInterval <- ConversionFactor
BackbonePixels <- 1
DataPath.mf <- paste(DataPath.ngroup.mf, 'intensities_inca34_', BackbonePixels, 'pixel/', sep='')

bp.loc <- fn_load_MF_bploc(ConversionFactor=ConversionFactor)
Today <- Sys.Date()

########################################################################
## Enter Fragment 1 Details
########################################################################
##First read in the arguments listed at the command line
Args <- (commandArgs(TRUE))
for(i in 1:length(Args)){
  eval(parse(text = Args[[i]]))
}

# PairwiseSimilarityThreshold <- 0.80 ##Command Line Argument

# NPermute <- 100             ## Command Line Argument
# NCores <- 20                ## Command Line Argument
# QScore_Qtl <- 0.70          ## Command Line Argument
# QScore_Cutoff <- NULL       ## Command Line Argument

Chr <- 'chr1'
ChrNum <- 1
# FragIndex_1 <- 36           ## Command Line Argument
# Regist_PixelFrom_1 <- 6     ## Command Line Argument
# Regist_PixelTo_1 <- 45      ## Command Line Argument

########################################################################
## Enter Fragment 2 Details
########################################################################
Chr <- 'chr1'
ChrNum <- 1
# FragIndex_2 <- 37           ## Command Line Argument
# Regist_PixelFrom_2 <- 6     ## Command Line Argument
# Regist_PixelTo_2 <- 45      ## Command Line Argument

NPixels <- Regist_PixelTo_1 - Regist_PixelFrom_1 + 1

DataPath.Out <- paste(DataPath.ngroup, 'MF_PermutationPlots/', 'Frag', FragIndex_1, 'P', 
                      Regist_PixelFrom_1, '-', Regist_PixelTo_1, '_AND_', 'Frag', 
                      FragIndex_2, 'P', Regist_PixelFrom_2, '-', Regist_PixelTo_2, '/', sep='')

dir.create(path = DataPath.Out, showWarnings = TRUE, recursive = FALSE, mode = "0755")

########################################################################
## Get Data to register
########################################################################
Data1_toRegist <- fn_prepDataForRegist_QualityScore(
  Chr                 = 'chr1', 
  FragIndex           = FragIndex_1, 
  DataPath.ngroup.mf  = DataPath.ngroup.mf,
  Regist_PixelFrom    = Regist_PixelFrom_1, 
  Regist_PixelTo      = Regist_PixelTo_1, 
  QScore_Cutoff       = QScore_Cutoff, 
  QScore_Qtl          = QScore_Qtl
)

# length(Data1_toRegist)
# names(Data1_toRegist)
# save(Data1_toRegist, file = 'Temp/Data_toRegist.RData')

## Extract all elements of the list
for(Name in names(Data1_toRegist)){
  Object <- Data1_toRegist[[Name]]
  assign(x=paste('Data1', Name, sep='_'), value=Object, envir=.GlobalEnv)
}
dim(Data1_Intensity_toRegist)

Data2_toRegist <- fn_prepDataForRegist_QualityScore(
  Chr                 = 'chr1', 
  FragIndex           = FragIndex_2, 
  DataPath.ngroup.mf  = DataPath.ngroup.mf,
  Regist_PixelFrom    = Regist_PixelFrom_2, 
  Regist_PixelTo      = Regist_PixelTo_2, 
  QScore_Cutoff       = QScore_Cutoff, 
  QScore_Qtl          = QScore_Qtl
)

## Extract all elements of the list
for(Name in names(Data2_toRegist)){
  Object <- Data2_toRegist[[Name]]
  assign(x=paste('Data2', Name, sep='_'), value=Object, envir=.GlobalEnv)
}
dim(Data2_Intensity_toRegist)

########################################################################
## Perform Iterated Registration of the actual Nmaps
########################################################################
Lambdas_ConstrainedWarping <- c(1, 0.3, 0.05)
IterRegData_1 <- fn_iteratedRegistration(Data_toRegist = Data1_toRegist, 
                                         Lambdas_ConstrainedWarping = Lambdas_ConstrainedWarping, 
                                         FragIndex = FragIndex_1, 
                                         Regist_PixelFrom = Regist_PixelFrom_1,
                                         Regist_PixelTo = Regist_PixelTo_1
)
Regfd_Initial_1 <- IterRegData_1[['Regfd_Initial']]
Regfd_Final_1 <- IterRegData_1[['Regfd_Final']]
RegisteredNmaps_1 <- IterRegData_1[['RegisteredCurves']]
RegisteredNmaps.D1_1 <- IterRegData_1[['RegisteredCurves.D1']]

IterRegData_2 <- fn_iteratedRegistration(Data_toRegist=Data2_toRegist, 
                                         Lambdas_ConstrainedWarping=Lambdas_ConstrainedWarping, 
                                         FragIndex=FragIndex_2, 
                                         Regist_PixelFrom=Regist_PixelFrom_2,
                                         Regist_PixelTo=Regist_PixelTo_2
)
Regfd_Initial_2 <- IterRegData_2[['Regfd_Initial']]
Regfd_Final_2 <- IterRegData_2[['Regfd_Final']]
RegisteredNmaps_2 <- IterRegData_2[['RegisteredCurves']]
RegisteredNmaps.D1_2 <- IterRegData_2[['RegisteredCurves.D1']]

TitleText_1 <- paste('MF Frag', FragIndex_1, '; Pixels', Regist_PixelFrom_1, 'To', 
                     Regist_PixelTo_1, ';', length(Data1_Molecules_forConsensus), 'Nmaps')
ForOrigPlot_1 <- plotAll_regist_fda(registOutput = Regfd_Initial_1, Lambda = 1, 
                                  TitleText = TitleText_1,
                                  Xlabel = paste('After Iteration', 1),
                                  Ylabel = 'Normalized, Centered & Scaled Intensity',
                                  Xarg_fine = Regist_PixelFrom_1:Regist_PixelTo_1,
                                  PlotBeforeRegist = FALSE,
                                  BeforeAfterDist = FALSE, saveToPDF = FALSE)
# ForOrigPlot_1[['AllPlots']]$Plot.Orig
Index <- length(Lambdas_ConstrainedWarping)
PlotsAndData_1 <- plotAll_regist_fda(registOutput = Regfd_Final_1, Lambda = 1, 
                                   TitleText = TitleText_1,
                                   Xlabel = paste('After Iteration', Index),
                                   Ylabel = 'Normalized, Centered & Scaled Intensity',
                                   Xarg_fine = Regist_PixelFrom_1:Regist_PixelTo_1,
                                   PlotBeforeRegist = FALSE,
                                   BeforeAfterDist = FALSE, saveToPDF = FALSE)
AllPlots_1 <- PlotsAndData_1[['AllPlots']]
# try(AllPlots_1[['Plot.Orig']])
# try(AllPlots_1[['Plot.Regist']])

TitleText_2 <- paste('MF Frag', FragIndex_2, '; Pixels', Regist_PixelFrom_2, 'To', 
                     Regist_PixelTo_2, ';', length(Data2_Molecules_forConsensus), 'Nmaps')
Index <- 1
ForOrigPlot_2 <- plotAll_regist_fda(registOutput = Regfd_Initial_2, Lambda = 1, 
                                    TitleText = TitleText_2,
                                    Xlabel = paste('After Iteration', 1),
                                    Ylabel = 'Normalized, Centered & Scaled Intensity',
                                    Xarg_fine = Regist_PixelFrom_2:Regist_PixelTo_2,
                                    PlotBeforeRegist = FALSE,
                                    BeforeAfterDist = FALSE, saveToPDF = FALSE)
# ForOrigPlot_2[['AllPlots']]$Plot.Orig

Index <- length(Lambdas_ConstrainedWarping)
PlotsAndData_2 <- plotAll_regist_fda(registOutput = Regfd_Final_2, Lambda = 1, 
                                     TitleText = TitleText_2,
                                     Xlabel = paste('After Iteration', 3),
                                     Ylabel = 'Normalized, Centered & Scaled Intensity',
                                     Xarg_fine = Regist_PixelFrom_2:Regist_PixelTo_2,
                                     PlotBeforeRegist = FALSE,
                                     BeforeAfterDist = FALSE, saveToPDF = FALSE)
AllPlots_2 <- PlotsAndData_2[['AllPlots']]
# try(AllPlots_2[['Plot.Orig']])
# try(AllPlots_2[['Plot.Regist']])

permute_T_Results_Before <- permute_pointwiseT(Mat1=Data1_Intensity_toRegist, 
                                        Mat2=Data2_Intensity_toRegist, 
                                        Nperm=10000, argvals=1:NPixels,   
                                        q=0.05, returnPlot=T)

# plot(permute_T_Results_Before$TPlot)
# plot(permute_T_Results_Before$Plot_pval)

# permute_T_Results_After <- permute_pointwiseT(Mat1=RegisteredNmaps_1, 
#                                         Mat2=RegisteredNmaps_2, 
#                                         Nperm=10000, argvals=1:NPixels,   
#                                         q=0.05, returnPlot=T)
# plot(permute_T_Results_After$TPlot)
# plot(permute_T_Results_After$Plot_pval)

TStats_Obs <- fn_absTStat(Mat1 = RegisteredNmaps_1, Mat2 = RegisteredNmaps_2)
########################################################################
## Perform Iterated Registration of the permuted Nmaps
########################################################################
NPermute <- NPermute
cl <- makeCluster(NCores)
registerDoParallel(cl)
PermutationOutput <- foreach(permutationIndex=1:NPermute, .inorder=FALSE, .combine= 'rbind', .packages=Packages_Par) %dopar% fn_permute_iterate(permutationIndex, Data1_toRegist, Data2_toRegist, DataPath.out=DataPath.Out, 
                     Lambdas_ConstrainedWarping = Lambdas_ConstrainedWarping)
stopCluster(cl)

PermutationOutput <- as.data.frame(PermutationOutput)
colnames(PermutationOutput) <- c('PermutationIndex', 'T_Sup', 'T_Median', 'T_Mean')
qval_TSup <- quantile(PermutationOutput$T_Sup, 0.95)
pval_TSup <- mean(TStats_Obs[['TSup']] < PermutationOutput$T_Sup)

Maintitle <- paste('Permutation test of iterated registration process \n', 
                   NPermute, 'permutations, after', length(Lambdas_ConstrainedWarping), 'iterations')

PermutePlot <- qplot(T_Sup, data = PermutationOutput) + geom_histogram() +
  geom_vline(xintercept = TStats_Obs[['TSup']], colour = 'royalblue1', size = 1) +
  geom_vline(xintercept = qval_TSup, colour = 'gray60', size = 1) +
  ggtitle(label=Maintitle) + xlab(label = paste('pvalue', pval_TSup)) +
  annotate(geom = 'text', x = TStats_Obs[['TSup']], y=Inf, vjust=-0.5, hjust=2, angle = 90,
           label='Obs', col='royalblue1') +
  annotate(geom='text', x = qval_TSup, y=Inf, vjust=-0.5, hjust=2, angle = 90, label='Crit', col='gray60')

########################################################################
## Save Output 
########################################################################
Filename.out <- paste(DataPath.Out, 'Permute_Iter_Regist_', 'Frag', FragIndex_1, 'P', 
                      Regist_PixelFrom_1, '-', Regist_PixelTo_1, '_AND_', 'Frag', 
                      FragIndex_2, 'P', Regist_PixelFrom_2, '-', Regist_PixelTo_2, '_', 
                      'QSqtl', (100*QScore_Qtl), '_', Today, '.pdf', sep='')
pdf(file=Filename.out, onefile=TRUE)
try(ForOrigPlot_1[['AllPlots']]$Plot.Orig)
try(AllPlots_1[['Plot.Regist']])

try(ForOrigPlot_2[['AllPlots']]$Plot.Orig)
try(AllPlots_2[['Plot.Regist']])

plot(permute_T_Results_Before$TPlot)
plot(permute_T_Results_Before$Plot_pval)

# plot(permute_T_Results_After$TPlot)
# plot(permute_T_Results_After$Plot_pval)

plot(PermutePlot)
dev.off()
######### Quit R without saving RData ###########
quit(save = 'no')
#################################################
