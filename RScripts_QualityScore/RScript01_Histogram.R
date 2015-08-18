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

########################################################################
## Read in data for a groupNum, frameNum & MoleculeID
########################################################################
groupNum <- '2433096'
frameNum <- 25
MoleculeID <-  99

Filename <- paste(ProjectPath, groupNum, '/', frameNum, '/', 'molecule', MoleculeID, '.txt', sep='')
Data <- read.table(Filename, sep=' ', header=T, stringsAsFactors=F)
Data <- Data[,1:3]
Xlim <- range(Data[Data>0])
########################################################################
## Draw histogram
########################################################################
## 1 pixel
Pixel1 <- subset(Data, range1 > 0)[,'range1']
MainTitle <- paste('Group', groupNum, 'Frame', frameNum, 'Molecule', MoleculeID)
Xlabel <- '1 pixel dilation'

Hist1 <- qplot() + geom_histogram(aes(x = Pixel1)) + xlim(Xlim) + 
	ggtitle(label=MainTitle) + xlab(label=Xlabel) + 
  theme(plot.title = element_text(size = 10, colour = "gray20"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)
  )

Density1 <- qplot() + geom_density(aes(x = Pixel1), kernel = 'epanechnikov', 
                       fill = 'turquoise2', col = 'turquoise2') + xlim(Xlim) + 
	ggtitle(label=MainTitle) + xlab(label=Xlabel) + 
  theme(plot.title = element_text(size = 10, colour = "gray20"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)
  )

## 2 pixel
Pixel2 <- subset(Data, range1 > 0)[,'range2']
MainTitle <- paste('Group', groupNum, 'Frame', frameNum, 'Molecule', MoleculeID)
Xlabel <- '2 pixel dilation'

Hist2 <- qplot() + geom_histogram(aes(x = Pixel2)) + xlim(Xlim) + 
	ggtitle(label=MainTitle) + xlab(label=Xlabel) + 
  theme(plot.title = element_text(size = 10, colour = "gray20"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)
  )

Density2 <- qplot() + geom_density(aes(x = Pixel2), kernel = 'epanechnikov', 
                       fill = 'turquoise2', col = 'turquoise2') + xlim(Xlim) + 
	ggtitle(label=MainTitle) + xlab(label=Xlabel) + 
  theme(plot.title = element_text(size = 10, colour = "gray20"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)
  )

## 3 pixel
Pixel3 <- subset(Data, range1 > 0)[,'range3']
MainTitle <- paste('Group', groupNum, 'Frame', frameNum, 'Molecule', MoleculeID)
Xlabel <- '3 pixel dilation'

Hist3 <- qplot() + geom_histogram(aes(x = Pixel3)) + xlim(Xlim) + 
	ggtitle(label=MainTitle) + xlab(label=Xlabel) + 
  theme(plot.title = element_text(size = 10, colour = "gray20"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)
  )

Density3 <- qplot() + geom_density(aes(x = Pixel3), kernel = 'epanechnikov', 
                       fill = 'turquoise2', col = 'turquoise2') + xlim(Xlim) + 
	ggtitle(label=MainTitle) + xlab(label=Xlabel) + 
  theme(plot.title = element_text(size = 10, colour = "gray20"), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)
  )

Filename.out <- paste(ProjectPath, groupNum, '/', frameNum, '/', 'molecule', MoleculeID, '.pdf', sep='')
pdf(file=Filename.out, pointsize=6)
grid.arrange(Hist1, Density1, Hist2, Density2, Hist3, Density3, ncol=2)
dev.off()

