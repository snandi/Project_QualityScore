## This function library will contain functions exclusively used for the Quality Score project

##########################################################################################
## This function clusters the intensities around the molecules and plots the histograms
## This is called by RScript03
##########################################################################################
fn_ClusterAndPlot <- function(PixelData, Molecule, FragIndex, ...){
  PixelData_Norm_Long <- melt(data = PixelData[,c('MoleculeID', 'Pixel1_Norm', 'Pixel2_Norm', 'Pixel3_Norm')], 
                              id.vars = 'MoleculeID', 
                              measure.vars = c('Pixel1_Norm', 'Pixel2_Norm', 'Pixel3_Norm'))
  Model3 <- Mclust(data = log(PixelData[,'Pixel3_Norm']), G=1:4)
  NClust3 <- Model3$G
  clValid::connectivity(clusters = Model3$classification, Data = log(PixelData[,'Pixel3_Norm']))
  Model2 <- Mclust(data = log(PixelData[,'Pixel2_Norm']), G=1:4)
  NClust2 <- Model2$G
  Model1 <- Mclust(data = log(PixelData[,'Pixel1_Norm']), G=1:4)
  NClust1 <- Model1$G
  clValid::connectivity(clusters = Model1$classification, Data = log(PixelData[,'Pixel1_Norm']))
  
  PixelData$Cluster1 <- as.factor(Model1$classification)
  PixelData$Cluster2 <- as.factor(Model2$classification)
  PixelData$Cluster3 <- as.factor(Model3$classification)
  
  ClusterMeans3 <- unlist(lapply(X=split(x = PixelData[,'Pixel3_Norm'], f = PixelData$Cluster3), FUN=mean))
  MaxDiff3 <- max(ClusterMeans3) - min(ClusterMeans3)
  ClusterMeans2 <- unlist(lapply(X=split(x = PixelData[,'Pixel2_Norm'], f = PixelData$Cluster2), FUN=mean))
  MaxDiff2 <- max(ClusterMeans2) - min(ClusterMeans2)
  ClusterMeans1 <- unlist(lapply(X=split(x = PixelData[,'Pixel1_Norm'], f = PixelData$Cluster1), FUN=mean))
  MaxDiff1 <- max(ClusterMeans1) - min(ClusterMeans1)
  
  Discard <- ifelse(test = (MaxDiff3 > 0.5) && (MaxDiff2 > 0.3) , yes = 'Discard', no = 'Keep')
  
  SurfaceNoiseScore <- exp(-MaxDiff3)*exp(-MaxDiff2)*exp(-MaxDiff1)
  if(SurfaceNoiseScore < 0.6) Discard <- 'Discard'
  
  Maintitle <- paste('Reference Fragment', FragIndex, 'Molecule', Molecule, '\n', 
                     'Surface Noise Score', round(SurfaceNoiseScore, 4), 
                     ', Num of Clusters', NClust3, ',', Discard)
  
  Hist1_Norm <- ggplot(data = PixelData_Norm_Long, aes(x = value, col = variable, fill = variable)) + 
    geom_histogram(binwidth = 0.01) + 
    facet_wrap(~variable, ncol = 1) +
    geom_density(kernel = 'epanechnikov', col = 'gray50', lwd = 1) + 
    ggtitle(label = Maintitle) + xlab(label = 'Normalized pixel intensities')
  
  print(Hist1_Norm)
}
##########################################################################################

##########################################################################################
## This function outputs cluster metrics, in addition to the histogram plots. This is 
## the only change to fn_ClusterAndPlot. This is called by RScript04
##########################################################################################
fn_ClusterPlotOutput <- function(PixelData, Molecule, FragIndex, StraightScore, ...){
  PixelData_Norm_Long <- melt(data = PixelData[,c('MoleculeID', 'Pixel1_Norm', 'Pixel2_Norm', 'Pixel3_Norm')], 
                              id.vars = 'MoleculeID', 
                              measure.vars = c('Pixel1_Norm', 'Pixel2_Norm', 'Pixel3_Norm'))
  
  Model3 <- Mclust(data = log(PixelData[,'Pixel3_Norm']), G=1:4)
  NClust3 <- Model3$G
  Conn3 <- clValid::connectivity(clusters = Model3$classification, Data = log(PixelData[,'Pixel3_Norm']))
  Dunn3 <- clValid::dunn(clusters = Model3$classification, Data = log(PixelData[,'Pixel3_Norm']))
  
  Model2 <- Mclust(data = log(PixelData[,'Pixel2_Norm']), G=1:4)
  NClust2 <- Model2$G
  Conn2 <- clValid::connectivity(clusters = Model2$classification, Data = log(PixelData[,'Pixel2_Norm']))
  Dunn2 <- clValid::dunn(clusters = Model2$classification, Data = log(PixelData[,'Pixel2_Norm']))
  
  Model1 <- Mclust(data = log(PixelData[,'Pixel1_Norm']), G=1:4)
  NClust1 <- Model1$G
  Conn1 <- clValid::connectivity(clusters = Model1$classification, Data = log(PixelData[,'Pixel1_Norm']))
  Dunn1 <- clValid::dunn(clusters = Model1$classification, Data = log(PixelData[,'Pixel1_Norm']))
  
  PixelData$Cluster1 <- as.factor(Model1$classification)
  PixelData$Cluster2 <- as.factor(Model2$classification)
  PixelData$Cluster3 <- as.factor(Model3$classification)
  
  ClusterMeans3 <- unlist(lapply(X=split(x = PixelData[,'Pixel3_Norm'], f = PixelData$Cluster3), FUN=mean))
  MaxDiff3 <- max(ClusterMeans3) - min(ClusterMeans3)
  ClusterMeans2 <- unlist(lapply(X=split(x = PixelData[,'Pixel2_Norm'], f = PixelData$Cluster2), FUN=mean))
  MaxDiff2 <- max(ClusterMeans2) - min(ClusterMeans2)
  ClusterMeans1 <- unlist(lapply(X=split(x = PixelData[,'Pixel1_Norm'], f = PixelData$Cluster1), FUN=mean))
  MaxDiff1 <- max(ClusterMeans1) - min(ClusterMeans1)
  
  Discard <- ifelse(test = (MaxDiff3 > 0.5) && (MaxDiff2 > 0.3) , yes = 'Discard', no = 'Keep')
  
  SurfaceNoiseScore <- exp(-MaxDiff3)*exp(-MaxDiff2)*exp(-MaxDiff1)
  ## SurfaceNoiseScore <- exp(-MaxDiff3*(3/9))*exp(-MaxDiff2*(3/9))*exp(-MaxDiff1*(3/9))
  ## Different layers can be given different weights
  
  if(SurfaceNoiseScore < 0.65) Discard <- 'Discard'
  
  ClusterMetrics <- list(MoleculeID = Molecule, 
                         NClust1 = NClust1, MaxDiff1 = MaxDiff1, Conn1 = Conn1, Dunn1 = Dunn1, 
                         NClust2 = NClust2, MaxDiff2 = MaxDiff2, Conn2 = Conn2, Dunn2 = Dunn2, 
                         NClust3 = NClust3, MaxDiff3 = MaxDiff3, Conn3 = Conn3, Dunn3 = Dunn3, 
                         SurfaceNoiseScore = SurfaceNoiseScore, StraightScore = StraightScore,
                         Discard = Discard
  )
  
  Maintitle <- paste('Ref Frag', FragIndex, 'Molecule', Molecule, '\n', 
                     'Surface Noise Score:', round(SurfaceNoiseScore, 4), 
                     'Straight Score:', StraightScore, '\n', 
                     ', # of Clusters in 3rd layer,', NClust3, 'Decision:', Discard)
  
  Hist1_Norm <- ggplot(data = PixelData_Norm_Long, aes(x = value, col = variable, fill = variable)) + 
    geom_histogram(binwidth = 0.01) + 
    facet_wrap(~variable, ncol = 1) +
    geom_density(kernel = 'epanechnikov', col = 'gray50', lwd = 1) + 
    ggtitle(label = Maintitle) + xlab(label = 'Normalized pixel intensities') +
    theme(plot.title = element_text(hjust = 0))
  
  # print(Hist1_Norm)
  return(list(ClusterMetrics = ClusterMetrics, HistPlot = Hist1_Norm))
}

# The Dunn Index is the ratio of the smallest distance between observations not in the same cluster to
# the largest intra-cluster distance. The Dunn Index has a value between zero and infinity, and should
# be maximized. For details see the package vignette.


# The connectivity indicates the degree of connectedness of the clusters, as determined by the knearest
# neighbors. The neighbSize argument specifies the number of neighbors to use. The connectivity
# has a value between 0 and infinity and should be minimized. For details see the package
# vignette.

##########################################################################################
## This function formats & saves the ClusterMetrics data
##########################################################################################
fn_saveClusterMetrics <- function(ClusterMetrics, DataPath.mf_Quality, FragIndex){
  ClusterMetrics <- as.data.frame(ClusterMetrics, StringsAsFactors = F)
  
  ClusterMetrics <- within(data = ClusterMetrics,{
    NClust1 <- round(as.numeric(as.vector(NClust1)), 4)
    MaxDiff1 <- round(as.numeric(as.vector(MaxDiff1)), 4)
    Conn1 <- round(as.numeric(as.vector(Conn1)), 4)
    Dunn1 <- round(as.numeric(as.vector(Dunn1)), 4)
    
    NClust2 <- round(as.numeric(as.vector(NClust2)), 4)
    MaxDiff2 <- round(as.numeric(as.vector(MaxDiff2)), 4)
    Conn2 <- round(as.numeric(as.vector(Conn2)), 4)
    Dunn2 <- round(as.numeric(as.vector(Dunn2)), 4)
    
    NClust3 <- round(as.numeric(as.vector(NClust3)), 4)
    MaxDiff3 <- round(as.numeric(as.vector(MaxDiff3)), 4)
    Conn3 <- round(as.numeric(as.vector(Conn3)), 4)
    Dunn3 <- round(as.numeric(as.vector(Dunn3)), 4)
    
    SurfaceNoiseScore <- round(as.numeric(as.vector(SurfaceNoiseScore)), 4)
    StraightScore <- round(as.numeric(as.vector(StraightScore)), 4)
  })
  
  Below75 <- subset(ClusterMetrics, SurfaceNoiseScore < 0.75)
  
  Filename.ClustMetrics <- paste0(DataPath.mf_Quality, 'refFrag_', FragIndex, '/refFrag', FragIndex, 
                                  '_SurfaceNoiseScore.txt')
  write.table(x=ClusterMetrics, file=Filename.ClustMetrics, quote=FALSE, sep='\t', row.names=F)
  
  Filename.ClustMetrics <- paste0(DataPath.mf_Quality, 'refFrag_', FragIndex, '/refFrag', FragIndex, 
                                  '_SurfaceNoiseScoreBelow75.txt')
  write.table(x=Below75, file=Filename.ClustMetrics, quote=FALSE, sep='\t', row.names=F)  
  
  N_Discard <- nrow(subset(ClusterMetrics, Discard == 'Discard'))
  N_Total <- nrow(ClusterMetrics)
  
  Maintitle <- paste(N_Discard, 'Discarded, out of', N_Total, 'Molecules, in Ref Frag', FragIndex)
  
  ScoreHist <- qplot() + geom_histogram(aes(x = SurfaceNoiseScore), data = ClusterMetrics, binwidth = 0.05) +
    ggtitle(label=Maintitle)
  Filename.pdf <- paste0(DataPath.mf_Quality, 'refFrag_', FragIndex, '/refFrag', FragIndex, 
                         '_SurfaceNoiseDist.pdf')
  pdf(file=Filename.pdf)
  print(ScoreHist)
  dev.off()
  
  return(ClusterMetrics)
}

##########################################################################################
## This function reads the straight score file
##########################################################################################
fn_readStraighScoreFile <- function(DataPath.mf_Quality, FragIndex, ...){
  Folderpath_Quality <- paste0(DataPath.mf_Quality, 'refFrag_', FragIndex, '/')
  Filename_Straight <- paste0(Folderpath_Quality, 'refFrag', FragIndex, '_StraightScore.txt')
  File_Straight <- try(read.table(Filename_Straight, sep = ' ', header = T, stringsAsFactors = F))
  File_Straight$X <- NULL
  return(File_Straight)  
}

##########################################################################################
## This function reads the Cluster metrics dataset
##########################################################################################
fn_readClusterMetrics <- function(DataPath.mf_Quality, FragIndex, ...){
  Folderpath_Quality <- paste0(DataPath.mf_Quality, 'refFrag_', FragIndex, '/')
  Filename <- paste0(Folderpath_Quality, 'refFrag', FragIndex, '_SurfaceNoiseScore.txt')
  ClusterMetrics <- try(read.table(Filename, sep = '\t', header = T, stringsAsFactors = F))
  ClusterMetrics$QualityScore <- sqrt(ClusterMetrics$SurfaceNoiseScore * ClusterMetrics$StraightScore)
  return(ClusterMetrics)  
}

##########################################################################################
## This function returns molecule list based on some cutoff, for a given fragment
##########################################################################################
fn_returnGoodMolecules <- function(Molecules, FragIndex, DataPath.mf_Quality, 
                                   QScore_Cutoff = NULL, QScore_Qtl = NULL){
  
  if(is.null(QScore_Cutoff) & is.null(QScore_Qtl)){
    stop('Both QScore_Cutoff & QScore_Qtl are NULL')
  }
  
  ClusterMetrics <- fn_readClusterMetrics(DataPath.mf_Quality, FragIndex)
  ClusterMetrics$MoleculeUsed <- ClusterMetrics$MoleculeID %in% Molecules
  #table(ClusterMetrics$Discard, ClusterMetrics$MoleculeUsed)
  ClusterMetrics_Used <- subset(ClusterMetrics, MoleculeUsed == TRUE)
  
  if(is.null(QScore_Cutoff)){
    QScore_Cutoff <- quantile(x = ClusterMetrics$QualityScore, probs = (1 - QScore_Qtl))
  }
  P1 <- qplot() + geom_histogram(aes(x = QualityScore), data = ClusterMetrics, binwidth = 0.05) +
    ggtitle(label='Quality Score of all molecules') + xlim(c(0, 1.05)) + 
    geom_vline(xintercept = QScore_Cutoff, col = 'blue', size = 1)
  P2 <- qplot() + geom_histogram(aes(x = QualityScore), data = ClusterMetrics_Used, binwidth = 0.05) +
    ggtitle(label='Quality Score of molecules within Stretch limit') + xlim(c(0, 1.05)) +
    geom_vline(xintercept = QScore_Cutoff, col = 'blue', size = 1)
  
  P3 <- qplot(x = MoleculeUsed, y = SurfaceNoiseScore, data = ClusterMetrics) + 
    geom_boxplot(aes(fill = MoleculeUsed)) + ggtitle(label=paste('ref Fragment', FragIndex)) +
    geom_jitter() + ylab(label='Surface Noise Score') +
    theme(legend.position="top")
  P4 <- qplot(x = MoleculeUsed, y = QualityScore, data = ClusterMetrics) + 
    geom_boxplot(aes(fill = MoleculeUsed)) +
    geom_jitter() + ylab(label='Quality Score') +
    theme(legend.position="top")
  
  Plot1234 <- arrangeGrob(P1, P3, P2, P4, ncol = 2)
  Molecules_Select <- subset(ClusterMetrics_Used, QualityScore >= QScore_Cutoff)[,'MoleculeID']
  return(list(Molecules_Select = Molecules_Select, QScore_Cutoff = QScore_Cutoff, 
              QScore_Plot = Plot1234))
}


