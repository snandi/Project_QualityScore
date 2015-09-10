## This function library will contain functions exclusively used for the Quality Score project

fn_ClusterAndPlot <- function(PixelData, Molecule){
  PixelData_Norm_Long <- melt(data = PixelData[,c('MoleculeID', 'Pixel1_Norm', 'Pixel2_Norm', 'Pixel3_Norm')], 
                              id.vars = 'MoleculeID', 
                              measure.vars = c('Pixel1_Norm', 'Pixel2_Norm', 'Pixel3_Norm'))
  Model3 <- Mclust(data = log(PixelData[,'Pixel3_Norm']), G=1:4)
  NClust3 <- Model3$G
  Model2 <- Mclust(data = log(PixelData[,'Pixel2_Norm']), G=1:4)
  NClust2 <- Model2$G
  Model1 <- Mclust(data = log(PixelData[,'Pixel1_Norm']), G=1:4)
  NClust1 <- Model1$G
  
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
