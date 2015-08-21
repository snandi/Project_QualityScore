library("pdfCluster")
data("wine", package = "pdfCluster")
winesub <- wine[, c(2, 5, 8)]
winevector <- c(winesub$Alcohol, winesub$Alcalinity, winesub$Flavanoids)

ggplot() + geom_histogram(aes(x = winevector, y = ..density..), fill='gray50') + 
  geom_density(aes(x = winevector))

## pdf.ke: Kernel estimate of a pdf
pdf.ke <- pdfCluster::kepdf(x=winesub)
plot.(pdf.ke)

## Clustering
cl.winesub <- pdfCluster::pdfCluster(x=winesub)
summary(cl.winesub)

#plot(cl.winesub)


cl.Pixels <- pdfCluster::pdfCluster(x = PP, graphtype='unidimensional', bwtype='adaptive')
summary(cl.Pixels)
#plot(cl.Pixels)
PP.group <- groups(cl.Pixels)
PP.1 <- PP[which(PP.group==1)]
PP.2 <- PP[which(PP.group==2)]

pdf.ke <- pdfCluster::kepdf(PP)
plot(pdf.ke)

Hist_Dens <- ggplot(aes(x = PP), data=as.data.frame(PP)) + 
  geom_histogram(aes(x = PP, y = ..density..), fill='gray60') + 
  geom_density(aes(x = PP), kernel = 'epanechnikov', col = 'gray20', lwd = 1) + 
  stat_function(fun = dgamma, args=c(shape=Shape, rate=Rate), col='red', size = 1) +
  geom_line(aes(x = pdf.ke@x, y = pdf.ke@estimate), col = 'royalblue1', size = 1)

#   geom_density(aes(x = PP.1), kernel = 'epanechnikov', col = 'red', lwd=1, lty=2) + 
#   geom_density(aes(x = PP.2), kernel = 'epanechnikov', col = 'red', lwd=1, lty=2)

Hist_Dens
