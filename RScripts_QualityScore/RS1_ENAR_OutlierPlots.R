rm(list = ls(all.names = TRUE))
rm(list = objects(all.names = TRUE))
#dev.off()

library(ggplot2)

Y1 <- rgamma(n = 100, shape = 1, rate = 0.5)
Y2 <- rgamma(n = 100, shape = 1, rate = 1)

Data <- as.data.frame(rbind(cbind(Y = Y1, Num = 1), cbind(Y = Y2, Num = 2)))

set.seed(100)
Plot1 <- qplot() + geom_boxplot(aes(x = factor(Num), y = Y), data = Data, outlier.size = 5, fill = 'gray70', color = 'red') +
  xlab(label = '') + ylab(label = '') +
  coord_flip()
Filename1 <- paste0('/ua/snandi/Seminar_Conference_Workshop/ENAR_2016/Presentation/Plots/', 'Plot1.bmp')
bmp(filename = Filename1)
plot(Plot1)
dev.off()

library(MASS)
Sigma <- matrix(data = c(1, 0.8, 0.8, 2), nrow = 2)
set.seed(105)
XY <- mvrnorm(n = 50, mu = c(0, 0.5), Sigma = Sigma)

XY <- as.data.frame(XY)
XY$Dist <- 0
Center <- c(0, 0.5)
for(i in 1:nrow(XY)){
  X <- XY[i,]
  Y <- Center
  XY[i,'Dist'] <- dist(rbind(X, Y), method = 'euclidean')
}

#qplot() + geom_point(aes(x = V1, y = V2, size = 1/Dist), data = XY, col = 'red') +
Plot2 <- qplot() + geom_point(aes(x = V1, y = V2), data = XY, col = 'red', size = 5) +
  xlab(label = 'x') + ylab(label = 'y') +
  theme(legend.position = '') +
  geom_point(aes(x = 0, y = 0.5), col = 'black', size = 7) 

Filename1 <- paste0('/ua/snandi/Seminar_Conference_Workshop/ENAR_2016/Presentation/Plots/', 'Plot2.bmp')
bmp(filename = Filename1)
plot(Plot2)
dev.off()

Plot3 <- qplot() + geom_point(aes(x = V1, y = V2), data = XY, col = 'red', size = 4) +
  xlab(label = 'x') + ylab(label = 'y') +
  theme(legend.position = '') +
  geom_point(aes(x = 0, y = 0.5), col = 'black', size = 7, shape = 18, fill = 'black') +
  stat_ellipse(aes(x = V1, y = V2), data = XY, geom = 'polygon', level = 0.8, alpha = 0.5) +
  stat_ellipse(aes(x = V1, y = V2), data = XY, geom = 'polygon', level = 0.95, alpha = 0.3)
Filename1 <- paste0('/ua/snandi/Seminar_Conference_Workshop/ENAR_2016/Presentation/Plots/', 'Plot3.bmp')
bmp(filename = Filename1)
plot(Plot3)
dev.off()

