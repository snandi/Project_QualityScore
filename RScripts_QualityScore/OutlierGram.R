##############################################################
## This was downloaded from https://github.com/dapr12/Classifficationfda/blob/master/Classiffda/R/outgram.R
##############################################################

#' OutlierGram - OutlierGram 
#' @param x: matrix with n rows and p columns. Each row represents a curve.
#' @param t: vector of length p with the time points at which curves are measured.
#' @param mag: whether to look for magnitude outliers (mag=T) or not (mag=F).
#' @param factorsh: value of the factor for the detection of shape outlier.
#' @param factormg: value of the factor for the detection of magnitude outlier.
#' @param p1: at each time point, curves with values in the p1 percentile while be shifted downwards and the (mbd,mei)
#' will be recalculated in the new position.
#' @param p2: at each time point, curves with values in the p2 percentile while be shifted upwards and the (mbd,mei)
#' will be recalculated in the new position.
#' @param plotting: whether to produce a plot with the observed curves and detected outliers and the outliergram 
#' (plotting=T) or not (plotting=F).
#' @param title: main title for the graphics with the observed curves and detected outliers.
#' @param ylabel: ylabel for the graphics with the observed curves and detected outliers.
#' @param linecol: color for line in graphics 1= gray and black, 2= multicolor.
#' @param legend: whether to produce a legend in the graphics with the observed curves and detected outliers
#' to distinguish shape and magnitude outliers.    
#' @return The component list is: 
#'\itemize{
#'  \item{ modified.band.depth}{}
#'  \item{ modified.epigraph.index}{}
#'  \item{ dist}{}
#'  \item{ magnitude.outliers}{}
#'  \item{ shape.outliers}{}
#'  \item{ outliers }
#'}
#' @references {Shape Outlier Detection and Visualisation for Functional Data: the Outliergram. Ana Arribas-Gil and Juan Romo. 
#' Biostatistics, (DOI) 10.1093/biostatistics/kxu006, 2014}
#' @examples \dontrun{
#' ## OutlierGram
#' fdaobjMale<-fdaclass(Gwd$Male,Gwd$Age,c(1,18))
#' OG<-outliergram(fdaobjMale)
#' #' OG$shape.outliers
#' OG$magnitude.outliers
#' OG$outliers
#' NewMat<-Mat[,-OG$outliers]
#' matplot(time, NewMat,type="l",col =color,xlab="",ylab="x(t)",main="Observations",lty=1)
#' }
#' @export
outliergram<- function(fdaobj, mag=T, factorsh=1.5, factormg=1.5, p1=1, p2=0, ...)
{
  
  if (!(inherits(fdaobj, "FdaClass"))) 
    stop("Argument FD  not a functional data object.")
  
  data<- fdaobj$data
  
  argvals<- fdaobj$argvals
  
  rangevals<- fdaobj$rangeval
  
  time<- argvals
  
  x<- as.matrix(t(data))
  
  n <- nrow(x)
  
  p <- ncol(x) 
  
  if (n>1) {
    
    if (p == 1) {x <- t(x)}   
    
    rmat=apply(t(x),1,rank,ties.method="max") 
    
    down=rmat-1
    
    up=n-rmat
    
    mbd=(rowSums(up*down)/p+n-1)/combinat(n,2)
    
    epi=rowSums(up+1)/(n*p)
    
    s.epi<-sort(epi,index.return=T)
    
    if (mag) { 
      
      s.mbd<-sort(mbd,index.return=T)
      
      m=ceiling(n*0.5)
      
      center=x[s.mbd$ix[(m+1):n],]
      
      inf=apply(center,2,min)
      
      sup=apply(center,2,max)
      
      dist=factormg*(sup-inf)
      
      upper=sup+dist
      
      lower=inf-dist
      
      upper<-rep(upper,n)
      
      dim(upper)<-c(p,n)
      
      upper<-t(upper)
      
      lower<-rep(lower,n)
      
      dim(lower)<-c(p,n)
      
      lower<-t(lower)
      
      outly=(x<=lower)+(x>=upper)
      
      outrow=rowSums(outly)
      
      mag.out<-which(as.vector(outrow>0))
    }
    
    else{mag.out<-c()}
    
    a0=-2/(n*(n-1))
    
    a1=2*(n+1)/(n-1)
    
    a2=a0
    
    P=a0+a1*epi+a2*n^2*epi^2
    
    d<-P-mbd
    
    d.sum<-fivenum(d)
    
    limit<-d.sum[4]+factorsh*(d.sum[4]-d.sum[2])
    
    sh.out<-which(as.vector(d>=limit))
    
    x.or<-x
    
    sh.out2<-c()
    
    me.mb<-c()
    
    non.out<-setdiff(1:n,sh.out)
    
    for (i in non.out) {
      
      x<-x.or
      
      if (epi[i]<0.5){         
        
        s<-sort(x[i,]-apply(x[-i,,drop=F],2,quantile,p=p1),index.return=T)   
        
        if (s$x[p] >0){
          
          x[i,]<-x[i,]-s$x[p]  
          
          rmat=apply(t(x),1,rank,ties.method="max")  
          
          down=rmat-1
          
          up=n-rmat
          
          mbd2=(rowSums(up*down)/p+n-1)/combinat(n,2)
          
          epi2=rowSums(up+1)/(n*p)
          
          s.epi2=sort(epi2,index.return=T)         
          
          if (mbd2[i]<(a0+a1*epi2[i]+a2*n^2*epi2[i]^2)-limit) {
            
            sh.out2<-append(sh.out2,i)
            
            me.mb<-cbind(me.mb,c(epi2[i],mbd2[i]))
            
            d[i]<-mbd2[i]-(a0+a1*epi2[i]+a2*n^2*epi2[i]^2)
            
          }
        }
      }
      
      if (epi[i]>=0.5){
        
        s<-sort(x[i,]-apply(x[-i,,drop=F],2,quantile,p=p2),index.return=T)   
        
        if (s$x[1] <0){
          
          x[i,]<-x[i,]-s$x[1]  
          
          rmat=apply(t(x),1,rank,ties.method="max") 
          
          down=rmat-1
          
          up=n-rmat
          
          mbd2=(rowSums(up*down)/p+n-1)/combinat(n,2)         
          
          epi2=rowSums(up+1)/(n*p)
          
          s.epi2=sort(epi2,index.return=T)        
          
          if (mbd2[i]<(a0+a1*epi2[i]+a2*n^2*epi2[i]^2)-limit) {
            
            sh.out2<-append(sh.out2,i)
            
            me.mb<-cbind(me.mb,c(epi2[i],mbd2[i]))
            
            d[i]<-mbd2[i]-(a0+a1*epi2[i]+a2*n^2*epi2[i]^2)
            
          }
        }
      }
    }
    
    sh.out<-sort(append(sh.out,sh.out2))
    
    d.sh<-d[sh.out]
    
    sd.sh<-sort(d.sh,decreasing=T,index.return=T)
    
    sh.out<-sh.out[sd.sh$ix]
    
    x<-x.or
    
    title="Observations"
    
    ylabel="x(t)"
    
    linecol=2
    
    type="l"
    
    lty=2
    
    lwd=1
    
    cex=1
    
    dev.new(width=12,height=6)
    
    layout(matrix(c(1,2), 1, 2, byrow = T))
    
    if (length(t)<1) {t<-1:p}
    
    matplot(time,t(x),type="l",xlab="",ylab=ylabel,main=title,lty=1)
    
    for (ou in mag.out) {
      lines(time,x[ou,],col=1,lty=2)      
    }
    for (ou in sh.out) {
      lines(time,x[ou,],col=1,lty=1)      
    }
    
    
    plot(epi,mbd,type="n",xlab="Modified Epigrahp Index",ylab="Modified Band Depth",main="Outliergram")
    
    text(epi,mbd)
    
    for (i in sh.out2) {
      
      j<-which(sh.out2==i)
      
      points(me.mb[1,j],me.mb[2,j],col ="red",pch=20,cex=3)
      
      text(me.mb[1,j],me.mb[2,j])
      
    }     
    
    Ps=a0+a1*s.epi$x+a2*n^2*s.epi$x^2
    
    lines(s.epi$x,Ps,col="black")     
    
    lines(s.epi$x,Ps-limit,col="black",lty=2)
    
    #} 
    
  }
  return(list(modified.band.depth=mbd,
              modified.epigrhap.index=epi,
              dist=d,
              magnitude.outliers=mag.out,
              shape.outliers=sh.out,
              outliers=sort(unique(c(sh.out,mag.out))) ) )
}     

combinat<-function(n,p){
  
  if (n<p){combinat=0}
  
  else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
  
}
