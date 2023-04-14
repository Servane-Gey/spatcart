#' Chess simulated data set
#'
#' @description Simulation of a bivariate marked spatial point process with marks distribued as a chessboard.
#'
#' @param n if model = "unif", number of observations. If model = "Poisson", intensity of the Poisson point process used to simulate data (see spatstat).
#' @param model nominal, model = c("unif", "Poisson"). If model = "unif", generates points uniformly on the window. If model = "Poisson", generates points as a Poisson point process in the window. Default is model = "unif".
#' @param h margin parameter defining the mixing proportion of marks inside windows. Default is h = 0.4.
#' @param d set the number of squares in the chessboard as dxd. Default is d = 3.
#' @param graph logical. If TRUE (default), plot the data scatter plot and the chessboard's squares.
#'
#' @return a list of data.frame objects with data points, marks distribution and Bayes predictions.
#'
#' @import spatstat
#' @importFrom graphics points
#' @importFrom graphics lines
#' @importFrom stats runif
#' @importFrom stats rbinom
#'
#' @export
#'
#' @examples
#' library(spatstat)
#' damier(1000, model = "unif")
#' damier(1000, model = "unif", h = 0.1)
#' damier(1000, model = "Poisson")
damier=function(n, model = "unif", h=0.4,d=3, graph=TRUE){

  stopifnot(d/2-floor(d/2)>0)

  if(model=="unif"){
    x=runif(2*n)
    x=matrix(x,nrow=n)
  }
  if(model=="Poisson"){
    x=rpoispp(n, win=square(1))
    x=data.frame(x)
  }

  nx=nrow(x)
  p = numeric(nx)
  p[x[,1]<1/d & x[,2]<1/d]=1/2+h

  l=numeric(d^2)
  l=matrix(l,nrow=d)
  l[1,1]=1

  if(d>1){
    for(j in 1:(d-1)){
      if(l[j,1]==1){p[x[,1]>=j/d & x[,1]<(j+1)/d & x[,2]<1/d]=1/2-h
      l[j+1,1]=0}
      else {p[x[,1]>=j/d & x[,1]<(j+1)/d & x[,2]<1/d]=1/2+h
      l[j+1,1]=1}
    }

    for(i in 1:d){
      for(k in 2:d){
        if(l[i,k-1]==1){ p[x[,2]>=(k-1)/d & x[,2]<k/d & x[,1]>=(i-1)/d & x[,1]<i/d]=1/2-h
        l[i,k]=0}
        else {p[x[,2]>=(k-1)/d & x[,2]<k/d & x[,1]>=(i-1)/d & x[,1]<i/d]=1/2+h
        l[i,k]=1}
      }
    }
  }

  y=rbinom(nx,size=1,prob=p)
  bayes = integer(nx)
  bayes[p>=1/2]=1
  bayes[p<1/2]=0

  data = data.frame(cbind(y,x))
  colnames(data)=c("label","x1","x2")
  data$label=factor(data$label,levels=0:1,labels=c("blue","red"))
  bayes=factor(bayes,levels=0:1,labels=c("blue","red"))

  # graph ------------
  if(graph){
    class1=data[data$label=="blue",]
    class2=data[data$label=="red",]

    plot(0:1,0:1,type="n", xlab="", ylab="")
    points(class1[,2],class1[,3],pch=19, cex = 0.4, col="blue")
    points(class2[,2],class2[,3],pch=21, cex = 0.4, col="red")
    for(i in 1:(d-1)){
      lines(c(i/d,i/d),c(0,1))
      lines(c(0,1),c(i/d,i/d))
    }
  }
  sortie=list(data, p, bayes)
  names(sortie)=c("data", "prob", "bayes")
  return(sortie)
}
