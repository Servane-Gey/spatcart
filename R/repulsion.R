#' Repulsive simulated data set
#'
#' @description Simulation of a locally repulsive spatial marked point process
#'
#' @param n intensity of the Poisson point process used to simulate data (see spatstat).
#' @param r value of the local repulsive scale. 0< r <1.
#' @param h margin parameter defining the mixing proportion of marks inside window. Default is h = 0.4.
#' @param ncase number of marks of the marked point process. Default is ncase = 2.
#' @param graph logical. If TRUE (Default), plot the data scatter plot.
#'
#' @return a data.frame object with variables (x, y, marks).
#'
#' @import spatstat
#' @importFrom  graphics points
#' @importFrom graphics lines
#' @importFrom  stats runif
#' @importFrom stats rbinom
#'
#' @export
#'
#' @examples
#' repulsion(1000, 0.05)
#' repulsion(1000, 0.1, h = 0.1)
repulsion = function(n, r, h=0.4, ncase = 2, graph = TRUE){

  data=damier(n, model = "Poisson", h=h, d=(ncase-1), graph = FALSE)$data
  data.repuls =  data
  index.x1 = findInterval(data$x1,seq(0,1,length=ncase+1))

  for (i_ in seq(1,ncase, by=2)) {
    ##obs de la case
    data.sub=data[index.x1==i_,]
    ##et on decoupe en deux
    minoritaire =  names(which.min(table(data.sub$label)))
    majoritaire =  names(which.max(table(data.sub$label)))

    minoritaire = data.sub[data.sub$label==minoritaire,]
    majoritaire = data.sub[data.sub$label==majoritaire,]

    ##et on ne garde que les majoritaires qui sont trop pres d'un minoritaire
    proche = 0
    if (nrow(minoritaire)>0){ ##si un seul label
      dist = crossdist(minoritaire$x1, minoritaire$x2, majoritaire$x1, majoritaire$x2)
      proche = apply(dist,2,function(k_)min(k_)<r)
    }

    k=0 ##pour limiter le nombre de boucles
    while ( (sum(proche)>0)&(k<100) ){
      k = k+1
      dist = crossdist(minoritaire$x1, minoritaire$x2, majoritaire$x1, majoritaire$x2)
      proche = apply(dist,2,function(k_)min(k_)<r)
      majoritaire$x1[proche] = runif(sum(proche),min=(i_-1)/ncase,max=i_/ncase) ##on simule les x
      majoritaire$x2[proche] = runif(sum(proche)) ##on simule les y
    }

    data.repuls[index.x1==i_,] = rbind(majoritaire,minoritaire)

  }
  if(graph == TRUE){
    class1=data.repuls[data.repuls$label=="blue",]
    class2=data.repuls[data.repuls$label=="red",]

    plot(0:1,0:1,type="n", xlab="", ylab="")
    points(class1[,2],class1[,3], pch=19, cex = 0.4, col="blue")
    points(class2[,2],class2[,3], pch=21, cex = 0.4, col="red")
    for(i in 1:(ncase-1)){
      lines(c(i/ncase,i/ncase),c(0,1))
    }
  }
  return(data.repuls)
}
