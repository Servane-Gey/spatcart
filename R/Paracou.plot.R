#' Heatmap of interaction between species of the Paracou data set
#'
#' @description  Plot the heatmap of interaction between marks with altitudes information for the Paracou data set.
#'
#' @param a an object of class "tree". Currently provided by the functions spatcart or spattree.
#' @param ypp the Paracou bivariate marked spatial point process used to generate tree a.
#' @param K the values of interactions between marks of the point process ypp restricted to the leaves of tree a. Interactions are computed through the Ripley's Kcross function. They can be found in output values of functions spatcart or spattree.
#' @param d a scalar parameter defining the dxd grid resolution to produce the heatmap associated to the leaves of tree a and values of K. Default value is d = 100, leading to a 10e+4 grid of points.
#'
#' @return an object of class "ggplot".
#'
#' @references *Spatial Classification Trees*, by A. Bar-Hen, S. Gey and J.-M. Poggi (2020).
#'
#' @import spatstat
#' @import dplyr
#' @import tidyr
#' @import reshape2
#' @import metR
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' library(spatcart)
#' # Loading Paracou bivariate spatial point process
#' ypp = Paracoudata()
#'
#' # SpatCART trees
#' t = spatcart(ypp, 15, graph = FALSE)
#'
#' a = t$opt.tree.min
#' K = t$K[rownames(a$frame)][a$frame$var=="<leaf>"]
#'
#' Paracou.plot(a,ypp,K)
Paracou.plot = function(a, ypp, K, d = 100){
  # Initialize global variables for package
  pointsx = pointsy = Specie = NULL
  x = y = z = Interaction = NULL

  grille = make.grid(a,ypp,K,d)$grid
  nuage = tibble(
    pointsx = ypp$x,
    pointsy = ypp$y,
    Specie = ypp$marks
  )

  if(length(K)>1){
    p = ggplot(grille, aes(x,y, z=Interaction)) + geom_tile(aes(fill = Interaction)) +
      scale_fill_gradient(low="white", high="grey")+
      geom_point(data = nuage, mapping = aes(pointsx, pointsy, shape = Specie), inherit.aes = FALSE) +
      scale_shape_manual(breaks = c("Oxandra", "Vouacapoua"), values=c(19, 0)) +
      theme_bw() + xlab("") + ylab("") +
      gg.partition.tree(a, label="marks", color = "black", ordvars= c("x","y"))
  } else {
    p = ggplot(grille, aes(x,y, z=Interaction)) + geom_tile(aes(fill = Interaction)) +
      scale_fill_gradient(low="white", high="grey")+
      geom_point(data = nuage, mapping = aes(pointsx, pointsy, shape = Specie), inherit.aes = FALSE) +
      scale_shape_manual(breaks = c("Oxandra", "Vouacapoua"), values=c(19, 0)) +
      theme_bw() + xlab("") + ylab("")
  }
  p = p+stat_contour(data=plot1.melt, aes(x = (x-1)*5, y = (y-1)*5, z = z), color="black")
  p = p+geom_text_contour(data=plot1.melt, aes(x = (x-1)*5, y = (y-1)*5, z=z))
  return(p)
}
