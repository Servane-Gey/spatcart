#' Heatmap of interaction between marks
#'
#' @description Plot the heatmap induced by a tree partition of interaction between marks for a bivariate marked spatial point process
#'
#' @param a an object of class "tree". Currently provided by the functions spatcart or spattree.
#' @param ypp the bivariate marked spatial point process used to generate tree a, of class 'ppp' (see spatstat).
#' @param K the values of interactions between marks of the point process ypp restricted to the leaves of tree a. Interactions are computed through the Ripley's Kcross function. They can be found in output values of functions spatcart or spattree.
#' @param d a scalar parameter defining the dxd grid resolution to produce the heatmap associated to the leaves of tree a and values of K. Default value is d = 100, leading to a 10e+4 grid of points.
#'
#' @return an object of class "ggplot".
#'
#' @references *Spatial Classification Trees*, by A. Bar-Hen, S. Gey and J.-M. Poggi (2020).
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' library(spatcart)
#' # Loading Paracou data set
#' ypp = Paracoudata()
#' r = 15
#' t = spatcart(ypp, r, graph = FALSE)
#'
#' a = t$opt.tree.min
#' K = t$K[rownames(a$frame)][a$frame$var=="<leaf>"]
#'
#' partition.spattree(a, ypp, K)
partition.spattree = function(a, ypp, K, d = 100){

  # Initialize global variables
  pointsx = pointsy = Mark = NULL
  x = y = z = Interaction = NULL
  nuage = tibble(
    pointsx = ypp$x,
    pointsy = ypp$y,
    Mark = ypp$marks
  )
  grille = make.grid(a, ypp, K, d)$grid


  if(length(K)>1){
    p = ggplot(grille, aes(x,y, z=Interaction)) + geom_tile(aes(fill = Interaction)) +
      scale_fill_gradient(low="white", high="grey") +
      geom_point(data = nuage, mapping = aes(pointsx, pointsy, shape = Mark), inherit.aes = FALSE) +
      scale_shape_manual(breaks = levels(ypp$marks), values=c(19, 0)) +
      theme_bw() + xlab("") + ylab("") +
      gg.partition.tree(a, label="marks", color = "black", ordvars= c("x","y"))
  } else {
    p = ggplot(grille, aes(x,y, z=Interaction)) + geom_tile(aes(fill = Interaction)) +
      scale_fill_gradient(low="white", high="grey") +
      geom_point(data = nuage, mapping = aes(pointsx, pointsy, shape = Mark), inherit.aes = FALSE) +
      scale_shape_manual(breaks = levels(ypp$marks), values=c(19, 0)) +
      theme_bw() + xlab("") + ylab("")
  }
  return(p)
}
