#' Spatial Classification Trees
#'
#' @description Spatial Classification Trees using Ripley's intertype K-function to partition bivariate marked spatial point process, with pruning.
#'
#' @param ypp bivariate marked spatial point process of class 'ppp' (see spatstat).
#' @param r initial resolution scale to use to construct spatcart maximal tree. Decreasing with splitting..
#' @param method optional. Parameter method = c("deviance", "misclass") sets the pruning method to use. "deviance" is related to class probability trees, while "misclass" is related to classical misclassification trees. Default is method = "deviance".
#' @param ties logical. If TRUE (default), in case of ties in splitting rule, the first split achieving maximum loss in impurity is taken. If FALSE, scale adaptation inside nodes is made to avoid ties.
#' @param offset parameter related to the tree() function of package tree.
#' @param wt vector of non-negative observational weights; fractional weights are allowed.
#' @param parms other parameters related to the tree function of package tree.
#' @param minsplit minimal size for the node to be split. Default is minsplit = 10.
#' @param minleaf minimal size of each resulting leaf for a node to be split. Default is minleaf = 5.
#' @param mindev minimal value of node impurety to control node splitting. Default is mindev = 0.
#' @param graph Logical. If TRUE (Default), plot the heatmaps associated to the point process y and optimal trees produced by spatcart.
#'
#' @return a list of objects containing the maximal, largest and smallest optimal trees; the sequence of pruned subtrees; the values of scale resolution and impurity inside nodes of the maximal tree.
#'
#' @import dplyr
#' @import tidyr
#' @import spatstat
#' @import tree
#' @import ggplot2
#'
#' @details parameter minsplit must be at least 2xminleaf.
#'
#' @references *Spatial Classification Trees*, by A. Bar-Hen, S. Gey and J.-M. Poggi (2020).
#'
#' @export
#'
#' @examples
#' library(spatcart)
#' library(spatstat)
#' # Simulate bivariate marked spatial point process
#' chess = damier(500, h=0.45, model="Poisson")
#' ypp = ppp(x=chess$data$x1,
#' y=chess$data$x2,
#' marks = chess$data$label,
#' range(chess$data$x1),
#' range(chess$data$x2))
#'
#' # Checking initial resolution to construct maximal tree via Ripley's Kcross
#' major =  names(which.max(intensity(ypp)))
#' minor = names(which.min(intensity(ypp)))
#' K0=Kcross(ypp, major, minor, correction="none")
#'
#' K01 = data.frame(
#'   r = K0$r,
#'   difference = K0$un-K0$theo
#' )
#'
#' library(ggplot2)
#' DiffK = ggplot(K01, aes(x=r, y=difference))+geom_line()+xlab("r")+
#'   ylab("Kest - Ktheo")+
#'   ggtitle("Difference between estimated and theoretical K")
#' DiffK
#'
#' r = K0$r[which.min(K0$un-K0$theo)]
#'
#' # Keeping ties in splitting rules
#' t = spatcart(ypp, r)
#'
#' # Optimizing ties in splitting rules
#' t = spatcart(ypp, r, ties = FALSE)
#'
#' # Setting stopping rule for splitting
#' t = spatcart(ypp, r, minsplit = 50, minleaf = 25)
spatcart = function(ypp, r, method = "deviance", ties = TRUE, offset=NULL, wt=NULL, parms=NULL, minsplit=10, minleaf=5, mindev=0, graph = TRUE){
  if (!is.ppp(ypp))
    stop("data must be of class ppp")

  # Maximal Tree
  spatt = ypp %>% spattree(r, ties, offset, wt, parms, minsplit, minleaf, mindev)

  # Pruning
  if(sum(spatt$tree$frame$var=="<leaf>")>1){
    class(spatt$tree)="tree"
    seq = spatt$tree %>% prune.tree(method=method)

    nbleav = choice.tree(seq)
    t1 = spatt$tree %>% prune.tree(method=method, best=nbleav$maxgap)
    t2 = spatt$tree %>% prune.tree(method=method, best=nbleav$plateau)

  } else {
    t1 = spatt$tree
    t2 = spatt$tree
  }

  if(graph){
    t = spatt

    #### Maximal tree
    a = t$tree
    K = t$K[rownames(a$frame)][a$frame$var=="<leaf>"]

    pmax = partition.spattree(a, ypp, K) +
      ggtitle("SpatCART maximal tree")+
      theme(plot.title = element_text(size=16))+
      theme(legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))


    #### Largest optimal tree
    a = t1
    K = t$K[rownames(a$frame)][a$frame$var=="<leaf>"]

    p1 = partition.spattree(a, ypp, K) +
      ggtitle("SpatCART largest optimal tree")+
      theme(plot.title = element_text(size=16))+
      theme(legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))

    #### Smallest optimal tree
    a = t2
    K = t$K[rownames(a$frame)][a$frame$var=="<leaf>"]

    p2 = partition.spattree(a, ypp, K) +
      ggtitle("SpatCART smallest optimal tree")+
      theme(plot.title = element_text(size=16))+
      theme(legend.title = element_text(size = 14),
            legend.text = element_text(size = 12))

    print(pmax)
    print(p1)
    print(p2)
  }
  list(max.tree = spatt$tree, pruned.seq = seq, opt.tree.max=t1, opt.tree.min = t2, cp=spatt$cp, K=spatt$K)
}
