#' CART and SpatCART on the Paracou data set
#'
#' @description Produce graphical results of CART and SpatCART algorithms at fixed initial resolution for the Paracou data set.
#'
#' @param r the initial resolution used by SpatCART to construct the maximal tree (see spattree).
#' @param minsplit minimal size for the node to be split. Default is minsplit = 10 and.
#' @param minleaf minimal size of each resulting leaf for a node to be split. Default is minleaf = 5.
#' @param graph logical. If TRUE (Default), also produce the graphs of the number of leaves of the pruned subtrees sequence with respect to the complexity parameter used in the pruning algorithm. The triangle symbol represents the tree selected via the modified largest jump method slpe heuristics, while the diamond symbol represents the tree selected via the largest plateau method one.
#'
#' @return a list of 7 objects of class "ggplot".
#'
#' @references *Spatial Classification Trees*, by A. Bar-Hen, S. Gey and J.-M. Poggi, HAL 01837065 (2018).
#'
#' @import spatstat
#' @import tree
#' @import dplyr
#' @import tidyr
#' @import reshape2
#' @import metR
#' @import ggplot2
#' @importFrom graphics lines
#' @importFrom stats stepfun
#'
#' @export
#'
#' @examples
#' Paracou(15)
Paracou = function(r, minsplit = 10, minleaf = 5, graph = TRUE){
  ypp = Paracou.pp

  #################################################################
  # Ripley's Kcross
  #################################################################
  # Initialize global variables for package
  scale = difference = NULL
  major =  names(which.max(intensity(ypp)))
  minor = names(which.min(intensity(ypp)))
  K0=Kcross(ypp, major, minor, correction="none")

  K01 = tibble(
    scale = K0$r,
    difference = K0$un-K0$theo
    )
  DiffK = ggplot(K01, aes(x=scale, y=difference))+geom_line()+xlab("r")+
    ylab("Kest - Ktheo")+
    ggtitle("Difference between estimated and theoretical K")+
    geom_vline(xintercept = r, colour = "blue", linetype = "dashed")+
    theme(axis.title.x = element_text(size=16))+
    theme(axis.title.y = element_text(size=16))+
    theme(plot.title = element_text(hjust = 0.5, size=16))

  #################################################################
  # SpatCART
  #################################################################


  # Maximal tree
  t = spattree(ypp, r, minsplit=minsplit,minleaf=minleaf)

  a = t$tree
  K = t$K[rownames(a$frame)][a$frame$var=="<leaf>"]
  p1 = Paracou.plot(a,ypp,K)+
    ggtitle("Paracou - SpatCART maximal tree")+
    theme(plot.title = element_text(size=16))+
    theme(legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))

  # Smallest probability tree
  a = spatprune(t$tree)$opt.tree.min
  K = t$K[rownames(a$frame)][a$frame$var=="<leaf>"]
  p2 = Paracou.plot(a,ypp,K)+
    ggtitle("Paracou - SpatCART class probability tree")+
    theme(plot.title = element_text(size=16))+
    theme(legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))

  # Smallest misclassification tree
  a = spatprune(t$tree, method = "misclass")$opt.tree.min
  K = t$K[rownames(a$frame)][a$frame$var=="<leaf>"]
  p3 = Paracou.plot(a,ypp,K)+
    ggtitle("Paracou - SpatCART classification tree")+
    theme(plot.title = element_text(size=16))+
    theme(legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))

  #################################################################
  # CART
  #################################################################

  # Maximal tree
  # Initialize global variables for package
  x = y = z = Specie = NULL
  nuage = data.frame(
    x = ypp$x,
    y = ypp$y,
    Specie = ypp$marks
    )
  colnames(nuage) = c("x","y", "Specie")
  t = tree(Specie~.,data=nuage,split="gini",model=T, minsize=minsplit,mincut=minleaf)

  a = t
  p4 = ggplot(nuage, aes(x,y)) + geom_point(aes(x=x, y=y, shape=Specie)) +
    scale_shape_manual(values=c(19, 0)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    xlab("") + ylab("") +
    gg.partition.tree(a, label="Specie", color = "black", ordvars= c("x","y")) +
    ggtitle("Paracou - CART maximal tree")+
    theme(plot.title = element_text(size=16))+
    theme(legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
  p4 = p4+stat_contour(data=plot1.melt, aes(x = (x-1)*5, y = (y-1)*5, z = z), color="black")
  p4 = p4+geom_text_contour(data=plot1.melt, aes(x = (x-1)*5, y = (y-1)*5, z=z))

  # Smallest probabilty tree
  seq_dev = prune.tree(t)
  nbleav = choice.tree(seq_dev)

  alpha.maxgap = seq_dev$k[seq_dev$size==nbleav$maxgap]
  alpha.plateau = seq_dev$k[seq_dev$size==nbleav$plateau]
  plot(stepfun(seq_dev$k[-1],seq_dev$size), main=c("CART class probability trees", "Nb leaves vs complexity"), xlab="Complexity", ylab = "Number of leaves")
  lines(c(alpha.maxgap,alpha.maxgap), c(0, nbleav$maxgap), type = "b", pch = 2, lty="dashed", col="blue")
  lines(c(-30,alpha.maxgap), c(nbleav$maxgap, nbleav$maxgap), type = "b", pch = 2, lty="dashed", col="blue")
  lines(c(alpha.plateau,alpha.plateau), c(0, nbleav$plateau), type = "b", pch = 8, lty="dotdash", col="red")
  lines(c(-30,alpha.plateau), c(nbleav$plateau, nbleav$plateau), type = "b", pch = 8, lty="dotdash", col="red")

  t_dev = prune.tree(t, best = nbleav$plateau)
  a = t_dev

  p5 = ggplot(nuage, aes(x,y)) + geom_point(aes(x=x, y=y, shape=Specie)) +
    scale_shape_manual(values=c(19, 0)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    xlab("") + ylab("") +
    gg.partition.tree(a, label="Specie", color = "black", ordvars= c("x","y")) +
    ggtitle("Paracou - CART class probability tree")+
    theme(plot.title = element_text(size=16))+
    theme(legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
  p5 = p5+stat_contour(data=plot1.melt, aes(x = (x-1)*5, y = (y-1)*5, z = z), color="black")
  p5 = p5+geom_text_contour(data=plot1.melt, aes(x = (x-1)*5, y = (y-1)*5, z=z))

  # Smallest misclassification tree
  seq_mcr = prune.misclass(t)
  nbleav = choice.tree(seq_mcr)

  alpha.maxgap = seq_mcr$k[seq_mcr$size==nbleav$maxgap]
  alpha.plateau = seq_mcr$k[seq_mcr$size==nbleav$plateau]
  plot(stepfun(seq_mcr$k[-1],seq_mcr$size), main=c("CART classification trees", "Nb leaves vs complexity"), xlab="Complexity", ylab = "Number of leaves")
  lines(c(alpha.maxgap,alpha.maxgap), c(0, nbleav$maxgap), type = "b", pch = 2, lty="dashed", col="blue")
  lines(c(-30,alpha.maxgap), c(nbleav$maxgap, nbleav$maxgap), type = "b", pch = 2, lty="dashed", col="blue")
  lines(c(alpha.plateau,alpha.plateau), c(0, nbleav$plateau), type = "b", pch = 8, lty="dotdash", col="red")
  lines(c(-30,alpha.plateau), c(nbleav$plateau, nbleav$plateau), type = "b", pch = 8, lty="dotdash", col="red")

  t_mcr = prune.misclass(t, best=nbleav$plateau)
  a = t_mcr

  p6 = ggplot(nuage, aes(x,y)) + geom_point(aes(x=x, y=y, shape=Specie)) +
    scale_shape_manual(values=c(19, 0)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    xlab("") + ylab("") +
    gg.partition.tree(a, label="Specie", color = "black", ordvars= c("x","y")) +
    ggtitle("Paracou - CART classification tree")+
    theme(plot.title = element_text(size=16))+
    theme(legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
  p6 = p6+stat_contour(data=plot1.melt, aes(x = (x-1)*5, y = (y-1)*5, z = z), color="black")
  p6 = p6+geom_text_contour(data=plot1.melt, aes(x = (x-1)*5, y = (y-1)*5, z=z))

  return(
    list(DiffK = DiffK,
         SC.maxtree = p1,
         SC.Classprobtree = p2,
         SC.Classtree = p3,
         C.maxtree = p4,
         C.Classprobtree= p5,
         C.Classtree = p6)
    )
}

