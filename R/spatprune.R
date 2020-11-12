#' Spatial Classification Trees Pruning
#'
#' @description Spatial Classification Trees Pruning using slope heuristics. Requires package "tree".
#'
#' @param t maximal tree obtained from the spattree() or tree() function.
#' @param method optional. Parameter method = c("deviance", "misclass") sets the pruning method to use. "deviance" is related to class probability trees, while "misclass" is related to classical misclassification trees. Default is method = "deviance".
#' @param graph logical. If TRUE (Default), also produce the graphs of the number of leaves of subtrees of the pruned subtrees sequence with respect to the complexity parameter used in the pruning algorithm. The triangle symbol represents the tree selected via the modified largest jump method, while the diamond symbol represents the tree selected via the largest plateau method.
#'
#' @return a list of objects containing the sequence of subtrees pruned from tree t, the largest and smallest optimal trees obtained from the slope heuristics.
#'
#' @references *Spatial Classification Trees*, by A. Bar-Hen, S. Gey and J.-M. Poggi, HAL 01837065 (2018).
#'
#' @import dplyr
#' @import tidyr
#' @import tree
#' @importFrom stats stepfun
#' @importFrom graphics lines
#'
#' @export
#'
#' @examples
#' library(tree)
#' t = tree(Species~Sepal.Length+Petal.Length, data = iris)
#' pruned = spatprune(t)
spatprune = function(t, method = "deviance", graph = TRUE){
  if(sum(t$frame$var=="<leaf>")>1){
    seq = prune.tree(t, method=method)

    nbleav = choice.tree(seq)
    t1 = prune.tree(t, method=method, best=nbleav$maxgap)
    t2 = prune.tree(t, method=method, best=nbleav$plateau)

    if(graph){
      alpha.maxgap = seq$k[seq$size==nbleav$maxgap]
      alpha.plateau = seq$k[seq$size==nbleav$plateau]

      if(method == "deviance"){
        plot(stepfun(seq$k[-1],seq$size), main=c("Class probability trees", "Nb leaves vs complexity"), xlab="Complexity", ylab = "Number of leaves")
      }
      if(method == "misclass"){
        plot(stepfun(seq$k[-1],seq$size), main=c("Classification trees", "Nb leaves vs complexity"), xlab="Complexity", ylab = "Number of leaves")
      }
      lines(c(alpha.maxgap,alpha.maxgap), c(0, nbleav$maxgap), type = "b", pch = 2, lty="dashed", col="blue")
      lines(c(-50,alpha.maxgap), c(nbleav$maxgap, nbleav$maxgap), type = "b", pch = 2, lty="dashed", col="blue")
      lines(c(alpha.plateau,alpha.plateau), c(0, nbleav$plateau), type = "b", pch = 8, lty="dotdash", col="red")
      lines(c(-50,alpha.plateau), c(nbleav$plateau, nbleav$plateau), type = "b", pch = 8, lty="dotdash", col="red")
    }
  } else {
    t1 = t
    t2 = t
  }
  list(pruned.seq = seq, opt.tree.max=t1, opt.tree.min = t2)
}
