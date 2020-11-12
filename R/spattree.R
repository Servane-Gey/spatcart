#' Spatial Classification Maximal Tree
#'
#' @description Spatial Classification Maximal Tree construction using the Ripley's intertype K-function to partition spatial bivariate marked point process. Requires packages "spatstat" and "tree".
#'
#' @param y spatial bivariate marked point process of class 'ppp' (see spatstat).
#' @param r initial resolution scale to use to construct spatcart maximal tree. Decreasing with splitting..
#' @param ties logical. If TRUE (default), in case of ties in splitting rule, the first split achieving maximum loss in impurity is taken. If FALSE, scale adaptation inside nodes is made to avoid ties.
#' @param offset parameter related to the tree() function of package tree.
#' @param wt vector of non-negative observational weights; fractional weights are allowed.
#' @param parms other parameters related to the tree function of package tree.
#' @param minsplit minimal size for the node to be split. Default is minsplit = 10.
#' @param minleaf minimal size of each resulting leaf for a node to be split. Default is minleaf = 5.
#' @param mindev minimal value of node impurety to control node splitting. Default is mindev = 0.
#'
#' @return a list of objects containing the maximal tree; the sequence of pruned subtrees and the corresponding complexities; the values of scale resolution and impurity inside nodes of the maximal tree.
#'
#' @import dplyr
#' @import tidyr
#' @import spatstat
#' @import tree
#'
#' @details parameter minsplit must be at least 2xminleaf.
#'
#' @references *Spatial Classification Trees*, by A. Bar-Hen, S. Gey and J.-M. Poggi, HAL 01837065 (2018).
#'
#' @export
#'
#' @examples
#' library(spatcart)
#' # Loading Paracou data set
#' ypp = Paracoudata()
#'
#' # Keeping ties in splitting rules
#' t = spattree(ypp, 15)
#'
#' # optimizing ties in splitting rules
#' t = spattree(ypp, 15, ties = FALSE)
#'
#' # Changing stopping rule for splitting
#' t = spattree(ypp, 15, minsplit = 20, minleaf = 10)
spattree = function(y, r, ties = TRUE, offset=NULL, wt=NULL, parms=NULL, minsplit=10, minleaf=5, mindev=0){
  if (!is.ppp(y))
    stop("data must be of class ppp")

  # Maximal Tree
  var=""
  n = 0L
  dev = 0L
  yval = ""
  splits = matrix(c("",""), ncol=2L)
  colnames(splits) = c("cutleft","cutright")
  yprob = matrix(ncol=length(levels(y$marks)))
  colnames(yprob) = levels(y$marks)

  # Define initialization and outpout to function fit as "Global"
  tree.frame <<- data.frame(var=var,n=n,dev=dev,yval=yval,splits=I(splits),yprob=I(yprob))
  where <<- integer()
  cp <<- numeric()
  K <<- numeric()

  sub = rep(TRUE,y$n)
  nodeindex = 1L
  spatt = y %>% spfit(r, ties, offset, wt, parms, sub, nodeindex, minsplit, minleaf, mindev)

  spatt$frame$n = spatt$frame$n %>% as.numeric
  spatt$frame$var = spatt$frame$var %>% as.factor
  for(i in 1:length(spatt$where)){
    spatt$where[i] = (rownames(spatt$frame)==spatt$where[i]) %>% which
  }
  spatt$where = spatt$where %>% as.integer
  names(spatt$where) = 1L:length(spatt$where)

  donnees = data.frame(y)
  t = tree(marks~., data=donnees, minsize=minsplit, mincut=minleaf, mindev=mindev, model=T)

  t$frame = spatt$frame
  t$where = spatt$where
  t$y = y$marks
  t$weigths = wt
  # t$model = data.frame(y)
  if(sum(t$frame$var=="<leaf>")>1) class(t)="tree"
  list(tree = t, cp=spatt$cp, K=spatt$K)
}
