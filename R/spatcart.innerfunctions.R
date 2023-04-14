#' spatcart inner function
#' @description Node initialisation.
#'
#' @param y bivariate marked point process of class 'ppp'.
#' @param offset inner parameter
#' @param parms inner parameter
#' @param wt inner parameter
#'
#' @return node initialisation.
#'
#' @import spatstat
#'
#' @export
#'
#'
spitemp = function(y, offset, parms, wt) {
  if (!is.ppp(y))
    stop("data must be of class ppp")
  if (!missing(parms) && length(parms) > 0L)
    warning("parameter argument ignored")
  if (length(offset)) y = y - offset
  sfun = function(yval, dev, wt, ylevel, digits ) {
    paste(" Label = ", format(signif(yval, digits)),
          ", Deviance = " , format(signif(dev, digits)),
          sep = '')
  }
  environment(sfun) = .GlobalEnv
  list(y = y, parms = NULL, numresp = 1L, numy = NULL, summary = sfun)
}


#' spatcart inner function
#' @description Node evaluation
#'
#' @param y bivariate marked point process of class 'ppp'.
#' @param wt inner parameter
#' @param parms inner parameter
#'
#' @return node evaluation
#'
#' @import dplyr
#' @import tidyr
#' @import spatstat
#'
#' @export
#'
#'
spetemp = function(y, wt, parms){
  # y must be of class ppp
  marque = y$marks %>% levels %>% as.factor
  yval = y %>% intensity %>% which.max
  rate = ((y$marks!=marque[yval]) %>% sum)/y$n
  nmarque = marque %>% length %>% integer
  for(i in 1:length(marque)){
    nmarque[i] = (y$marks==marque[i]) %>% sum
  }
  rss = (1-(((nmarque/y$n)^2) %>% sum))*y$n

  if (is.na(rss)) rss=0 ## if no mark

  list(label = marque[yval] , deviance = rss, rate = rate, nmarks=nmarque)
}

#' spatcart inner function
#' @description Impurity loss function on ordered variable x, at scales between 0 and r, starting form parent nodes's impurity Kt.
#'
#' @param y spatial bivariate marked point process of class 'ppp'.
#' @param r scale resolution parameter.
#' @param numvar index of variable used.
#' @param wt inner parameter.
#' @param x ordered variable numvar.
#' @param Kt parent node impurity.
#' @param parms inner parameters.
#'
#' @return a data.frame object giving impurity loss at each scale value for each split along x ordered values.
#'
#' @import dplyr
#' @import tidyr
#' @import spatstat
#'
#' @export
#'
#'
impurity = function(y, r, numvar, wt, x, Kt, parms){
  # y doit etre de classe ppp
  n = y$n
  areaNode=y$window %>% area
  # continuous x variable: do all the Ripley's K estimations
  marque = y$marks %>% levels %>% as.factor
  # i0 = y %>% intensity %>% which.max
  # largeur = min(range(y$x), range(y$y))
  # if(i0==1L) i1=2L else i1=1L
  # Kt = y %>% Kcross(marque[i0], marque[i1], seq(0L, r, by=r/500L), correction="none")
  l1 = intensity(y)[1L]
  l2 = intensity(y)[2L]

  DeltaI = data.frame(matrix(nrow=length(Kt$r), ncol=n-1L))
  temp <- rep(0L, n)
  for (i in 1L:(n-1L)) {
    temp[i] = 1L
    if(x[i+1L] - x[i]>.Machine$double.eps) {
      # left and right Ripley's K computation
      yleft=y[1L:i,]
      if(numvar==1){
        left=ppp(x=yleft$x,y=yleft$y,marks=yleft$marks,c(y$window$xrange[1],(x[i+1L]+x[i])/2),y$window$yrange)
      } else {
        left=ppp(x=yleft$x,y=yleft$y,marks=yleft$marks,y$window$xrange,c(y$window$yrange[1],(x[i+1L]+x[i])/2))
      }
      areaLeft = left$window %>% area
      if(!all(left$marks==marque[1L]) && !all(left$marks==marque[2L])){
        marque = left$marks %>% levels %>% as.factor
        i0 = left %>% intensity %>% which.max
        if(i0==1L) i1=2L else i1=1L
        Kleft = left %>% Kcross(marque[i0], marque[i1], r=Kt$r, correction="none")
      } else {
        Kleft=0L
      }

      yright=y[(i+1L):n,]
      if(numvar==1){
        right = ppp(x=yright$x,y=yright$y,marks=yright$marks,c((x[i+1L]+x[i])/2,y$window$xrange[2]),y$window$yrange)
      } else {
        right = ppp(x=yright$x,y=yright$y,marks=yright$marks,y$window$xrange,c((x[i+1L]+x[i])/2,y$window$yrange[2]))
      }
      areaRight = right$window %>% area
      if(!all(right$marks==marque[1L]) && !all(right$marks==marque[2L])){
        marque = right$marks %>% levels %>% as.factor
        i0 = right %>% intensity %>% which.max
        if(i0==1L) i1=2L else i1=1L
        Kright = right %>% Kcross(marque[i0], marque[i1], r=Kt$r, correction="none")
      } else {
        Kright=0L
        # print(i)
      }

      # Intensity and area repartition
      alphaLeft = areaLeft/areaNode
      alphaRight = areaRight/areaNode

      l1left = intensity(left)[1L]
      l2left = intensity(left)[2L]
      l1right = intensity(right)[1L]
      l2right = intensity(right)[2L]

      # Impurity loss computation
      if(!is.fv(Kleft) && !is.fv(Kright)) {
        DeltaI[,i]= Kt$un
      }
      if(!is.fv(Kleft) && is.fv(Kright)) {
        # DeltaI[,i] = Kt$un-(1L-alpha)*l1right*l2right/(l1*l2)*Kright$un
        DeltaI[,i] = Kt$un-alphaRight*l1right*l2right/(l1*l2)*Kright$un
      }
      if(!is.fv(Kright) && is.fv(Kleft)) {
        # DeltaI[,i] = Kt$un-alpha*l1left*l2left/(l1*l2)*Kleft$un
        DeltaI[,i] = Kt$un-alphaLeft*l1left*l2left/(l1*l2)*Kleft$un
      }
      if(is.fv(Kright) && is.fv(Kleft)) {
        # DeltaI[,i] = Kt$un-alpha*l1left*l2left/(l1*l2)*Kleft$un-(1L-alpha)*l1right*l2right/(l1*l2)*Kright$un
        DeltaI[,i] = Kt$un-alphaLeft*l1left*l2left/(l1*l2)*Kleft$un-alphaRight*l1right*l2right/(l1*l2)*Kright$un
      }
    } else {
      DeltaI[,i] = 0L
    }
  }
  return(DeltaI)
}


#' spatcart inner function
#' @description Goodness of splitting rules on ordered variable x
#'
#' @param y spatial bivariate marked point process of class 'ppp'.
#' @param r scale resolution
#' @param numvar index of splitting variable
#' @param wt inner parameter.
#' @param x orderd values of splitting variable numvar
#' @param Kt parent node impurity.
#' @param parms inner parameters.
#'
#' @return a list with the goodness and imputity losses at each split along ordered variable x.
#'
#' @import dplyr
#' @import tidyr
#' @import spatstat
#'
#' @export
#'
#'
spstemp = function(y, r, numvar, wt, x, Kt, parms)
{
  # y must be of class ppp
  DeltaI = y %>% impurity(r, numvar, wt, x, Kt, parms)
  Krange = which(Kt$r >= r)[1L]
  if(is.na(Krange)) {
      Krange = (Kt$un-Kt$theo) %>% which.min
  }

  goodness = DeltaI %>% slice(Krange)
  list(goodness = goodness, DeltaI=DeltaI)
}


#' spatcart inner function
#'
#' @description Maximal tree computation. Produces list of "frame" and "where" to put in class "tree".
#'
#' @param ynode spatial bivariaye marked point process restricted to node, of class 'ppp'.
#' @param r scale resolution at node
#' @param ties logical. If TRUE, in case of ties in splitting rule, the first split achieving maximum loss in impurity is taken. If FALSE, scale adaptation inside nodes is made to avoid ties.
#' @param offset inner parameter.
#' @param wt inner parameter.
#' @param parms inner parameters.
#' @param sub subset indexes of observations in node.
#' @param nodeindex node index.
#' @param minsplit minimal size of each resulting leaf for a node to be split.
#' @param minleaf minimal size of each resulting leaf for a node to be split.
#' @param mindev minimal value of node impurety to control node splitting.
#'
#' @return a list of objects : the maximal tree of class 'tree', and 2 vectors giving inside nodes resolution scales and impurity.
#'
#' @import dplyr
#' @import tidyr
#' @import spatstat
#'
#' @export
#'
#'
spfit = function(ynode, r, ties, offset, wt, parms, sub, nodeindex, minsplit, minleaf, mindev)
{
  if (!is.ppp(ynode))
    stop("data must be of class ppp")
  if(ynode$n!=sum(sub))
    stop("Response and data in node must have the same length (maximal tree)")
  if(nodeindex>=.Machine$integer.max)
    stop("Maximum depth reached in maximal tree")
  # cat(nodeindex, " n =" , sum(sub))
  # Current node initialisation
  indices = sub %>% which
  initnode = ynode %>% spitemp(offset, parms, wt)
  evalnode = ynode %>% spetemp(wt, parms)

  # Compute inner node K-function
  major =  ynode %>% intensity %>% which.max %>% names
  minor = ynode %>% intensity %>% which.min %>% names
  rvalues = seq(0L, r, by=r/500L)
  if(sum(ynode$marks==major)!=0 && sum(ynode$marks==minor)!=0){
    Kt = ynode %>% Kcross(major, minor, rvalues, correction="none")
  } else {
    Kt = data.frame(matrix(nrow=length(rvalues), ncol=2))
    colnames(Kt) = c("r", "un")
    Kt$r = rvalues
    Kt$un = 0L
  }
  # cat(" r = ", r, " K(rt) = ", Kt$un[Kt$r>=r][1L])

  # Split current node if conditions are fulfilled
  if(evalnode$deviance>mindev && sum(sub)>minsplit){

    # Current node attributes
    tree.frame$yval = evalnode$label
    tree.frame$dev = evalnode$deviance
    tree.frame$n = ynode$n
    tree.frame$yprob = matrix(evalnode$nmarks/ynode$n, ncol=length(levels(ynode$marks)))
    colnames(tree.frame$yprob) = levels(ynode$marks)
    rownames(tree.frame) = nodeindex

    # Current node splitting
    var = data.frame(cbind(ynode$x,ynode$y))
    colnames(var) = c("x","y")
    best.split = data.frame(matrix(nrow=2L,ncol=ncol(var)))

    Delta = data.frame(matrix(nrow=length(Kt$r),ncol=3L))
    colnames(Delta) = c("r",colnames(var))
    Delta[,1L] = Kt$r
    rownames(best.split) = c("threshold","Impurety.loss")
    colnames(best.split) = colnames(var)
    for(i in 1L:ncol(var)){
      # print(yorder)
      index = order(var[,i])
      x = sort(var[,i])
      yorder = ynode[index,]
      split = yorder %>% spstemp(r, i, wt, x, Kt, parms)
      ngoodness = split$goodness %>% length
      if(ngoodness>2*minleaf){
        goodness.max = max(split$goodness[(minleaf+1L):(ngoodness-minleaf-1L)])
        if(ties){
          # cat(" nb concurrents = ", length(which(split$goodness==goodness.max)))
          num = which(split$goodness==goodness.max)[1L]
        } else {
          # Computation of optimal scale r and optimal split s
          tie.index = (split$goodness==goodness.max) %>% which
          if(length(tie.index)==1){
            num = tie.index
          } else {
            Delta.tie = split$DeltaI[,tie.index]
            j = which(Kt$r>=r)[1L]
            tie.index.r = tie.index
            while(length(tie.index.r)>1 & j<nrow(Delta.tie)){
              j = j+1
              max.deltar = Delta.tie %>% slice(j) %>% max(na.rm = TRUE)
              tie.index.r = tie.index.r[Delta.tie[j,] == max.deltar]
              Delta.tie = Delta.tie[,Delta.tie[j,] == max.deltar]
            }
            num = tie.index.r[1L]
          }
        }
        if(x[num]<max(x) && x[num]>min(x)){
          best.split[1L,i] = (x[num]+x[num+1L])/2L
        } else {
          best.split[1L,i] = x[num]
        }
        best.split[2L,i] = split$goodness[num]
        Delta[,i+1L] = split$DeltaI[,num]
      } else {
        best.split[1L,i] = NA
        best.split[2L,i] = 0L
        Delta[,i+1L] = 0L
      }
    }
    i0 = which.max(best.split[2L,])[1L]
    x0 = best.split[1L,i0]
    delta = best.split[2L,i0]
    # cat(" best split = ", colnames(var)[i0], "<", x0)

    if(!is.na(x0)){
      tree.frame$var = colnames(var)[i0]
      tree.frame$splits=matrix(c(paste("<",round(x0, digits = 3)), paste(">",round(x0, digits = 3))),ncol=2L)
      colnames(tree.frame$splits) = c("cutleft","cutright")

      # Define left and right nodes
      vartemp = rep(NA,length(sub))
      vartemp[sub]=var[,i0]
      subleft = sub & vartemp<x0
      subright = sub & vartemp>=x0
    } else {
      subleft = sub
      subright = sub
    }

    critere = Delta[,i0+1]
    seuil = Delta$r[which.max(critere)]
    cp = seuil
    K = Kt$un[Kt$r>=cp][1L]
    # cat(" rt = ", cp, " K(rt) = ", K)

    # If one empty child node, or one child node having less than minleaf observations,
    # define current node as a leaf
    if(sum(subleft)==sum(sub) || sum(subright)==sum(sub) || delta<=.Machine$double.eps || sum(subleft)<minleaf || sum(subright)<minleaf){
      tree.frame$var = "<leaf>"
      tree.frame$splits = matrix(c("",""), ncol=2L)
      tree.frame$yval = evalnode$label
      tree.frame$dev = evalnode$deviance
      tree.frame$n = ynode$n
      tree.frame$yprob = matrix(evalnode$nmarks/ynode$n, ncol=length(levels(ynode$marks)))
      colnames(tree.frame$yprob) = levels(ynode$marks)
      rownames(tree.frame) = nodeindex
      wh = rep(nodeindex,sum(sub))
      names(wh)=which(sub)
      names(cp) = nodeindex
      names(K) = nodeindex
      # Else, recusively split current node
    } else {
      left = ynode[var[,i0]<x0,]
      right = ynode[var[,i0]>=x0,]
      # Window's adjustment to splitting rule
      if(i0==1){
        yleft = ppp(x=left$x,y=left$y,marks=left$marks,c(ynode$window$xrange[1],x0),ynode$window$yrange)
        yright = ppp(x=right$x,y=right$y,marks=right$marks,c(x0,ynode$window$xrange[2]),ynode$window$yrange)
        rleft = min(seuil, diff(yleft$window$xrange)/2)
        rright = min(seuil, diff(yright$window$xrange)/2)
      } else {
        yleft = ppp(x=left$x,y=left$y,marks=left$marks,ynode$window$xrange,c(ynode$window$yrange[1],x0))
        yright = ppp(x=right$x,y=right$y,marks=right$marks,ynode$window$xrange,c(x0,ynode$window$yrange[2]))
        rleft = min(seuil, diff(yleft$window$yrange)/2)
        rright = min(seuil, diff(yright$window$yrange)/2)
      }
      # cat(" seuil = ", seuil, " rleft = ", rleft, " rright = ", rright)
      names(cp) = nodeindex
      names(K) = nodeindex
      # cat(" >> ")
      tleft = yleft %>% spfit(rleft, ties, offset, wt, parms, subleft, 2*nodeindex, minsplit, minleaf, mindev)
      # cat(" >> ")
      tright = yright %>% spfit(rright, ties, offset, wt, parms, subright, 2*nodeindex+1, minsplit, minleaf, mindev)

      tree.frame=rbind(tree.frame,tleft$frame,tright$frame)
      wh=c(wh,tleft$where,tright$where)
      cp = c(cp, tleft$cp, tright$cp)
      K = c(K, tleft$K, tright$K)
    }
    # If conditions are not fulfilled, define current node as a leaf
  } else {
    tree.frame$var = "<leaf>"
    tree.frame$splits = matrix(c("",""), ncol=2L)
    tree.frame$yval = evalnode$label
    tree.frame$dev = evalnode$deviance
    tree.frame$n = sum(sub)
    tree.frame$yprob = matrix(evalnode$nmarks/ynode$n, ncol=length(levels(ynode$marks)))
    colnames(tree.frame$yprob) = levels(ynode$marks)
    rownames(tree.frame) = nodeindex
    wh = rep(nodeindex,sum(sub))
    names(wh)=which(sub)
    cp = 0L
    K = Kt$un[Kt$r>=r][1L]
    # cat(" Stop ! K(rt) = ", K)
    names(cp) = nodeindex
    names(K) = nodeindex
  }
  # Output : tree to current node and index of nodes
  list(frame = tree.frame, where = wh[order(as.integer(names(wh)))], cp = cp, K = K)
}

#' spatcart inner function
#'
#' @description Computation of the expanded heatmap grid associated to the leaves of a tree and the interactions between marks of a bivariate marked point process.
#'
#' @param tree an object of class "tree". Currently provided by the functions spatcart or spattree.
#' @param ypp the spatial bivaraitte marked point process used to generate tree a, of class 'ppp' (see spatstat).
#' @param K the values of interactions between marks of the point process ypp restricted to the leaves of tree a. Interactions are computed through the Ripley's Kcross function. They can be found in output values of functions spatcart or spattree.
#' @param d a scalar parameter defining the dxd grid resolution to produce the heatmap associated to the leaves of tree a and values of K. Default value is d = 100, leading to a 10e+4 grid of points.
#'
#' @return a list with 2 objects : the grid data.frame and the mesh = d resolution
#'
#' @import dplyr
#' @importFrom stats predict
#'
#' @export
#'
#'
make.grid = function(tree, ypp, K, d=100){

  if (sum(tree$frame$var=="<leaf>")!=length(K))
    stop("scales and tree leaves do not agree")

  x = seq(ypp$window$xrange[1], ypp$window$xrange[2], length.out = d)
  y = seq(ypp$window$yrange[1], ypp$window$yrange[2], length.out = d)

  grille = expand.grid(x=x-1/d,y=y+1/d)

  if(length(K)>1){

    grille$Interaction=numeric(d^2)
    test = data.frame(grille)
    b = predict(tree, newdata = test, type="where")

    for(i in unique(b)){
      grille$Interaction[b==i]=K[names(K)==rownames(tree$frame)[i]]
    }

  } else {
    grille$Interaction = rep(K, d^2)
  }

  list(grid = grille, mesh = d)
}
