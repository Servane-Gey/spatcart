#' Optimal subtree choice
#'
#' @description Optimal subtree choice using slope heuristics proposed in "Spatial Classification Trees" by Bar-Hen, Gey, Poggi.
#'
#' @param seq sequence of subtrees pruned from a maximal tree. Currently provided by functions prune.tree from package tree, spatprune and spatcart.
#' @param coeff a scalar parameter adjusting the size of the modified maximal gap slope heuristic. Default is to take 0.1x(size of the largest plateau).
#'
#' @return a tibble object with 2 arguments
#'
#' @references *Spatial Classification Trees*, by A. Bar-Hen, S. Gey and J.-M. Poggi (2021).
#'
#' @import dplyr
#' @import tidyr
#'
#' @export
#'
#' @examples
#' library(tree)
#' t = tree(Species~., data = iris)
#' seq = prune.misclass(t)
#' choice.tree(seq)
#'
choice.tree = function(seq, coeff = 0.1){

  if(length(seq$size[-1])==1){
    best.plat = seq$size[1]
    best.saut = seq$size[1]
  } else {
    comp = c(0,seq$k[-1])
    diff.k = diff(comp)
    diff.size = -diff(seq$size)
    ind.k = which.max(diff.k)
    ind.size = which.max(diff.size)

    # Selection with Plateau method
    best.plat = seq$size[ind.k]

    # Selection with regularized largest jump
    plateau = max(diff.k)
    jump = max(diff.size)
    saut = max(which(diff.size==jump))
    dff= diff.k<=(plateau*coeff)

    if(which.max(dff)>=ind.k || ind.size>=ind.k){
      best.saut = best.plat
    } else {
      plage = -((ind.k+1):length(seq$size))
      sauts.som = -diff(seq$size[plage])
      seq.som = seq$size[plage]
      diff.k = diff(seq$k[plage][-1])
      dff = diff.k<=(plateau*coeff)
      if (sum(dff)>=1){
        temp = sauts.som
        for (j in 1:length(dff)){
          if(dff[j]){
            temp[j+1]= temp[j]+temp[j+1]
            temp[j]=NA
          }
        }
        sauts.som = temp[!is.na(temp)]
        seq.som = c(seq.som[1],seq.som[2:length(seq.som)][!is.na(temp)])
      }
      saut.grand = max(which(sauts.som==max(sauts.som)))

      if(jump>max(sauts.som)){
        best.saut = seq$size[saut+1]
      } else {
        best.saut=seq.som[saut.grand+1]
      }
    }
  }
  sortie = tibble(
    plateau = best.plat,
    maxgap = best.saut
    )
  return(sortie)
}

