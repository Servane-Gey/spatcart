#' Partition of a tree (with ggplot2)
#'
#' @description Plot the partition of a tree involving one or two variables. Requires package ggplot2.
#'
#' @param tree an object of class "tree". Currently provided by the functions spatcart or spattree.
#' @param ordvars the ordering of the variables to be used in a 2D plot. Specify the names in a character string of length 2; the first will be used on the x axis.
#' @param ... graphical parameters.
#'
#' @return an object of class "ggplot".
#'
#' @importFrom stats model.frame
#' @import dplyr
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' library(tree)
#' tmax = tree(Species~Sepal.Length+Petal.Length, data = iris)
#' t = prune.misclass(tmax, best = 3)
#'
#' library(dplyr)
#' library(ggplot2)
#' iris %>% select(Sepal.Length, Petal.Length, Species) %>%
#' ggplot(aes(Sepal.Length, Petal.Length, color = Species))+geom_point()+
#' gg.partition.tree(t)
gg.partition.tree <- function (tree, ordvars, ...)
{
  ptXlines <- function(x, v, xrange, xcoord = NULL, ycoord = NULL,
                       tvar, i = 1L) {
    if (v[i] == "<leaf>") {
      y1 <- (xrange[1L] + xrange[3L])/2
      y2 <- (xrange[2L] + xrange[4L])/2
      return(list(xcoord = xcoord, ycoord = c(ycoord, y1,
                                              y2), i = i))
    }
    if (v[i] == tvar[1L]) {
      xcoord <- c(xcoord, x[i], xrange[2L], x[i], xrange[4L])
      xr <- xrange
      xr[3L] <- x[i]
      ll2 <- Recall(x, v, xr, xcoord, ycoord, tvar, i +
                      1L)
      xr <- xrange
      xr[1L] <- x[i]
      return(Recall(x, v, xr, ll2$xcoord, ll2$ycoord, tvar,
                    ll2$i + 1L))
    }
    else if (v[i] == tvar[2L]) {
      xcoord <- c(xcoord, xrange[1L], x[i], xrange[3L],
                  x[i])
      xr <- xrange
      xr[4L] <- x[i]
      ll2 <- Recall(x, v, xr, xcoord, ycoord, tvar, i +
                      1L)
      xr <- xrange
      xr[2L] <- x[i]
      return(Recall(x, v, xr, ll2$xcoord, ll2$ycoord, tvar,
                    ll2$i + 1L))
    }
    else stop("wrong variable numbers in tree.")
  }
  if (inherits(tree, "singlenode"))
    stop("cannot plot singlenode tree")
  if (!inherits(tree, "tree"))
    stop("not legitimate tree")
  frame <- tree$frame
  leaves <- frame$var == "<leaf>"
  var <- unique(as.character(frame$var[!leaves]))
  if (length(var) > 2L || length(var) < 1L)
    stop("tree can only have one or two predictors")
  nlevels <- sapply(attr(tree, "xlevels"), length)
  if (any(nlevels[var] > 0L))
    stop("tree can only have continuous predictors")
  x <- rep(NA, length(leaves))
  x[!leaves] <- as.double(substring(frame$splits[!leaves, "cutleft"],
                                    2L, 100L))
  m <- model.frame(tree)
  if (length(var) == 1L) {

    a <- attributes(attr(m, "terms"))
    x <- x[!leaves]

    if (!missing(ordvars)) {
      ind <- match(var, ordvars)
      if (any(is.na(ind)))
        stop("unmatched names in vars")
      a$term.labels = ordvars
    }

    if(var == a$term.labels[1L]) var2 = a$term.labels[2L]
    if(var == a$term.labels[2L]) var2 = a$term.labels[1L]
    ry = range(m[[var2]])
    ry = ry + c(-0.025, 0.025) * diff(ry)
    xx = matrix(nrow = 4L, ncol = length(x))
    xx[1L,] = xx[3L, ] = x
    xx[2L, ] = rep(ry[1L], length(x))
    xx[4L, ] = rep(ry[2L], length(x))

    if(var == a$term.labels[1L]){
    return(
      list(
        annotate(geom="segment", x=xx[1L, ], y=xx[2L, ], xend=xx[3L, ], yend=xx[4L, ])
      )
    )
    }
    if(var == a$term.labels[2L]){
      return(
        list(
          annotate(geom="segment", x=xx[2L, ], y=xx[1L, ], xend=xx[4L, ], yend=xx[3L, ])
        )
      )
    }
    invisible(list(x = x))
  } else {
    if (!missing(ordvars)) {
      ind <- match(var, ordvars)
      if (any(is.na(ind)))
        stop("unmatched names in vars")
      var <- ordvars[sort(ind)]
    }

    rx <- range(m[[var[1L]]])
    rx <- rx + c(-0.025, 0.025) * diff(rx)
    rz <- range(m[[var[2L]]])
    rz <- rz + c(-0.025, 0.025) * diff(rz)
    xrange <- c(rx, rz)[c(1, 3, 2, 4)]
    xcoord <- NULL
    ycoord <- NULL
    xy <- ptXlines(x, frame$var, xrange, xcoord, ycoord,
                   var)
    xx <- matrix(xy$xcoord, nrow = 4L)
    yy <- matrix(xy$ycoord, nrow = 2L)
    return(
      list(
        annotate(geom="segment", x=xx[1L, ], y=xx[2L, ], xend=xx[3L, ], yend=xx[4L, ])
      )
    )
  }
}
