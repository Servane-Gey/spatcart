#' Partition of a tree
#'
#' @description  Plot the partition of a tree involving one or two variables. Corrected version of the {tree} function partition.tree to be able to handle trees partitioning on only one variable.
#'
#' @param tree an object of class "tree". Currently provided by the functions spatcart or spattree.
#' @param ordvars the ordering of the variables to be used in a 2D plot. Specify the names in a character string of length 2; the first will be used on the x axis.
#' @param label a character string giving the column of the frame component of tree to be used to label the regions.
#' @param add if true, add to existing plot, otherwise start a new plot.
#' @param ordvars the ordering of the variables to be used in a 2D plot. Specify the names in a character string of length 2; the first will be used on the x axis.
#' @param ... graphical parameters.
#'
#' @return None
#'
#' @import dplyr
#' @importFrom stats model.frame
#' @importFrom graphics segments
#' @importFrom graphics text
#' @importFrom stats predict
#'
#' @export
#'
partition.tree.corrected = function (tree, label = "yval", add = FALSE, ordvars, ...)
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

    # Segments
    if(var == a$term.labels[1]) var2 = a$term.labels[2L]
    if(var == a$term.labels[2]) var2 = a$term.labels[1L]
    ry = range(m[[var2]])
    ry = ry + c(-0.025, 0.025) * diff(ry)
    xx = matrix(nrow = 4L, ncol = length(x))
    xx[1L,] = xx[3L, ] = x
    xx[2L, ] = rep(ry[1L], length(x))
    xx[4L, ] = rep(ry[2L], length(x))

    # Labels
    resp = frame$yval[leaves]

    rx = range(m[[var]])
    rx = rx + c(-0.025, 0.025) * diff(rx)
    limits = sort(c(rx,x))
    yy = matrix(nrow = 2L, ncol = length(resp))
    yy[1L, ] = sort(c(rx[1L],x)) + diff(limits)/2
    yy[2L, ] = rep(ry[2L]/2, length(resp))

    if (!add){
      rz = range(m[[var2]])
      rz = rz + c(-0.025, 0.025) * diff(rz)

      plot(rx, rz, xlab = a$term.labels[1L], ylab = a$term.labels[2L], type = "n",
           xaxs = "i", yaxs = "i", ...)
    }

    if(var == a$term.labels[1L]){
      # Segments
      segments(xx[1L, ], xx[2L, ], xx[3L, ], xx[4L, ])

      # Labels
      b = data.frame(
        x1 = yy[1L, ],
        y1 = yy[2L, ]
      )
      colnames(b) = c(a$term.labels[1L], a$term.labels[2L])
      nodes = predict(tree, newdata = b, type = "where")
      lab = factor(rep(NA, length(resp)), levels = levels(m[[1L]]))

      for(k in unique(nodes)){
      lab[nodes == k] = resp[which(leaves) == k]
      }

      text(yy[1L, ], yy[2L, ], as.character(lab), ...)
    }
    if(var == a$term.labels[2L]){
      # Segments
      segments(xx[2L, ], xx[1L, ], xx[4L, ], xx[3L, ])

      # Labels
      b = data.frame(
        x1 = yy[2L, ],
        y1 = yy[1L, ]
      )
      colnames(b) = c(a$term.labels[1L], a$term.labels[2L])
      nodes = predict(tree, newdata = b, type = "where")
      lab = factor(rep(NA, length(resp)), levels = levels(m[[1L]]))

      for(k in unique(nodes)){
        lab[nodes == k] = resp[which(leaves) == k]
      }

      text(yy[2L, ], yy[1L, ], as.character(lab), ...)
    }

    invisible(list(x = x))
  }
  else {
    if (!missing(ordvars)) {
      ind <- match(var, ordvars)
      if (any(is.na(ind)))
        stop("unmatched names in vars")
      var <- ordvars[sort(ind)]
    }
    lab <- frame$yval[leaves]
    if (is.null(frame$yprob))
      lab <- format(signif(lab, 3L))
    else if (match(label, attr(tree, "ylevels"), nomatch = 0L))
      lab <- format(signif(frame$yprob[leaves, label],
                           3L))
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
    if (!add)
      plot(rx, rz, xlab = var[1L], ylab = var[2L], type = "n",
           xaxs = "i", yaxs = "i", ...)
    segments(xx[1L, ], xx[2L, ], xx[3L, ], xx[4L, ])
    text(yy[1L, ], yy[2L, ], as.character(lab), ...)
  }
}
