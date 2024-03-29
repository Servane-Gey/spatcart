% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Paracou.plot.R
\name{Paracou.plot}
\alias{Paracou.plot}
\title{Heatmap of interaction between species of the Paracou data set}
\usage{
Paracou.plot(a, ypp, K, d = 100)
}
\arguments{
\item{a}{an object of class "tree". Currently provided by the functions spatcart or spattree.}

\item{ypp}{the Paracou bivariate marked spatial point process used to generate tree a.}

\item{K}{the values of interactions between marks of the point process ypp restricted to the leaves of tree a. Interactions are computed through the Ripley's Kcross function. They can be found in output values of functions spatcart or spattree.}

\item{d}{a scalar parameter defining the dxd grid resolution to produce the heatmap associated to the leaves of tree a and values of K. Default value is d = 100, leading to a 10e+4 grid of points.}
}
\value{
an object of class "ggplot".
}
\description{
Plot the heatmap of interaction between marks with altitudes information for the Paracou data set.
}
\examples{
library(spatcart)
# Loading Paracou bivariate spatial point process
ypp = Paracoudata()

# SpatCART trees
t = spatcart(ypp, 15, graph = FALSE)

a = t$opt.tree.min
K = t$K[rownames(a$frame)][a$frame$var=="<leaf>"]

Paracou.plot(a,ypp,K)
}
\references{
\emph{Spatial Classification Trees}, by A. Bar-Hen, S. Gey and J.-M. Poggi (2021).
}
