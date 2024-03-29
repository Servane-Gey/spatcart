% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatprune.R
\name{spatprune}
\alias{spatprune}
\title{Spatial Classification Trees Pruning}
\usage{
spatprune(t, method = "deviance", graph = TRUE)
}
\arguments{
\item{t}{maximal tree obtained from the spattree() or tree() function.}

\item{method}{optional. Parameter method = c("deviance", "misclass") sets the pruning method to use. "deviance" is related to class probability trees, while "misclass" is related to classical misclassification trees. Default is method = "deviance".}

\item{graph}{logical. If TRUE (Default), also produce the graphics of the number of leaves of subtrees of the pruned subtrees sequence with respect to the complexity parameter used in the pruning algorithm. The triangle symbol represents the tree selected via the modified largest jump method, while the diamond symbol represents the tree selected via the largest plateau method.}
}
\value{
a list of objects containing the sequence of subtrees pruned from tree t, the largest and smallest optimal trees obtained from the slope heuristics.
}
\description{
Spatial Classification Trees Pruning using slope heuristics. Requires package "tree".
}
\examples{
library(tree)
t = tree(Species~Sepal.Length+Petal.Length, data = iris)
pruned = spatprune(t)
}
\references{
\emph{Spatial Classification Trees}, by A. Bar-Hen, S. Gey and J.-M. Poggi (2021).
}
