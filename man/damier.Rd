% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/damier.R
\name{damier}
\alias{damier}
\title{Chess simulated data set}
\usage{
damier(n, model = "unif", h = 0.4, d = 3, graph = TRUE)
}
\arguments{
\item{n}{if model = "unif", number of observations. If model = "Poisson", intensity of the Poisson point process used to simulate data (see spatstat).}

\item{model}{nominal, model = c("unif", "Poisson"). If model = "unif", generates points uniformly on the window. If model = "Poisson", generates points as a Poisson point process in the window. Default is model = "unif".}

\item{h}{margin parameter defining the mixing proportion of marks inside windows. Default is h = 0.4.}

\item{d}{set the number of squares in the chessboard as dxd. Default is d = 3.}

\item{graph}{logical. If TRUE (default), plot the data scatter plot and the chessboard's squares.}
}
\value{
a list of data.frame objects with data points, marks distribution and Bayes predictions.
}
\description{
Simulation of a bivariate marked spatial point process with marks distribued as a chessboard.
}
\examples{
library(spatstat)
damier(1000, model = "unif")
damier(1000, model = "unif", h = 0.1)
damier(1000, model = "Poisson")
}
