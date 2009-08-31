\name{distances}
\alias{distances}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Distances }
\description{
  Calculates the difference between each centroid and each point.
}
\usage{
distances(distance_values, kernels, kernels_s, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance\_values}{ distance table with precalculated distance point-centroid values  }
  \item{kernels}{ kernels }
  \item{kernels\_s}{ s value per each kernel }
  \item{config}{ configuration }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns a distance table between points and classes, matrix row number is the point row index.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
