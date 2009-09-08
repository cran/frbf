\name{buildDistanceTable}
\alias{buildDistanceTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Build distance table }
\description{
   Build the distance table between each centroid and each point.
}
\usage{
buildDistanceTable(data_matrix, kernels, model_lambda, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data\_matrix}{ data matrix to classify }
  \item{kernels}{ kernels }
  \item{model\_lambda}{ previously calculated lambda function values }
  \item{config}{ configuration }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns the distance table between points and classes, matrix row number is the point row index.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Fernando Martins }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }


