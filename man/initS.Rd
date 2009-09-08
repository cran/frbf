\name{initS}
\alias{initS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Initialize s values }
\description{
  Initializes the s values.
}
\usage{
initS(distance_table, kernels, accuracy_matrix, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance\_table}{ distance table with precalculated distance point-centroid values }
  \item{kernels}{ kernel clusters }
  \item{accuracy\_matrix}{ the matrix with the accurate classification }
  \item{config}{ configuration }
}
\details{
Starts with a user given s value and modifies it until the results start getting worst.
The initial value is the same for all clusters.
}
\value{
  Returns the initial s parameter.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Fernando Martins }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}}
\code{\link[frbf:findLambda]{findLambda}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
