\name{classificationDistanceTable}
\alias{classificationDistanceTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Classification distance table }
\description{
  Builds the classification distance table for the data matrix points.
}
\usage{
classificationDistanceTable(remora_model, data_matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{remora\_model}{ \code{\link[frbf:RemoraModel-class]{Remora Wmodel}} }
  \item{data\_matrix}{ data matrix that holds the points to classify }
}
\value{
  Returns a distance table.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Fernando Martins }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraModel-class]{RemoraModel}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
