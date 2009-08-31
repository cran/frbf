\name{buildClassification}
\alias{buildClassification}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Build classification }
\description{
  Builds the classification for each point through the distance matrix. 
  It can be done using a class sum centroid distances or just by using the centroid distance.
}
\usage{
buildClassification(distance_matrix, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance\_matrix}{ matrix with the distances between each point and each kernel }
  \item{config}{ configuration }
}
\value{
  Returns a list of points and respective class to which it belongs.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}} }% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

