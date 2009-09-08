\name{remoraDistance}
\alias{remoraDistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Remora distance constructor }
\description{
  Builds a distance object.
}
\usage{
remoraDistance(index, distance = NULL, point = NULL, class_name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{index}{ point index in the matrix }
  \item{distance}{ distance between the point and the cluster }
  \item{point}{ data point }
  \item{class\_name}{ name of the classe to which the point belongs to }
}
\details{
  Null values will result in no slot value creation.
}
\value{
  Returns the result distance in a \code{\link[frbf:RemoraDistance-class]{RemoraDistance}} object
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Fernando Martins }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraKernels-class]{RemoraKernels}}
\code{\link[stats:kmeans]{K-Means}}
\code{\link[stats:prcomp]{prcomp}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }