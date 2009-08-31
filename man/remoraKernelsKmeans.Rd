\name{remoraKernelsKmeans}
\alias{remoraKernelsKmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Remora kernels KMeans constructor }
\description{
  Remora kernels constructor for kmeans input.
}
\usage{
remoraKernelsKmeans(class_name, k_means, cluster_points)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{class\_name}{ name of the class }
  \item{k\_means}{ \code{\link[stats:kmeans]{K-Means}} result object }
  \item{cluster\_points}{ data points that belong to the clusters }
}
\details{
  This constructs a remora \code{\link[frbf:RemoraKernels-class]{kernels}} object from the \code{\link[stats:kmeans]{K-Means}} object by splitting its information and calling the \code{\link[frbf:remoraKernels]{remoraKernels}} constructor with specific arguments.
}
\value{
  Returns a remora \code{\link[frbf:RemoraKernels-class]{kernels}} object.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraKernels-class]{RemoraKernels}}
\code{\link[stats:kmeans]{K-Means}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
