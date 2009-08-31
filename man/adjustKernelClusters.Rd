\name{adjustKernelClusters}
\alias{adjustKernelClusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Adjust kernel clusters }
\description{
Adjustes the number of kernel clusters to fit the number of clusters  specified by the user. 
When that is not possible, the number of kernel clusters will increase to force at least one cluster per kernel.
This function issues warnings when adjustments have been  made to the kernel clusters.
}
\usage{
adjustKernelClusters(kernels, clusters, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kernels}{ kernel clusters }
  \item{clusters}{ number of desired clusters }
  \item{config}{ configuration }
}
\details{
  Warnings will be issued if the kernels are adjusted.
}
\value{
  Returns the clusters adjusted.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}} }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

