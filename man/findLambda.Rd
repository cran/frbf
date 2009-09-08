\name{findLambda}
\alias{findLambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Finds the lambda function values }
\description{
 Calculates the lambda function values, i.e. the kernel function values. 
}
\usage{
findLambda(training_matrix, kernels, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training\_matrix}{ training data matrix }
  \item{kernels}{ kernels }
  \item{config}{ configuration }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns the calculated lambda function values.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Fernando Martins }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

