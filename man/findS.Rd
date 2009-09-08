\name{findS}
\alias{findS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Find s }
\description{
  Finds the s value, i.e. the spread parameter, for each cluster.
}
\usage{
findS(training_matrix, kernels, model_lambda, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training\_matrix}{ training data matrix }
  \item{kernels}{ kernels }
  \item{model\_lambda}{ previously calculated \code{\link[frbf:findLambda]{lambda}} function values }
  \item{config}{ configuration }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns the s parameter per cluster.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Fernando Martins }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}}
\code{\link[frbf:findLambda]{findLambda}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
