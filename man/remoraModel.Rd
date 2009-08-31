\name{remoraModel}
\alias{remoraModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Remora model constructor }
\description{
  Builds a remora model object.
}
\usage{
remoraModel(config, model_matrix, model_lambda, model_kernels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{config}{ configuration }
  \item{model\_matrix}{ matrix with the model data information }
  \item{model\_lambda}{ lambda function values }
  \item{model\_kernels}{ list of the kernels }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns the model in \code{\link[frbf:RemoraModel-class]{RemoraModel}} object.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}}
\code{\link[frbf:RemoraModel-class]{RemoraModel}}
\code{\link[frbf:findLambda]{findLambda}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }