\name{remora.model}
\alias{remora.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Remora Model}
\description{
Builds the Remora model.
Recieves the training matrix and a configuration and returns a model.
}
\usage{
remora.model(training_matrix, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training\_matrix}{ data matrix for training }
  \item{config}{ configuration }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns remora model.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Fernando Martins }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}}
\code{\link[frbf:RemoraModel-class]{RemoraModel}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
