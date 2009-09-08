\name{train}
\alias{train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Trainging step }
\description{
  Trains an FRBF network. Performs the learning step of the algorithm.
}
\usage{
train(training_matrix, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{training\_matrix}{ training data matrix }
  \item{config}{ configuration }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns the kernels from the training step.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Fernando Martins }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
