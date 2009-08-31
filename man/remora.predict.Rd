\name{remora.predict}
\alias{remora.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Remora predict }
\description{
  Predicts a given classification using a specific model.
}
\usage{
remora.predict(model, data_matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ remora \code{\link[frbf:RemoraModel-class]{model}} }
  \item{data\_matrix}{ data matrix to classify }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns the list containing Remora classification prediction
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraModel-class]{RemoraModel}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
