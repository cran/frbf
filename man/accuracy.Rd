\name{accuracy}
\alias{accuracy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculates the accuracy of a matrix }
\description{
  Calculates the accuracy of a given matrix.
}
\usage{
accuracy(accuracy_matrix, classification_vector, config, return_hit_vector = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{accuracy\_matrix}{ data matrix with accurate classification }
  \item{classification\_vector}{ the classification to verify accuracy }
  \item{config}{ the configuration }
  \item{return_hit_vector}{ if true it will return accuracy ratio per cluster and total hit ration, otherwise it will return the total number of hits }
}
\details{
  Compares the obtained classification and the real classification to check the accuracy.
}
\value{
  Returns the total number of hits or accuracy ratio per cluster and total hit ration
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}} }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
