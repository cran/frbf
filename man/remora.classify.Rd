\name{remora.classify}
\alias{remora.classify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Remora classify }
\description{
  Predicts the classification.
}
\usage{
remora.classify(model, data_matrix, return_nice = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ remora \code{\link[frbf:RemoraModel-class]{model}} }
  \item{data\_matrix}{ data matrix to classify }
  \item{return\_nice}{ specifies if the result should be user pretty printable or just raw, default is False which means raw }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns classification prediction.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Fernando Martins }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraModel-class]{RemoraModel}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
