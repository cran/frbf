\encoding{latin1}
\name{RemoraModel-class}
\docType{class}
\alias{RemoraModel-class}

\title{Class "RemoraModel" }
\description{ Remora model class}
\section{Objects from the Class}{
Objects should not be created by the user, since they are for internal use only.}
\section{Slots}{
	 \describe{
    \item{\code{config}:}{\code{"RemoraConfiguration"}, the configuration used to build the model }
    \item{\code{model}:}{\code{"data.frame"}, the matrix with the model data information }
    \item{\code{lambda}:}{\code{"list"}, the function lambda }
    \item{\code{kernels}:}{\code{"list"}, is the list of the kernels }
  }
}
\section{Methods}{
No methods defined with class "remoraConfiguration" in the signature.
}
\author{ Fernando Martins }
\note{ No object of this class requires to be created since this is the result of the frbf function call. }
\seealso{
  \code{\linkS4class{RemoraConfiguration}}
  \code{\linkS4class{RemoraModel}}
  \code{\linkS4class{RemoraKernels}}
  \code{\link[frbf:frbf]{frbf}}
}
\examples{
showClass("RemoraModel")
}
\keyword{classes}
