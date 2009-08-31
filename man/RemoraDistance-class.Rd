\encoding{latin1}
\name{RemoraDistance-class}
\docType{class}
\alias{RemoraDistance-class}

\title{Class "RemoraDistance" }
\description{Remora distance class }
\section{Objects from the Class}{
Objects should not be created by the user, since they are for internal use only.
}
\section{Slots}{ 
    \describe{
    \item{\code{index}:}{\code{"numeric"}, the point index in the matrix }
    \item{\code{distance}:}{\code{"numeric"}, the distance between the point and the cluster }
    \item{\code{className}:}{\code{"character"}, the name of the classe to which the point belongs to }
    \item{\code{point}:}{\code{"list"}, the point data }
  }
}
\section{Methods}{
No methods defined with class "RemoraDistance" in the signature.
}
\author{ Fernando Martins }
\note{ No object of this class requires to be created since this is for frbf internal use only. }
\examples{
showClass("RemoraDistance")
}
\keyword{classes}
