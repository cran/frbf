\encoding{latin1}
\name{RemoraConfiguration-class}
\docType{class}
\alias{RemoraConfiguration-class}

\title{Class "RemoraConfiguration" }
\description{ Remora configuration class }
\section{Objects from the Class}{
Objects should not be created by the user, since they are for internal use only.
}
\section{Slots}{
	 \describe{
    \item{\code{number_clusters}:}{\code{"numeric"}, the number of clusters to use in the training part }
    \item{\code{class_name}:}{\code{"character"}, the name of the training data matrix class column }
    \item{\code{weighting_function}:}{\code{"character"}, the name of the weighting function to use in the classification process }
    \item{\code{clustering_algorithm}:}{\code{"character"}, the name of the weighting function to use in the classification process }
    \item{\code{scale_variance}:}{\code{"logical"}, specifies if the scale should be performed for the principal components analysis }
    \item{\code{s}:}{\code{"numeric"}, the initial s value to use to find the s value }
    \item{\code{d}:}{\code{"numeric"}, the d value to use to find the s value }
    \item{\code{epsilon}:}{\code{"numeric"}, the eta value for function, for thoese functions that require it }
    \item{\code{niter}:}{\code{"numeric"}, the maximum number of iterations to perform to find s }
    \item{\code{niter_changes}:}{\code{"numeric"}, the number of iteraction without changes that can occur }
    \item{\code{perform_sum}:}{\code{"logical"}, specifies if the sum should be applied or not }
    \item{\code{verbose}:}{\code{"logical"}, for development debug use only }
  }
}
\section{Methods}{
No methods defined with class "remoraConfiguration" in the signature.
}
\author{ Fernando Martins }
\note{ No object of this class requires to be created since this is for frbf internal use only. }
\examples{
showClass("RemoraConfiguration")
}
\keyword{classes}
