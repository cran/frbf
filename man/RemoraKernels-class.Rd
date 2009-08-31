\encoding{latin1}
\name{RemoraKernels-class}
\docType{class}
\alias{RemoraKernels-class}

\title{Class "RemoraKernels" }
\description{ Remora kernels class }
\section{Objects from the Class}{
Objects should not be created by the user, since they are for internal use only.
}
\section{Slots}{
	 \describe{
    \item{\code{class_name}:}{\code{"character"}, the name of the class }
    \item{\code{eigen_values}:}{class \code{"list"}, the list of cluster eigen values }
    \item{\code{eigen_vector}:}{\code{"list"}, the list of cluster eigen vector matrices }
    \item{\code{clusters}:}{\code{"numeric"}, vector of integers indicating the cluster to which each point is allocated }
    \item{\code{cluster_points}:}{\code{"data.frame"}, matrix with the point that belong to the clusters }
    \item{\code{points_per_cluster}:}{\code{"list"}, list with cluster points organized by clusters }
    \item{\code{centroids}:}{\code{"matrix"}, matrix of cluster centres }
    \item{\code{size}:}{\code{"numeric"}, number of points in each cluster }
  }
}
\section{Methods}{
No methods defined with class "remoraKernels" in the signature.
}
\author{ Fernando Martins }
\note{ No object of this class requires to be created since this is for frbf internal use only. }
\examples{
showClass("RemoraKernels")
}
\keyword{classes}
