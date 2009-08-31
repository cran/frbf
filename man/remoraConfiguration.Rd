\name{remoraConfiguration}
\alias{remoraConfiguration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Remora configuration constructor }
\description{
  Builds a remora configuration object.
}
\usage{
remoraConfiguration(number_clusters, class_name, weighting_function = FUNCTION_REMORA_EUCLIDEAN, scale_variance = TRUE, s_value = 0.2, d = 0.23, epsilon = 0.01, niter = -1, niter_changes = 5, perform_sum = TRUE, clustering_algorithm = "", verbose = VERBOSE_NO)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{number_clusters}{  the total number of clusters, may be adjusted during the execution and will be used by the \code{\link[stats:kmeans]{K-Means}} algorithm }
  \item{class_name}{ the name, or index, of the column that holds the class of the training data matrix }
  \item{weighting_function}{ the name of the kernel function from the options "euclidean", "one\_minus", "one\_minus\_sq", "mahalanobis", 
        "exp\_one\_minus", "exp\_one\_minus\_sq", "exp\_one\_log", "normalized\_difference", "normalized\_difference\_sq", if none is specified, the "euclidean" function will be used }
  \item{scale_variance}{ specifies if the scale should be performed for the principal components analysis, default is True (see \code{\link[stats:prcomp]{prcomp}}) }
  \item{s_value}{ the initial value to use to find the \code{s} value, to use to find the kernels sigma value, the spread parameter adjustment, it has a default value of 0.2 }
  \item{d}{ the initial \code{d} value to use to find the \code{s} value, it will use the default value of 0.23 if no value is specified }
  \item{epsilon}{ the \code{epsilon} for the functions that require it, if none is specified a default value of 0.01 will be used }
  \item{niter}{ the maximum number of iterations to perform to find \code{s}, if no value is provided, a default will be calculated based on the number of training data points}
  \item{niter_changes}{ the number of iteraction without changes that can occur, if this number is reached withou any change, the iteration will stop,  if no value is specified 5 will be used by default }
  \item{perform_sum}{ specifies if the sum of the centroids per cluster should be applied, or not, default is True }
  \item{clustering_algorithm}{ specifies which of the \code{\link[stats:kmeans]{K-Means}} algorithms should be used, if none is specified, the default \code{\link[stats:kmeans]{K-Means}} algorithm will be used }
  \item{verbose}{ specifies the algorithm verbosity during the execution, from "no", "yes", "detail", "debug" options, if nothing is specified it will be "no" }
}
\details{
  Missing values for parameters that do not have a default value will halt the \code{\link[frbf:frbf]{frbf}} execution.
}
\value{
  Returns the configuration in \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}} object.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraConfiguration-class]{RemoraConfiguration}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }