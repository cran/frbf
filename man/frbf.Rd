\encoding{latin1}
\name{frbf}
\alias{frbf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Flexible Radial Basis Function }
\description{  
Models kernels in Radial Basis Function networks. 
The function provides flexibility to the kernel structure. 
This flexibility comes through the use of modifier functions applied to the distance computation procedure, 
essential for all kernel evaluations.
}
\usage{
frbf(data_matrix, number_clusters, class_name, 
   weighting_function = "euclidean", 
   scale_variance = TRUE, 
   s_value = 0.2, 
   d = 0.23,
   epsilon = 0.01, 
   niter = -1, 
   niter_changes = 5, 
   perform_sum = TRUE, 
   clustering_algorithm = "",
   verbose="no")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_matrix}{  the data to use to train the algorithm, must be a matrix or a data frame }
  \item{number_clusters}{  the total number of clusters, may be adjusted during the execution and will be used by the \code{\link[stats:kmeans]{K-Means}} algorithm }
  \item{class_name}{ the name, or index, of the column that holds the class of the training data matrix }
  \item{weighting_function}{  the name of the kernel function from the options "euclidean", "one\_minus", "one\_minus\_sq", "mahalanobis", 
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
  
The \code{epsilon} will only used for the functions that require it.

The \code{scale_variance} is the \code{\link[stats:prcomp]{prcomp}} \code{scale} parameter.

%The kernel functions available for the \code{weighting\_function} parameter are:
% - FUNCTION\_REMORA\_EUCLIDEAN value "euclidean"
% - FUNCTION\_REMORA\_ONE\_MINUS value "one\_minus"
% - FUNCTION\_REMORA\_ONE\_MINUS\_SQ value "one\_minus\_sq"
% - FUNCTION\_REMORA\_MAHALANOBIS value "mahalanobis"
% - FUNCTION\_REMORA\_EXP\_ONE\_MINUS value "exp\_one\_minus"
% - FUNCTION\_REMORA\_EXP\_ONE\_MINUS\_SQ value "exp\_one\_minus\_sq"
% - FUNCTION\_REMORA\_EXP\_ONE\_LOG value "exp\_one\_log"
% - FUNCTION\_REMORA\_NORMALIZED\_DIFFERENCE value "normalized\_difference"
% - FUNCTION\_REMORA\_NORMALIZED\_DIFFERENCE\_SQ value "normalized\_difference\_sq"
%
%The \code{verbose} parameter may be set to the following values:
% - VERBOSE\_NO value "no", performs silently
% - VERBOSE\_YES value "yes", display some information about the execution steps
% - VERBOSE\_DETAIL value "detail", display detailed information about the execution steps
% - VERBOSE\_DEBUG value "debug", displays debug information
The verbosity set in the \code{frbf} will affect the \code{\link[frbf:predict.RemoraModel]{predict}} function verbosity as well.
}
\value{
  Returns a \code{RemoraModel} S4 object with the following slots:
  \item{config }{the configuration used on the model}
  \item{model }{the matrix with the model data information }
  \item{lambda }{the function lambda calculated per each cluster }
  \item{kernels }{the model kernels, the result from the kmeans }
}
\references{ 
Andre O. Falcao, Thibault Langlois and Andreas Wichert
(2006) \emph{Flexible kernels for RBF networks}. 
Jornal of Neurocomputing, volume 69, pp 2356-2359. Elsevier.
 }
\author{ Fernando Martins and Andre Falcao }
\note{ 
%\section{Warning }
Some user defined values may be overwriten by the function if they are considered unsuitable for the learning procedure.
In such cases, warnings will be issued describing what has been changed.
}
\seealso{ 
\code{\link[stats:kmeans]{kmeans}}
\code{\link[stats:prcomp]{prcomp}}
\code{\link[frbf:RemoraModel-class]{RemoraModel}}
}
\examples{
# Example 
# data will be splitted into training and classification groups
data(iris) 
samp <- sample(1:150, 75)
# the training matrix will be use 75 random points
training_matrix <- iris[samp, ]
# the matrix to classify will use all the other points 
classification_matrix <- iris[-samp, ]

# create the model
model <- frbf(training_matrix, weighting_function="mahalanobis", class_name = "Species", number_clusters = 3, scale_variance = TRUE)

# the configuration used on the model 
print(model@config) 
# the matrix with the model data information 
print(model@model) 
# the function lambda calculated per each cluster 
print(model@lambda) 
# the model kernels, the result from the kmeans 
print(model@kernels) 

# predict 
classification <- predict(model, classification_matrix) 
# the classification points for the last 
print(classification)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ cluster }
\keyword{ neural }
% __ONLY ONE__ keyword per line
