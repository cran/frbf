\encoding{latin1}
\name{frbf-package}
\alias{remora-package}
\docType{package}
\title{
Flexible Radial Basis Function Package
}
\description{
Package containing the implementation of the Flexible Radial Basis Function.
}
\details{
\tabular{ll}{
Package: \tab remora\cr
Type: \tab Package\cr
Version: \tab 0.9\cr
Date: \tab 2009-06-13\cr
License: \tab GPL 3.0\cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
}
The package provides the Flexible Radial Basis Function that will create a classification model usable on the predict function.
}
\author{
Fernando Martins
Andre Falcao

Maintainers: Fernando Martins <fmp.martins@gmail.com>, Andre Falcao <afalcao@di.fc.ul.pt>
}
\references{
Andre O. Falcao, Thibault Langlois and Andreas Wichert
(2006) \emph{Flexible kernels for RBF networks}. 
Jornal of Neurocomputing, volume 69, pp 2356-2359. Elsevier.
}
\keyword{ package }
\seealso{
\code{\link[stats:kmeans]{kmeans}}
\code{\link[stats:prcomp]{prcomp}}
}
\examples{
# infert data is composed by 248 points and will be splitted 
data(infert) 
# the training matrix will be use the first 100 points 
training_matrix <- infert[c(1:100) ,] 
# the matrix to classify will use all the other points 
classification_matrix <- infert[c(101:248) ,] 

# create the model
model <- frbf(training_matrix, weighting_function="mahalanobis", class_name = "education", number_clusters = 3, scale_variance = FALSE)

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
