\encoding{latin1}
\name{predict.RemoraModel}
\alias{predict.RemoraModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Predict Classification }
\description{
The predict.RemoraModel funcion overloads the predict function to support the trained
model contained in the \code{\link[frbf:RemoraModel-class]{RemoraModel}} class.
This function recieves the model, a data matrix to classify, and classifies, i.e. predicts, to which class each of the given data points belong to.
}
\usage{
predict.RemoraModel(object, data_matrix, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ the \code{\link[frbf:RemoraModel-class]{model}}, obtained from the learning procedure (see \code{\link[frbf:frbf]{frbf}}) } 
  \item{data_matrix}{ the data to classify }
  \item{...}{ additional arguments affecting the predictions produced }
}
\details{
  The \code{data_matrix} can be a \code{matrix} or \code{data.frame}. 
  It can have the class column if it has the same name, or index, as the training matrix used in \code{\link[frbf:frbf]{frbf}}). 
  In such case it will be automatically ignored, otherwise the class column cannot be present in the data set.
}
\value{
  The result is a prediction list containing the class name of each data point. The position of the data point in result is the same as the position in the matrix given for classification.
}
\references{ Andr\'{e} O. Falc{\~{a}}o, Thibault Langlois and Andreas Wichert
(2006) \emph{Flexible kernels for RBF networks}. 
Jornal of Neurocomputing, volume 69, pp 2356-2359. Elsevier. }
\author{ Fernando Martins and Andr\'{e} Falc{\~{a}}o }
\note{ }
\seealso{ \code{\link[frbf:frbf]{frbf}} 
\code{\link[frbf:RemoraModel-class]{RemoraModel}}
}
\examples{
# infert data is composed by 248 points and will be splitted 
data(infert) 
# the training matrix will be use the first 100 points 
training_matrix <- infert[c(1:100) ,] 
# the matrix to classify will use all the other points 
classification_matrix <- infert[c(101:248) ,] 

# create the model
model <- frbf(training_matrix, class_name = "education", number_clusters = 10, scale_variance = FALSE)

# predict 
classification <- predict(model, classification_matrix) 

# the classification points for the last 
print(classification) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classif }% __ONLY ONE__ keyword per line
