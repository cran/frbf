\name{remoraKernels}
\alias{remoraKernels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Remora kernels constructor }
\description{
   Remora kernels constructor.
}
\usage{
remoraKernels(class_name, clusters, cluster_points, points_per_cluster, centroids, size, eigen_values, eigen_vector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{class\_name}{ name of the class }
  \item{clusters}{ vector of integers indicating the cluster to which each point is allocated }
  \item{cluster\_points}{ matrix with the point that belong to the clusters }
  \item{points\_per\_cluster}{ list with cluster points organized by clusters }
  \item{centroids}{ matrix of cluster centres }
  \item{size}{ number of points in each cluster }
  \item{eigen\_values}{ list of cluster eigen values }
  \item{eigen\_vector}{ list of cluster eigen vector matrices }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns a remora \code{\link[frbf:RemoraKernels-class]{kernels}} object.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Fernando Martins }
\note{ This function is part of the frbf implementation and is for internal use only. }
\seealso{ \code{\link[frbf:RemoraKernels-class]{RemoraKernels}}
\code{\link[stats:kmeans]{K-Means}}
\code{\link[stats:prcomp]{prcomp}} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }