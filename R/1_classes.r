#
# Remora implementation.
# Global definitions.
#


#
# Class definitions.
# S4 implementation.
#

#
# Class distance definition.
#
CLASS_REMORA_DISTANCE <- "RemoraDistance"
class_distance <- setClass(CLASS_REMORA_DISTANCE, representation(index = "numeric", distance = "numeric", className = "character", point = "list"), prototype = list(index=numeric(), distance=numeric(), className=character(), point=list()))

#
# Class configuration definition.
#
CLASS_REMORA_CONFIGURATION <- "RemoraConfiguration"
class_configuration <- setClass(CLASS_REMORA_CONFIGURATION, 
  representation(number_clusters = "numeric", class_name = "character", weighting_function="character", clustering_algorithm="character",
    scale_variance="logical", s="numeric", d="numeric", epsilon = "numeric", niter="numeric", niter_changes="numeric", perform_sum="logical", verbose="character"), 
  prototype = list(number_clusters=numeric(), class_name=character(), weighting_function=character(), clustering_algorithm=character(),
    scale_variance=logical(), s=numeric(), d=numeric(), epsilon=numeric(), niter=numeric(), niter_changes=numeric(), perform_sum=logical(), verbose=character()))

#
# Class remora kernels
#
CLASS_REMORA_KERNELS <- "RemoraKernels"
class_kernels <- setClass(CLASS_REMORA_KERNELS, 
  representation(class_name="character", eigen_values = "list", eigen_vector = "list", clusters = "numeric", cluster_points = "data.frame", points_per_cluster = "list", centroids="matrix", size="numeric"), 
  prototype = list(class_name=character(), eigen_values=list(), eigen_vector=list(), clusters = numeric(), cluster_points = data.frame(), points_per_cluster = list(), centroids = matrix(), size = numeric()))


#
# Class model definition.
#
CLASS_REMORA_MODEL <- "RemoraModel"
class_kernels <- setClass(CLASS_REMORA_MODEL, 
  representation(config=CLASS_REMORA_CONFIGURATION, model="data.frame", lambda="list", kernels="list"), 
  prototype = list(config=new(CLASS_REMORA_CONFIGURATION), model=data.frame(), lambda=list(), kernels=list()))
    
    
#
# Remora weighting functions available
#
# 0 Euclidean (default)
# 1 One Minus   
# 2 One Minus Squared
# 3 Mahalanobis
# 4 Exp One Minus
# 5 Exp One Minus Sq
# 6 Exp One Log
# 7 - Unimplemented
# 8 Normalized Difference
# 9 Normalized Difference Sq
FUNCTION_REMORA_EUCLIDEAN <- "euclidean"
FUNCTION_REMORA_ONE_MINUS <- "one_minus"
FUNCTION_REMORA_ONE_MINUS_SQ <- "one_minus_sq"
FUNCTION_REMORA_MAHALANOBIS <- "mahalanobis"
FUNCTION_REMORA_EXP_ONE_MINUS <- "exp_one_minus"
FUNCTION_REMORA_EXP_ONE_MINUS_SQ <- "exp_one_minus_sq"
FUNCTION_REMORA_EXP_ONE_LOG <- "exp_one_log"
FUNCTION_REMORA_NORMALIZED_DIFFERENCE <- "normalized_difference"
FUNCTION_REMORA_NORMALIZED_DIFFERENCE_SQ <- "normalized_difference_sq"
FUNCTIONS_REMORA = c(FUNCTION_REMORA_EUCLIDEAN, FUNCTION_REMORA_ONE_MINUS, FUNCTION_REMORA_ONE_MINUS_SQ, 
          FUNCTION_REMORA_MAHALANOBIS, FUNCTION_REMORA_EXP_ONE_MINUS, FUNCTION_REMORA_EXP_ONE_MINUS_SQ, 
          FUNCTION_REMORA_NORMALIZED_DIFFERENCE, FUNCTION_REMORA_NORMALIZED_DIFFERENCE_SQ)


#
# Verbose
#
# No verbose, means silent
VERBOSE_NO <- "no"
# Display some information
VERBOSE_YES <- "yes"
# Displays detailed information
VERBOSE_DETAIL <- "detail"
# Displays debug information
VERBOSE_DEBUG <- "debug"
VERBOSE_OPTIONS <- c(VERBOSE_NO, VERBOSE_YES, VERBOSE_DETAIL, VERBOSE_DEBUG)

#
# Builds a distance object.
# Null values will result in no slot valur creation.
#
# @param index is the point index in the matrix
# @param distance is the distance between the point and the cluster
# @param point is the point data
# @param class_name is the name of the classe to which the point belongs to
# @return result distance in a ClassRemoraDistance object
#
remoraDistance <- function(index, distance = NULL, point = NULL, class_name = NULL) {
  
  result_distance <- new(CLASS_REMORA_DISTANCE)
  result_distance@index <- index
  if (!is.null(distance)) {
    result_distance@distance <- distance
  }
  if (!is.null(point)) {
    result_distance@point <- list(point)
  }
  if (!is.null(class_name)) {
    result_distance@className <- class_name
  }
  
  result_distance
}


#
# Get Matrix Without Classifier.
# Recieves a matrix and splits it by its classes. 
# @param data_matrix: the data matrix
# @param class_column: data matrix class column identifier (name or number)
# @return: data frame without the class column
#   
getUnclassedMatrix <- function(data_matrix, class_column) {
  as.data.frame(data_matrix[ , !(colnames(data_matrix) %in% class_column)])
}


#
# Get Classe Names.
# Recieves a matrix and splits it by its classes. 
# @param data_matrix: the data matrix, tipically the training matrix
# @param class_column: the column identifier (name or number)
# @return: array of classes   
#
getClassNames <- function(data_matrix, class_column) {
  data_matrix <- as.data.frame(data_matrix)
  classes <- split(getUnclassedMatrix(data_matrix,class_column), data_matrix[class_column])
  classes
}

#
# Returns the flat index of a class and its sub item index.
# @param classe_name is the class name 
# @param sub_index is the sub item index
# @return composed flat index
# @see getFlatIndexClass, getHierarquicalIndex
#
getFlatIndex <- function(classes, sub_index) {
  sprintf("%s.%i", classes, sub_index)
}

#
# Returns the class part of a composed flat index.
# @para flatIndex is the flat index
# @return the class part of a composed flat index
# @see getFlatIndex, getHierarquicalIndex 
#
getFlatIndexClass <- function(flatIndex) {
  unlist(strsplit(flatIndex, '\\.'))[1]
}

#
# Returns the hierarquical index of a composed flat index.
# @para flatIndex is the flat index
# @return the hierarquical index of a composed flat index
# @see getFlatIndex, getFlatIndexClass
#
getHierarquicalIndex <- function(flatIndex) {
  unlist(strsplit(flatIndex, '\\.'))
}

#
# Builds a Remora configuration object.
#
# @param number_clusters is the number of clusters to use in the training part
# @param class_name is the name, or index, of the training data matrix class column
# @param weighting_function is the name of the weighting function to use in the classification process
# @param scale_variance specifies if the scale should be performed for the principal components analysis, default is True (@see prcomp)                                                                                                  
# @param s_value is the initial s value to use to find the kernel sigma value
# @param d is the initial d value to use to find the s value
# @param epsilon is the epsilon value for function, for thoese functions that require it
# @param niter is the maximum number of iterations to perform to find s, if no value is provided, a default will be calculated based on the number of training data points
# @param niter_changes is the number of iteration without changes that can occur, if the number of niter_changes is reached without any change, the iteration will stop, a default value will be used if none is specified
# @param perform_sum specifies if the sum of the centroids per cluster should be applied, or not
# @param clustering_algorithm specifies which of the k-means algorithm should be used, if none specified, the default k-means algorithm will be used (@see kmeans)
# @param verbose specifies the algorithm verbosity during it's execution (runtime implementation specific parameter)
#
# @return configuration in a ClassRemoraConfiguration object
#
remoraConfiguration <- function(number_clusters, class_name, weighting_function = FUNCTION_REMORA_EUCLIDEAN, scale_variance=TRUE, s_value = 0.2, d = 0.23, epsilon = 0.01, niter=-1, niter_changes=5, perform_sum = TRUE, clustering_algorithm = "", verbose=VERBOSE_NO) {
  
  result_configuration <- new(CLASS_REMORA_CONFIGURATION)
  
  if (missing(number_clusters)) {
    stop('Number of clusters must be specified.')
  } else {
    result_configuration@number_clusters <- number_clusters
  }
  
  if (missing(class_name)) {
    stop('Class name must be specified.')
  } else{
    result_configuration@class_name <- class_name
  }
  
  result_configuration@weighting_function <- weighting_function  
  result_configuration@epsilon <- as.numeric(epsilon)
  result_configuration@scale_variance <- scale_variance
  result_configuration@niter <- as.numeric(niter)
  result_configuration@niter_changes <- as.numeric(niter_changes)
  result_configuration@perform_sum <- perform_sum  
  result_configuration@clustering_algorithm <- clustering_algorithm  
  result_configuration@verbose <- verbose  
  
  if (as.numeric(s_value) > 0) {
    result_configuration@s <- as.numeric(s_value)
  } else {
   stop("s value must be greater than zero", call. = TRUE)
  }
  
  if (as.numeric(d) > 0) {
    result_configuration@d <- as.numeric(d)
  } else {
    stop("d must be greater than zero", call. = TRUE)
  }
  
#  if (as.numeric(alpha) > 0) {
#    if (as.numeric(alpha) < 1) {
#      result_configuration@alpha <- as.numeric(alpha)
#    } else {
#      stop("alpha must be less than 1 ", call. = TRUE)
#    }
#  } else {
#    stop("alpha must be greater than 0 ", call. = TRUE)
#  }

  if (weighting_function %in% FUNCTIONS_REMORA) {
    result_configuration@weighting_function <- weighting_function
  } else {
    stop("Undefined weighting function '", weighting_function ,"'.", , call. = TRUE)
  }
  
  result_configuration
}


#
# Remora kernels constructor.
# @param class_name: the name of the class
# @param clusters: vector of integers indicating the cluster to which each point is allocated
# @param cluster_points: matrix with the point that belong to the clusters
# @param points_per_cluster: list with cluster points organized by clusters
# @param centroids: matrix of cluster centres
# @param size: number of points in each cluster
# @param eigen_values: the list of cluster eigen values
# @param eigen_vector: the list of cluster eigen vector matrices
#
# @see kmeans, prcomp
#
remoraKernels <- function(class_name, clusters, cluster_points, points_per_cluster, centroids, size, eigen_values, eigen_vector) {
  result_kernel <- new(CLASS_REMORA_KERNELS)

  if (missing(class_name)) {
    stop('Class name must be specified.')
  } else {
    result_kernel@class_name <- class_name
  }

  if (missing(clusters)) {
    stop('Clusters must be specified.')
  } else {
    result_kernel@clusters <- clusters
  }
  
  if (missing(cluster_points)) {
    stop('Cluster points must be specified.')
  } else {
    result_kernel@cluster_points <- cluster_points
  }

  if (missing(points_per_cluster)) {
    stop('Points per cluster must be specified.')
  } else {
    result_kernel@points_per_cluster <- points_per_cluster
  }
  
  if (missing(centroids)) {
    stop('Centroids must be specified.')
  } else {
    result_kernel@centroids <- centroids
  }

  if (missing(size)) {
    stop('Size must be specified.')
  } else {
    result_kernel@size <- size
  }
  
  if (missing(eigen_values)) {
    result_kernel@eigen_values=list()
    for (point_cluster in points_per_cluster) {
      result_kernel@eigen_values[point_cluster] <- integer()
    }
  } else {
    eigen_values=eigen_values
  }
  
  if (missing(eigen_vector)) {
    result_kernel@eigen_vector=list()
    for (point_cluster in points_per_cluster) {
      result_kernel@eigen_values[point_cluster] <- list()
    }
  } else {
    eigen_vector=eigen_vector
  }
  
  result_kernel
}

#
# Builds a model object.
# 
# @param config is the configuration used to build the model
# @param model_matrix is the matrix with the model data information
# @param model_lambda is the function lambda (calculated only once)
# @param model_kernels is the list of the kernels
# @return result model in a ClassRemoraModel object
#
remoraModel <- function(config, model_matrix, model_lambda, model_kernels) {
  
  result <- new(CLASS_REMORA_MODEL)
  result@config <- config
  result@model <- model_matrix
  result@lambda <- model_lambda
  result@kernels <- model_kernels
  
  result
}

#
# Remora kernels constructor for kmeans input
# @param class_name: the name of the class
# @param k_means: kmeans result object
# @param cluster_points: the data points that belong to the clusters
#
# @see kmeans, prcomp
#
remoraKernelsKmeans <- function(class_name, k_means, cluster_points) {
  result_kernel <- remoraKernels(class_name, k_means$cluster, cluster_points, getPointsPerCluster(k_means$cluster, cluster_points), k_means$centers, k_means$size)
  result_kernel
}


#
# Configuration print function.
# @param x configuration to print
#
print.RemoraConfiguration <- function(x, ...) {
  configuration <- x
	
  cat('Configuration')

  cat('\n  class_name (class name): ')  
  cat(configuration@class_name)

  cat('\n  number_clusters (number of clusters): ')
  cat(configuration@number_clusters)
  
  cat('\n  clustering_algorithm (k-means algorithm): ')
  cat(configuration@clustering_algorithm)

  cat('\n  weighting_function (weighting function): ')  
  cat(configuration@weighting_function)
  
  cat('\n  niter (number of iteractions): ')
  cat(configuration@niter)

  cat('\n  niter_changes (number of iteractions without change allowed): ')
  cat(configuration@niter_changes)

  cat('\n  perform_sum (perform cluster sum): ')
  cat(configuration@perform_sum)

  cat('\n  scale_variance (scale variance): ')
  cat(configuration@scale_variance)

  cat('\n  s: ')
  cat(configuration@s)

  cat('\n  d: ')
  cat(configuration@d)

  #cat('\n  alpha: ')
  #cat(configuration@alpha)

  cat('\n  epsilon: ')
  cat(configuration@epsilon)

  cat('\n  verbose: ')
  cat(configuration@verbose)
  
  cat('\n')
}


#
# Print function for remora kernels.
# @param x the kernels to print
#
print.RemoraKernels <- function(x, ...) {
  cat('Kernels')
  kernels <- x

  cat('\n  class_name: ')  
  cat(kernels@class_name)

  cat('\n  eigen_values:\n') 
  for (ev in kernels@eigen_values) {
	cat('    ')
	cat(ev)
	cat('\n')
  }

  cat('\n  eigen_vector:\n') 
  for (ev in kernels@eigen_vector) {
	cat('    ')
	cat(ev)
	cat('\n')
  }

  cat('\n  centroids:\n')  
  print(kernels@centroids)

  cat('\n  size: ')  
  cat(kernels@size)

  cat('\n\n  clusters:\n')  
  print(kernels@clusters)

  cat('\n\n cluster_points:\n')  
  print(kernels@cluster_points)

  cat('\n  points_per_cluster:\n') 
  for (ppc in kernels@points_per_cluster) {
	print(ppc)
	cat('\n')
  }

  cat('\n')
}

#
# Print function for remora model.
# @param x the model to print
#
print.RemoraModel <- function(x, ...) {
  model <- x
  print(model@config)
  cat('\nModel Table:\n')
  print(model@model)

  cat('\nLambda:\n')
  for(lmd in model@lambda) {
    print(lmd)
  }
  cat('\n')
  print(model@kernels)
}