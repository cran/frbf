#
# Main functions
#


#
# Remora model function.
# The learning procedure of the frbf algorithm.
# 
# @param data_matrix the data to use to train the algorithm 
# @param number_clusters is the number of clusters to use in the training part
# @param class_name is the name, or index, of the training data matrix class column
# @param weighting_function is the name of the weighting function to use in the classification process
# @param scale_variance specifies if the scale should be performed for the principal components analysis, default is True (@see prcomp)
# @param s_value is the initial s value to use to find the kernel sigma value
# @param d is the initial d value to use to find the s value
# @param epsilon is the epsilon value for function, only for functions that require it
# @param niter is the maximum number of iterations to perform to find s, if no value is provided, a default will be calculated based on the number of training data points
# @param niter_changes is the number of iteration without changes that can occur, if the number of niter_changes is reached without any change, the iteration will stop, a default value will be used if none is specified
# @param perform_sum specifies if the sum of the centroids per cluster should be applied, or not
# @param clustering_algorithm specifies which of the k-means algorithm should be used, if none specified, the default k-means algorithm will be used (@see kmeans)
# @param verbose specifies the algorithm verbosity during it's execution (runtime implementation specific parameter)
# 
# @return model
#
frbf <- function(data_matrix, number_clusters, class_name, weighting_function = "euclidean", scale_variance=TRUE, s_value = 0.2, d = 0.23, epsilon = 0.01, niter=-1, niter_changes=5, perform_sum = TRUE, clustering_algorithm = '', verbose="no") {
  # train
  if (verbose.show(verbose)) {
    cat('Model phase...\n')
    flush.console()
  }
  config <- remoraConfiguration(number_clusters, class_name, weighting_function, scale_variance, s_value, d, epsilon, niter, niter_changes, perform_sum, clustering_algorithm, verbose)
  model <- remora.model(data_matrix, config)
  model
}

#
# Remora predict function.
#
# @param object: remora model, obtained from the learning procedure
# @param data_matrix: the data to classify
# @return prediction
#
# @see frbf
#
predict.RemoraModel <- function(object, data_matrix, ...) {
  model <- object
  # classify
  if (verbose.show(model@config@verbose)) {
    cat('Classification phase...\n')
    flush.console()
  }
  data_matrix <- as.matrix(getUnclassedMatrix(data_matrix, model@config@class_name))
  if (verbose.showDebug(model@config@verbose)) {
    classification <- remora.predict(model, data_matrix)
  } else {
    classification <- remora.classify(model, data_matrix)  
  }    
  classification    
}


##
## Auto test function.
##
#auto_test <- function () {  
#  cat('Auto testing frbf...\n')
#
#  data(iris)
#  sample_size <- (5/6) * nrow(iris)
#  sample_idx <- sample(c(1:nrow(iris)), sample_size)
#  training_matrix <- iris[sample_idx,]
#  classify_matrix <- iris[-sample_idx,]
#  matrix_class_name <- "Species"  
#  cat('\tUsing infert as control data...\n')
#
#  # Phase 1: Training
#  model <- frbf(training_matrix, class_name = matrix_class_name, number_clusters = 3, verbose="debug")
#  cat('\n')
#  print(model@config)
#  cat('\n')
#  #print(model)
#
#  # Phase 2: Classification  
#  prediction <- predict(model, classify_matrix)
#
#
#  # Phase 3: Analysis
#  classe_names <- as.vector(unique(classify_matrix[,matrix_class_name]))
#  table_cross <- table(classify_matrix[,matrix_class_name], prediction, dnn=c('Iris Classes\t|', '\t|  Prediction')) 
#  cat('\n')
#  print(table_cross)
#  total_hit <- 0
#  cat('\nAccuracy:\n')
#  for(cls in classe_names) {
#    accuracy <- 0
#	cat('     ')
#	cat(cls)
#    cat(': ')
#	if (cls %in% rownames(table_cross)) {
#		total_class <- sum(table_cross[cls, ])
#		if (cls %in% colnames(table_cross)) {
#			hit <- table_cross[cls, cls]
#			accuracy <- hit/total_class
#			total_hit <- total_hit + hit
#		}
#	}
#    cat(accuracy * 100)
#    cat('%\n')
#  }
#  cat('     Total: ')
#  cat((total_hit / nrow(classify_matrix)) * 100)
#  cat('%\n')
#
#}
