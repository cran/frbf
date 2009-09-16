#
# Utilities
#

# 
# Version 1, September 2009
# Fernando Martins
# fmp.martins@gmail.com
# http://www.vilma-fernando.net/fernando
#

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
# Returns the flat index of a class and its centroid index.
# @param classe_name is the class name 
# @param sub_index is the sub item index
# @return composed flat index
# @see getFlatIndexClass, getFlatIndexCentroid
#
getFlatIndex <- function(classes, sub_index) {
  sprintf("%s.%i", classes, sub_index)
}

#
# Returns the class part of a composed flat index.
# @para flatIndex is the flat index
# @return the class part of a composed flat index
# @see getFlatIndex, getFlatIndexCentroid 
#
getFlatIndexClass <- function(flatIndex) {
  #unlist(strsplit(flatIndex, '\\.'))[1]
  result <- unlist(strsplit(flatIndex, '\\.'))
  if (length(result) > 2) {    
    for(i in c(1:(length(result)-2)+1)) {
      result[1] <- sprintf("%s.%s", result[1], result[i])
    }
  }
  result[1]
}

#
# Returns the centroid part of a composed flat index.
# @para flatIndex is the flat index
# @return the centroid part of a composed flat index
# @see getFlatIndex, getFlatIndexClass 
#
getFlatIndexCentroid <- function(flatIndex) {
  #unlist(strsplit(flatIndex, '\\.'))[1]
  result <- unlist(strsplit(flatIndex, '\\.'))
  result[length(result)]
}

#
# Returns the hierarquical index of a composed flat index.
# @para flatIndex is the flat index
# @return the hierarquical index of a composed flat index
# @see getFlatIndex, getFlatIndexClass
#
#getHierarquicalIndex <- function(flatIndex) {
#  unlist(strsplit(flatIndex, '\\.'))
#}


#
# Calculates the accuracy of a matrix
# @param accuracy_matrix data matrix with accurate classification
# @param classification_vector the classification to verify accuracy
# @param config the remora configuration 
# @param return_hit_vector if true it will return accuracy ratio per cluster 
#        and total hit ration, otherwise it will return the total number of hits
# @return total number of hits or accuracy ratio per cluster and total hit ration 
# 
accuracy <- function(accuracy_matrix, classification_vector, config, return_hit_vector = TRUE) {  

  TOTAL <- 'Total'
  #class_names <- names(classification_vector)
  class_column <- config@class_name
  cluster_column <- ncol(accuracy_matrix)
  total_vector <- length(classification_vector)
  total_correct <- 0
  
  points_class <- list()  
  cluster_correct_classify <- list()  
  cluster_hit_ratio <- list()
  flat_index <- list()
 
  for (point in c(1:nrow(accuracy_matrix))) {
    class_name <- accuracy_matrix[point, class_column]
    class_cluster <- accuracy_matrix[point, cluster_column]
    flatIndex <- getFlatIndex(class_name, class_cluster)
    flat_index[flatIndex] <- flatIndex
    if (is.null(points_class[[flatIndex]])) {
      points_class[flatIndex] <- 1
    } else {
      points_class[flatIndex] <- as.numeric(points_class[flatIndex]) + 1
    }
    cluster_correct_classify[flatIndex] <- 0
    cluster_hit_ratio[flatIndex] <- 0
  }
  cluster_hit_ratio[TOTAL] <- 0

  # holds the number of correct classifications per cluster
  for (index in c(1:total_vector)) {
    centroid_class <- classification_vector[[index]]
    #full_index <- getHierarquicalIndex(centroid_class)
    full_index <- list()
    full_index[1] <- getFlatIndexClass(centroid_class)
    full_index[2] <- getFlatIndexCentroid(centroid_class)
        
    if (accuracy_matrix[index, class_column] == full_index[1] && accuracy_matrix[index, cluster_column] == full_index[2]) {
      cluster_correct_classify[centroid_class] <- as.numeric(cluster_correct_classify[centroid_class]) + 1
      total_correct <- total_correct + 1
    }
  }

  # holds the hit ration per cluster
  for (cluster_name in flat_index) {
      cluster_hit_ratio[cluster_name] <- as.numeric(cluster_correct_classify[cluster_name]) / as.numeric(points_class[cluster_name])      
  }
  cluster_hit_ratio[TOTAL] <- total_correct / nrow(accuracy_matrix)
  
  if (return_hit_vector) {
    cluster_hit_ratio
  } else {
    total_correct
  }
}


#
# Difference between each centroid and each point.
#
# @param distance_values distance table with precalculated distance point-centroid values 
# @param kernels the kernels (found on stage one)
# @param kernels_s the s value per each kernel
# @param config the Remora algorithm configuration
# @return distance table between points and classes, matrix row number is the point row index
#
distances <- function(distance_values, kernels, kernels_s, config) {
  class_names <- colnames(distance_values)
#  total_data_points <- nrow(distance_values)
#  total_clusters <- ncol(distance_values)
#
#  distance_table <- matrix(nrow = total_data_points, ncol= total_clusters, byrow=TRUE, dimnames = list(c(1:total_data_points), class_names))
#  
#  cat('\nTests...\n')
#  # iterates for all data rows (each row is a point)
#  time_starta <- unclass(Sys.time())
#  for (data_column in class_names) {
#    #flat_index <- data_column
#    full_index <- getHierarquicalIndex(data_column)
#    class_name <- full_index[1]
#    centroid_index <- as.numeric(full_index[2])
#    
#    k <- kernels[[class_name]]
#      
#    kernel_s <- as.numeric(kernels_s[data_column])
#        
#    # wji é o total do tamanho do centroid em análise
#    wji <- k@size[centroid_index] # the number of patterns that are included in each cluster.
#       
#    for (data_row in c(1:total_data_points)) {
#      # formula (4) 
#      # sigma %*% pT %*% fn_lambda %*% p %*% sigmaT
#      # p * fn_lambda * pT has been pre calculated in fn_lambda        
#      # pmpT has been previously calculated
#      pmpT <- distance_values[data_row, data_column]
#
#      kij <- exp(-1 * pmpT * kernel_s)
#      point_distance <- wji * kij
#
#      # keep distance value in the distance table
#      distance_table[data_row, data_column] <- point_distance      
#    } 
#  }  
#  time_enda <- unclass(Sys.time())
#  cat(time_enda-time_starta)
#  cat(' seconds to calculate the distances\n')    
#  print(distance_table)
  
  # TODO: TESTAR COM A CRIAÇÃO DE UMA TABELA NOVA
  distance_table <- distance_values
#  time_startb <- unclass(Sys.time())
   for (data_column in class_names) {
    #flat_index <- data_column
    #full_index <- getHierarquicalIndex(data_column)
    #class_name <- full_index[1]
    #centroid_index <- as.numeric(full_index[2])
    class_name <- getFlatIndexClass(data_column)
    centroid_index <- as.numeric(getFlatIndexCentroid(data_column))
    
    k <- kernels[[class_name]]
    kernel_s <- as.numeric(kernels_s[data_column])
        
    # wji é o total do tamanho do centroid em análise
    #browser()
    wji <- k@size[centroid_index] # the number of patterns that are included in each cluster.
    distance_table[,data_column] = wji * exp(-distance_table[,data_column] * kernel_s)
  }
#  time_endb <- unclass(Sys.time())
#  cat(' [')
#  cat(time_endb-time_startb)
#  cat(' seconds] ')
#  print(distance_table)
  
  # return the distances
#  cat('distances | check distance_table\n')
#  browser()
    
  as.matrix(distance_table)
}

#
# Shows the accuracy
#
showAccuracy <- function(classification, correct_classification) {                                  
    size <- length(correct_classification)                                      
    tot <- sum (classification == correct_classification)
#    cat('Classification Hits:\n')
#    tot <- 0
#    for (idx in c(1:size)) {
#      if (classification[idx] == correct_classification[idx]) {
#        tot <- tot + 1
#        cat('1 ')
#      } else {
#        cat('0 ')
#      }
#    }
    cat('Accuracy: ')
    cat(tot/size)
    cat('\n')
}

#
# Checks vervosity is set
#
verbose.show <- function(verbosity){  
  if (verbosity == VERBOSE_NO) {
    FALSE
  } else {
    TRUE
  }
}

#
# Checks if the vervosity is detail
#
verbose.showDetail <- function(verbosity){  
  if (verbosity == VERBOSE_NO) {
    FALSE
  } else if (verbosity == VERBOSE_YES) {
    FALSE
  } else {
    TRUE
  }
}

#
# Checks if the vervosity is debug
#
verbose.showDebug <- function(verbosity){  
  if (verbosity == VERBOSE_DEBUG) {
    TRUE
  } else {
    FALSE
  }

}
