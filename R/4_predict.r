#
# Remora Prediction
#

# 
# Version 1, September 2009
# Fernando Martins
# fmp.martins@gmail.com
# http://www.vilma-fernando.net/fernando
#

#
# Builds the classification for each point through the distance matrix. 
# It can be done using a class sum centroid distances or just by using the centroid distance.
#
# @param distance_matrix is the matrix with the distances between each point and each kernel
# @param config is the configuration
# @returns list of points and respective class to which it belongs
#
buildClassification <- function(distance_matrix, config) {
  col_names <- colnames(distance_matrix)
  result <- list()

  if (config@perform_sum) {             
    # perform the centroid class sum    
    sum_matrix <- data.frame(distance_matrix)
    colnames(sum_matrix) <- colnames(distance_matrix)
    for (class_name in col_names) {
      #fullIndex <- getHierarquicalIndex(class_name)
      fullIndex <- list()
      fullIndex[1] <- getFlatIndexClass(class_name)
      fullIndex[2] <- getFlatIndexCentroid(class_name)
      if (fullIndex[2] > 1) {
        sum_matrix[getFlatIndex(fullIndex[1], 1)] <- sum_matrix[getFlatIndex(fullIndex[1], 1)] + sum_matrix[class_name]
        sum_matrix <- getUnclassedMatrix(sum_matrix, class_name)
      }
    }    
  } else {
    # no need to perform the centroid class sum, so just simulate it
    sum_matrix <- distance_matrix                            
  }
  
  # check the distance over the summed values
  result2 <- max.col(sum_matrix)
  result <- list()
  for (class_index in c(1:length(result2))) {
    result[class_index] <- col_names[result2[class_index]]
  }
  names(result) <- names(rownames(sum_matrix))
  
#  cat('buildClassification | result')
#  browser()
  result
}



#
# Builds the classification distance table for the data matrix points.
# @param remora_model is the remora model to use
# @param data_matrix is the data matrix that holds the points to classify
# @return distance table
#
classificationDistanceTable <- function(remora_model, data_matrix) {
  config <- remora_model@config
  number_points <- nrow(data_matrix)
  model = remora_model@model
  model_lambda = remora_model@lambda
  class_names <- unique(model$class)
  kernels <- remora_model@kernels
  kernels_s <- list()

  data_distance_table <- buildDistanceTable(as.matrix(data_matrix), kernels, model_lambda, config)
  
  for (current_class in class_names) {
    for (current_index in model[ model["class"] == current_class, "centroid"]) {
      flatIdx <- getFlatIndex(current_class, current_index)
      # get s values
      kernels_s[flatIdx] <- as.numeric(model[ model["class"] == current_class &  model["centroid"] == current_index, "s"])
    }
  }
  data_distance_table <- distances(data_distance_table, kernels, kernels_s, config)  # get distances using s values
  data_distance_table
}



#
# Predicts and the classifies the information.
#
# @param model: Remora model
# @param data_matrix: data matrix to classify
# @return: Remora classification prediction
#
# @see remora.predict function
#
remora.classify <- function(model, data_matrix, return_nice=FALSE) {

  prediction <- remora.predict(model, data_matrix)
  #classification <- prediction@prediction

  result_matrix <- as.data.frame(data_matrix)
  class_column <- ncol(data_matrix) + 1
  for(point_index in c(1:length(prediction))) {  
    result_matrix[point_index, class_column] = prediction[point_index]
  }
  names(result_matrix) <- c(names(data_matrix), "class")

#  if (verbose.showDebug(model@config@verbose)) {
#    cat('\nClassification predicted.\nPrediction is:\n')
#    print(result_matrix)
#  }
    
  if (return_nice) {
    result_matrix
  } else {
    prediction
  }
}
