#
# Remora
#
# Phase 1 - Model
#
# The number of kernels in each class is determined by the proportion of the 
# total variance in the training set associated with that class.
# Training data are separated according to each class, and the kernel�s centers 
# are determined using a clustering algorithm (K-means). 
# This algorithm operates over the normalized variables on each class.
# After the clustering stage, the kernel location parameters correspond 
# to the centroids of the resulting clusters
#


#
# Get Number of Clusters.
# Retrieves the number of clusters depending on 
# each class variance.
#
# @param training_matrix: the training data matrix
# @param classes: array of classes 
# @param clusters: configuration
# @return array of clusters per class 
#
getNumberOfClusters <- function(training_matrix, classes, config) {  
  clusters <- config@number_clusters
  class_column <- config@class_name
  #variance_matrix <- sd(getUnclassedMatrix(training_matrix, class_column))^2
  total_classes <- length(classes)
  variances <- c(1:total_classes)
  
  for (i in 1:total_classes) {
    if (nrow(classes[[i]]) < 2) {
      stop('Trainig class "',names(classes[i]),'" must have at least two data points.')
    }
    variances[i] <- sum((sd(classes[[i]])^2))
  }
  total_variance <- sum(variances)
  kernels <- c(1:total_classes)
  
  total_kernels <- 0
  for (i in 1:total_classes) {
    kernels[i] <- round((variances[i] * clusters) / total_variance)
    if (kernels[i] < 1) {
      kernels[i] <- 1
    }
    total_kernels <- total_kernels + kernels[i]
    names(kernels)[i] <- names(classes[i])    
  }  

  kernels <- adjustKernelClusters(kernels, clusters, config)
  kernels
}



#
# Adjust Kernel Clusters.
# Adjustes the number of kernel clusters to fit the number of clusters 
# specified by the user. When that is not possible, the number of kernel
# clusters will increase to force at least one cluster per kernel.
# This function issues warnings when adjustments have been 
# made to the kernel clusters.
#
# @param kernels the kernel clusters
# @param clusters the number of desired clusters
# @param config configuration
# @return kernel clusters adjusted
#
adjustKernelClusters <- function(kernels, clusters, config) {
  total_kernels <- sum(kernels)

  if (total_kernels > clusters) {	
  	while ((total_kernels > clusters) && (max(kernels) > 1)) { 
  		for (i in 1:length(kernels)) {
  			if (kernels[i] == max(kernels)) {
  				kernels[i] <- kernels[i] - 1
  				break;
  			}
  		}
  		total_kernels <- sum(kernels)
  	}
  	if (total_kernels > clusters) {	
  		# It was not possible to adjust to the number of clusters defined by the user, so the number of clusters has been updated
  		warning("Insuficient number of clusters. Number of clusters has been updated to a total of ", total_kernels, " clusters.", call. = TRUE)
  		if (verbose.show(config@verbose)) {
  		   cat("\tInsuficient number of clusters. Number of clusters has been updated to a total of ", total_kernels, " clusters.\n")
  		   flush.console()
      }
  	} else {
  		# Kernels have been adjusted to fit the number of clusters defined by the user
  		warning('Cluster distribution has been ajusted to fit the defined number of clusters.', call. = TRUE)
  		if (verbose.show(config@verbose)) {
  		   cat("\tCluster distribution has been ajusted to fit the defined number of clusters.\n")
  		   flush.console()
      }
  	}
  }

  kernels
}


#
# Get KMeans returns the kernels.
# Performs the k-means calculations for each class.
#
# @param training_matrix: the training data matrix
# @param config: the Remora algorithm configuration
# @return list of k-means calculations
#
getKMeans <- function(training_matrix, config) {
  if (verbose.showDetail(config@verbose)) {
    cat('\tPerforming Kmeans... ')      
    time_starta <- unclass(Sys.time())
    flush.console()
  }
  
  class_column <- config@class_name
  number_clusters <- config@number_clusters
  algorithm <- config@clustering_algorithm
  classes <- names(number_clusters)
  total_classes = length(classes)  
  k_means <- list()
  
  if (algorithm == "") {
    # apply kmeans default algorithm
    algorithm <- NULL
  }

  
  for (i in 1:total_classes) {
    d <- training_matrix[training_matrix[class_column] == classes[i],]
    d <- getUnclassedMatrix(d, class_column)

    n_row <- nrow(d)
    if (number_clusters[i] > 1) {      
      # maximum number of iterations
      iter.max <- n_row %/% 2
      iter.max <- ifelse(iter.max < 10, 10, iter.max)
  
      # number of random sets to choose
      nstart <- n_row %/% 4
      nstart <- ifelse(nstart < 4, 4, nstart)
      
      # K-Means 
      k <- kmeans(d, number_clusters[i], iter.max, nstart, algorithm)
      k_means[i] <- remoraKernelsKmeans(classes[i],k,d)
      
    } else {
      # build a remora kernel from a dummy kmeans result
      # simulate a kmeans result
      dummy_cluster <- vector()
      for (v in row.names(d)) {
        dummy_cluster[v] <- 1
      }      
      names(dummy_cluster) <- row.names(d)
    
      dummy_centers <- matrix(colMeans(d), ncol=length(names(d)), nrow=1, dimnames=list(1, names(d)))
      
      # simulate points per cluster   
      points_result <- getPointsPerCluster(dummy_cluster, d)
       
      k_means[i] <- remoraKernels(classes[i], dummy_cluster, d, points_result, dummy_centers, n_row)
    }
  }
  
  if (verbose.showDebug(config@verbose)) {
    time_enda <- unclass(Sys.time())
    cat('(')
    cat(time_enda-time_starta)
    cat(' seconds)')
  } 
  if (verbose.show(config@verbose)) {
    cat('\n')
    flush.console()
  }

  names(k_means) <- classes
  k_means
}

#
# Gets a list of points per cluster, where each cluster holds 
# the corresponfing point.
# @param point_cluster: list points and their cluster bind
# @param data_cluster: data frame with the data points
# @retuns list of cluster and their points
#
getPointsPerCluster <- function(point_cluster, data_points) {
  
  result <- list()
  
  for (point in unique(point_cluster)) {
    aux <- data_points[names(point_cluster[point_cluster == point]), ]
    names(aux) <- names(data_points)
    aux[,"point_index"] = rownames(aux)
    rownames(aux) <- c(1:length(rownames(aux)))
    aux <- aux[, c("point_index", names(data_points))]
    result[point] <- list(aux)
  }
  
  result
  
#  points_result <- list()
#  points_size <- length(point_cluster)
#  frame_names <- c( "point_index", names(data_points))
#  for (point in c(1:points_size)) {
#    
#    points_id <- names(point_cluster[point_cluster == point])        
#    if (length(points_id) > 0) {
#      points_data <- data.frame()
#      row_id <- 0
#      for (id in points_id) {
#        row_id <- row_id + 1
#        points_data[row_id, 1] <- id
#        for (column in names(data_points)) {
#          points_data[row_id, column] <- data_points[id, column]
#        }
#      }
#      names(points_data) <- frame_names
#      points_result[point] <- list(points_data)
#    }
#  }
#
#  cat('points_result')
#  browser()
#  points_result
}

#
# Finds the Principal Components recurring to PCA.
#
# @param kernels is the kernels 
# @param config is the configuration
# @return principal component analysis
#
getPCA <- function(kernels, config) {
  if (verbose.showDetail(config@verbose)) {
    cat('\tPerforming PCAs...\n')  
    flush.console()
  }
  
  pca_result <- list()

  # iterates for clusters/kernels
  for (cluster_name in names(kernels)) {
      
    # get cluster specific information
    k <- kernels[[cluster_name]]
    points_per_cluster <- k@points_per_cluster
    for (cluster_id in c(1:length(points_per_cluster))) {
      cluster_points <- points_per_cluster[[cluster_id]]      
      
      #
      # TEST WITH CENTER = TRUE
      #
      pca <- prcomp(cluster_points[, c(2:length(cluster_points))], scale. = config@scale_variance)
      stddev  <- pca$sdev
      rotation <- pca$rotation
      k@eigen_values[cluster_id] <- list(stddev)
      k@eigen_vector[cluster_id] <- list(rotation)
    }
    pca_result[cluster_name] <- k
  }

  pca_result
}


#
# Trainging step.
#
# @param training_matrix: the training data matrix
# @param config: the Remora algorithm configuration
# @returns kernels from the training step
#
train <- function(training_matrix, config) {
  
  training_matrix <- as.data.frame(training_matrix)
  classes <- getClassNames(training_matrix, config@class_name)
  config@number_clusters <- getNumberOfClusters(training_matrix, classes, config)
  kernels <- getKMeans(training_matrix, config)
  remora_kernels <- getPCA(kernels, config)
    
  remora_kernels
}


#
# Initialize s values.
# Starts with a user given s value and modifies it until the results start getting worst.
#
# @param distance_table: distance table with precalculated distance point-centroid values
# @param kernels: the kernel clusters
# @param accuracy_matrix: the matrix with the accurate classification
# @param config: the configuration
#
# @return initial s values 
#
initS <- function(distance_table, kernels, accuracy_matrix, config) {
  if (verbose.showDebug(config@verbose)) {
    cat('\t\tInitializing sigmas...')
  }
  sigma <- c(1:length(colnames(distance_table)))
  names(sigma) <- colnames(distance_table)
  for (classes in colnames(distance_table)) {
    sigma[classes] <- as.numeric(config@s)
  }
  
  # get distances for sigma values
  dst_last <- distances(distance_table, kernels, sigma, config)  
  # classification for the distances found
  cls_last <- buildClassification(dst_last, config)              
  # accuracy hit
  res2 <- accuracy(accuracy_matrix, cls_last, config, FALSE)     
  init_sigma <- sigma
  if (verbose.showDebug(config@verbose)) {
    cat('\n\tFinding initial kernel sigmas.\n\tInitial sigma: ')
    cat(sigma)
    cat(' got ')
    cat(res2);    
    cat(' hits \n')
  }    
  repeat {    
    init_sigma <- sigma
    res1 <- res2
    sigma <- sigma + 0.1
    # get distances for sigma values
    hit_current <- distances(distance_table, kernels, sigma, config)  
    # classification for the distances found
    cls_current <- buildClassification(hit_current, config)           
    # accuracy hit
    res2 <- accuracy(accuracy_matrix, cls_current, config, FALSE)     
   
    if (verbose.showDebug(config@verbose)) {
      cat('\tinit kernel sigma ')
      cat(sigma)
      cat(' got ')
      cat(res2);          
      cat(' hits\n')
    }    

    if (res2 <= res1) {
      break
    }
  }
  
  if (verbose.showDebug(config@verbose)) {
    cat('\n')
    flush.console()
  }
  init_sigma
}

#
# Finds the s value for each cluster.
#
# @param training_matrix training data matrix
# @param kernels the kernels (found on stage one)
# @param model_lambda the previously calculated lambda function values
# @param config the remora configuration 
# @return s parameter per cluster
#
findS <- function(training_matrix, kernels, model_lambda, config) {
  if (verbose.showDetail(config@verbose)) {
    cat('\tFinding S values:\n')  
    flush.console()
  }
  
  number_points <- nrow(training_matrix)
  test_matrix <- as.matrix(getUnclassedMatrix(training_matrix, config@class_name))
  class_names <- names(kernels)
  best_kernels_s <- list()
  test_kernels_s <- list()
  kernels_s_up <- list()
  kernels_s_down <- list()
  
  # if necessary calculates the niter parameter
  if (config@niter < 0) {
    config@niter <- round(nrow(training_matrix) * 0.05)
    warning("niter parameter has been calculated, value is ", config@niter)
  }
  # if necessary adjusts the niter parameter to a minimum value
  if (config@niter < 10) {    
    config@niter <- 10
    warning('niter was too low, it has been redefined to ', config@niter)
  }

  # initialize values
  distance_table <- buildDistanceTable(test_matrix, kernels, model_lambda, config)
  d <- config@d
  d_start <- d # prototype uses 0.23
  d_end <- 0.01
    
  iter <- 1
  random_index = vector()
  flat_index <- list()
  # structure is passed into a flat structure and point cluster index is added
  accuracy_matrix <- training_matrix
  new_cluster_column_index <- ncol(training_matrix) + 1
  for (classes in sample(names(model_lambda))) {
    for (lambda in sample(c(1:length(model_lambda[[classes]])))) {
      flatIndex <- getFlatIndex(classes, lambda)
      flat_index[flatIndex] <- flatIndex

      points_per_cluster <- kernels[[classes]]@points_per_cluster
      for (ppc_idx in (c(1:length(points_per_cluster)))) {
        ppc <- points_per_cluster[[ppc_idx]]
        for (idx in ppc["point_index"]) {
          accuracy_matrix[idx, new_cluster_column_index] <- ppc_idx;
        }
      }
    }
  }
  names(accuracy_matrix) <- c(names(training_matrix),"cluster_index")  
  best_kernels_s <- initS(distance_table, kernels, accuracy_matrix, config)
  kernels_s_up <- best_kernels_s
  kernels_s_down <- best_kernels_s
  last_change <- config@niter_changes
  best_hit <- 0  
  
  if (verbose.showDetail(config@verbose)) {
    cat('\t\tMaximum number of iteractions: ')
    cat(config@niter)
    cat('\n\t\tIteractions: ')  
    flush.console()
  }
  
  # iter
  
  time_start <- unclass(Sys.time())  
  while (last_change >= 0 && iter < config@niter) {  
    if (verbose.showDebug(config@verbose)) {
	    cat('\n\t\t#')
      cat(iter)
      cat(':\n')	
      flush.console()
    }  
	  time_start_iter <- unclass(Sys.time())
    # iterate over random indexs
    for (index in sample(flat_index)) {
      #       			
      # try s up      
      #
      test_kernels_s <- best_kernels_s
      kernels_s_up[index] <- as.numeric(kernels_s_up[index]) * (1 + d)
      test_kernels_s[index] <- kernels_s_up[index]
      
      # get distances for new s up value
      dst_up <- distances(distance_table, kernels, test_kernels_s, config)
              
      # classification for the distances found
      class_up <- buildClassification(dst_up, config) 

      # accuracy for the distances
      hit <- accuracy(accuracy_matrix, class_up, config, FALSE)             
  
      if (hit > best_hit) {
        # new, better s value found
    		if (verbose.showDebug(config@verbose)) {
    			cat('\t\t\ts upper value ')
    			cat(test_kernels_s[[index]])
    			cat(' is better for ')
    			cat(index)
    			cat(' with hit of ')
    			cat(hit)
    			cat(';\n')
    			flush.console()
    		}
        best_kernels_s[index] <- test_kernels_s[index]
        best_hit <- hit        
        last_change <- config@niter_changes
      } 

      #
      # try s down
      #
      test_kernels_s <- best_kernels_s
      kernels_s_down[index] <- as.numeric(kernels_s_down[index]) * (1 - d)
      test_kernels_s[index] <- kernels_s_down[index]      
      
      # get distances for new s up value
      dst_down <- distances(distance_table, kernels, test_kernels_s, config)  
      
      # classification for the distances found      
      class_down <- buildClassification(dst_down, config)                     
      
      # accuracy for the distances
      hit <- accuracy(accuracy_matrix, class_down, config, FALSE)                   
	    
      if (hit > best_hit) {
        # new, better s value found
    		if (verbose.showDebug(config@verbose)) {
    			cat('\t\t\ts lower value ')
    			cat(test_kernels_s[[index]])
    			cat(' is better for ')
    			cat(index) 			
     			cat(' with hit of ')
    			cat(hit)
    			cat(';\n')
    			flush.console()
    		}
        best_kernels_s[index] <- test_kernels_s[index]
        best_hit <- hit
        last_change <- config@niter_changes
      }
      
    }
    d <- d_start+iter/config@niter*(d_end-d_start);
    if (d < d_end) {
      cat("\t\tunexpected d < d_end\n")
      d <- d_end
    }
    last_change <- last_change - 1  
    iter <- iter + 1
	  time_end_iter <- unclass(Sys.time())
    if (verbose.showDebug(config@verbose)) {
      cat('\t\t#')
      cat(iter-1)
      cat(': ')
	    cat(time_end_iter - time_start_iter)
	    cat(' seconds\n')
      flush.console()
    }
  }  
  
  time_end <- unclass(Sys.time())

  if (verbose.showDetail(config@verbose)) {
	iter_time <- time_end - time_start
    cat('\n')
    cat('\tNumber of iteractions performed: ')
    cat(iter-1)
    cat('\n')
  	cat('\tTime spent: ')
  	cat(iter_time/60)
  	cat(' minutes.\n')
  	cat('\tAverage time per iteraction: ')
  	cat(iter_time/iter)
  	cat(' seconds.\n')
  	cat('\tFinal sigma values:\n')
  	print(best_kernels_s)
    flush.console()
  }
  
  best_kernels_s
}


#
# Build the table distance between each centroid and each point.
#
# @param data_matrix data matrix to classify
# @param kernels the kernels (found on stage one)
# @param model_lambda the previously calculated lambda function values
# @param config the Remora algorithm configuration
# @return distance table between points and classes, matrix row number is the point row index
#
buildDistanceTable <- function(data_matrix, kernels, model_lambda, config) {
  number_points <- nrow(data_matrix)
  class_names <- names(kernels)
  flat_index <- list()
  
  # iterates for clusters/kernels
  for (cluster_name in class_names) {
    k <- kernels[[cluster_name]]
    centroids <- k@centroids          
    for (centroid_index in c(1:nrow(k@centroids))) {
      flatIndex <- getFlatIndex(cluster_name, centroid_index)
      flat_index[flatIndex] <- flatIndex
    }
  }

  distance_table <- matrix(nrow = number_points, ncol=length(flat_index), byrow=TRUE, dimnames = list(c(1:number_points), names(flat_index)))

  # iterates for all data rows (each row is a point)
  for (data_row in c(1:number_points)) {
    # extract data point information
    point <- data_matrix[data_row,]
     
    # iterates for clusters/kernels
    for (cluster_name in class_names) {
      k <- kernels[[cluster_name]]
      centroids <- k@centroids      
      
      for (centroid_index in c(1:nrow(centroids))) {      
        # get cluster specific information        
        flatIndex <- getFlatIndex(cluster_name, centroid_index)
    
        # perform principal components analysis
        # all points that are part of this cluster
        fn_lambda <- model_lambda[[cluster_name]][[centroid_index]]

        # exatract cluster specific information        
        centroid <- k@centroids[centroid_index, ]

        # calculate formula values                
        sigma <- t(as.array(point - centroid))
        sigmaT <- t(sigma)                    
        
        # formula (4) result
        # sigma %*% pT %*% fn_lambda %*% p %*% sigmaT
        # p * fn_lambda * pT has been pre calculated in fn_lambda        
        pmpT <- sigma %*% fn_lambda %*% sigmaT
        point_distance <- pmpT

        # include distance value in the distance table
        distance_table[data_row, flatIndex] <- point_distance        
      }
    } 
  }  
  
  # return the distances    
  as.matrix(distance_table)
}

#
# Calculates the lambda function values, 
# i.e. the kernel function values. 
#
# @param training_matrix: training data matrix
# @param kernels: the kernels
# @param config: the remora configuration 
# @return calculated lambda function values
#
findLambda <- function(training_matrix, kernels, config) {
    result <- list()

    # In Remora paper fn 4 = 7, thus it has been eliminated
    if (length(grep (FUNCTION_REMORA_EUCLIDEAN, config@weighting_function)) > 0) {
      # 0, Eclidean
      use_function <- 0
    } else if (length(grep (FUNCTION_REMORA_ONE_MINUS, config@weighting_function)) > 0){
      # 1, One Minus
      use_function <- 1
    } else if (length(grep (FUNCTION_REMORA_ONE_MINUS, config@weighting_function)) > 0){
      # 2, One Minus Squared
      use_function <- 2
    } else if (length(grep (FUNCTION_REMORA_MAHALANOBIS, config@weighting_function)) > 0){
      # 3, Mahalanobis (Gaussian)
      use_function <- 3
    } else if (length(grep (FUNCTION_REMORA_EXP_ONE_MINUS, config@weighting_function)) > 0){
      # 4, Exp One Minus
      use_function <- 4
    } else if (length(grep (FUNCTION_REMORA_EXP_ONE_MINUS_SQ, config@weighting_function)) > 0){
      # 5, Exp One Minus Squared
      use_function <- 5
    } else if (length(grep (FUNCTION_REMORA_EXP_ONE_LOG, config@weighting_function)) > 0){
      # 6, Exp One Log
      use_function <- 6
    } else if (length(grep (FUNCTION_REMORA_NORMALIZED_DIFFERENCE, config@weighting_function)) > 0){
      # 7, Normalized Difference
      use_function <- 7
    } else if (length(grep (FUNCTION_REMORA_NORMALIZED_DIFFERENCE_SQ, config@weighting_function)) > 0){
      # 8, Normalized Difference Squared
      use_function <- 8
    } else {
      stop("Undefined weighting function '", config@weighting_function ,"'.", , call. = TRUE)
    }
  
    # iterates for clusters/kernels
    for (cluster_name in names(kernels)) {      
      k <- kernels[[cluster_name]]
      current_cluster <- list()

      for (cluster_id in c(1:length(k@eigen_values))) {
        p <- k@eigen_vector[[cluster_id]]
        lambda <- k@eigen_values[[cluster_id]]
        pT <- t(p)
                       
        # apply selected function
        # in Remora paper fn 4 = 7, thus it has been eliminated
        if (use_function == 0) {
          fn_lambda <- diag(1, ncol(p), ncol(p))
        } else if (use_function == 1) {
          fn_lambda <- diag(1 - lambda, ncol(p), ncol(p))
        } else if (use_function == 2) {
          fn_lambda <- diag((1 - lambda)^2, ncol(p), ncol(p))          
        } else if (use_function == 3) {
          fn_lambda <- diag(1 / (lambda + config@epsilon), ncol(p), ncol(p))
        } else if (use_function == 4) {
          fn_lambda <- diag(exp(1 - lambda), ncol(p), ncol(p))          
        } else if (use_function == 5) {
          fn_lambda <- diag((exp(1 - lambda))^2, ncol(p), ncol(p))          
        } else if (use_function == 6) {
          fn_lambda <- diag((1 - log(lambda + config@epsilon)), ncol(p), ncol(p))          
        } else if (use_function == 7) {
          fn_lambda <- diag((1 - lambda) / (1 + lambda), ncol(p), ncol(p))          
        } else if (use_function == 8) {
          fn_lambda <- diag(((1 - lambda) / (1 + lambda))^2, ncol(p), ncol(p))          
        }
        lambda_matrix <- p %*% fn_lambda %*% pT        
        current_cluster[cluster_id] <- list(lambda_matrix)
      }
      result[cluster_name] <- list(current_cluster)
  }   
  
  result     
}


#
# Builds the Remora model.
# Recieves the training matrix and a configuration and returns a model.
#
# @param training_matrix: data matrix for training
# @param config: Remora algorithm configuration
# @return remora model
#
remora.model <- function(training_matrix, config) {
  kernels <- train(training_matrix, config)
  model_lambda <- findLambda(training_matrix, kernels, config)  
  s_values <- findS(training_matrix, kernels, model_lambda, config)
  i <- 0

  model_result <- data.frame()
    
  # iterates for clusters/kernels
  for (cluster_name in names(kernels)) {

    # get cluster specific information
    k <- kernels[[cluster_name]]
    clusters <- k@centroids
    for (cluster_index in c(1:nrow(clusters))) {
      # exatract cluster specific information      
      centroid <- clusters[cluster_index, ]
      kernel_s_value <- as.numeric(s_values[getFlatIndex(cluster_name, cluster_index)])

      i <- i + 1

      model_result[i, 1] <- cluster_name      
      model_result[i, 2] <- k@size[cluster_index]
      model_result[i, 3] <- s_values[getFlatIndex(cluster_name, cluster_index)] #kernel_s_value
      model_result[i, 4] <- cluster_index
      for (centroid_columns in c(1:length(centroid))) {  
        model_result[i, 4 + centroid_columns] <- centroid[centroid_columns]
      }
     }
  }  
  
  names(model_result) <- c("class", "size", "s", "centroid", names(centroid))
  
  result = remoraModel(config, model_result, model_lambda, kernels)
  
  if (verbose.show(config@verbose)) {
    cat('\nModel phase concluded.\n')
  }
  
  result
}

#
# Predicts a given classification using a specific model.
#
# @param model: Remora model
# @param data_matrix: data matrix to classify
# @return: list containing Remora classification prediction
#
# @see remora.model function
#
remora.predict <- function(model, data_matrix) {
  prediction_table <- classificationDistanceTable(model, data_matrix)
  prediction <- buildClassification(prediction_table, model@config)  
  
  rbf_result <- c(1:length(prediction))

  for(idx in c(1:length(prediction))) {
    rbf_result[idx] <- getFlatIndexClass(prediction[[idx]])
  }
  
  rbf_result
}
