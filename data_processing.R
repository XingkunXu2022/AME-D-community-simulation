get_dif_matrix <- function(values){
  
  output_mat <- matrix(nrow = length(values),
                       ncol = length(values))
  
  for (i in 1:length(values)){
    for (j in 1:length(values)){
      output_mat[i,j] <- abs(values[i]-values[j])
    }
  }
  
  return(output_mat)
}


match_matrix_values <- function(matrix1, matrix2){
  
  matrix1_values <- c()
  matrix2_values <- c()
  
  if( nrow(matrix1) != nrow(matrix2) | ncol(matrix1) != ncol(matrix2)){
    print('Matrix1 and matrix2 are not matched!')
    return(NULL)
  }
  
  for (i in 1:nrow(matrix1)){
    for (j in 1:ncol(matrix2)){
      matrix1_values <- c(matrix1_values,matrix1[i,j])
      matrix2_values <- c(matrix2_values,matrix2[i,j])
    }
  }
  
  output_data_frame <- data.frame(matrix1 = matrix1_values,
                                  matrix2 = matrix2_values)
  
  return(output_data_frame)
}