#' entropy_weight
#' @description Calculate the entropy weight for each column in a data frame or matrix.
#' @param df A data table or matrix with numeric values.
#' @return A numeric vector of entropy weights for each column.
#' @export
entropy_weight <- function(df) {
  # Ensure this is a data frame or matrix
  if (!is.data.frame(df) && !is.matrix(df)) {
    stop("Input must be a data frame or matrix")
  }
  df <- as.matrix(df)
  # Call C++ function
  entropy_weight_cpp(df)
}
