# Elbow method function
# Credit to Jonas for original idea and Esben Eickhardt for R implementation
# https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve

.elbow_finder <- function(x_values, y_values, correct_decrease = FALSE) {
  if(correct_decrease){
    test_increase <- stats::lm(y_values ~ x_values)
    if (stats::coef(test_increase)[2] > 0) {
      y_values <- -y_values
    }
  }
  
  # Max values to create line
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x),
                       y = c(max_y_y, max_x_y))
  
  # Creating straight line between the max values
  fit <- stats::lm(max_df$y ~ max_df$x)
  
  # Distance from point to line
  distances <- c()
  for (i in seq_along(x_values)){
    distances <- c(
      distances,
      abs(stats::coef(fit)[2] * x_values[i] - y_values[i] +
            stats::coef(fit)[1]) / sqrt(stats::coef(fit)[2]^2 + 1^2))
  }
  
  # Max distance point
  x_max_dist <- x_values[which.max(distances)]
  y_max_dist <- y_values[which.max(distances)]
  
  if(correct_decrease){
    if (stats::coef(test_increase)[2] > 0) {
      y_max_dist <- -y_max_dist
    }
  }
  
  return(c(x_max_dist, y_max_dist))
}
