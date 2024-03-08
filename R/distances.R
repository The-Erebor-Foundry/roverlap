
distance_point_to_line <- function(x_point, y_point, x_line, y_line) {
  check_line_coords__(x_line, y_line)
  distance_point_to_line_cpp(x_point, y_point, x_line, y_line)
}

is_na_null__ <- function(x) {
  is.null(x) || is.na(x)
}




point_overlap_line <- function(x_point, y_point, distance, x_line, y_line) {
  if (is_na_null__(x_point) || is_na_null__(y_point) || is_na_null__(distance)) {
    stop("[ERROR]: The `x_point`, `y_point`, and the `distance` arguments should neither be NA or NULL values!")
  }
  if (length(x_point) != 1L || length(y_point) != 1L) {
    stop("[ERROR]: The length of `x_point` and `y_point` arguments are greater than 1L. Did you mean to use the function `roverlap::npoints_overlapping_line()` instead?")
  }
  check_line_coords__(x_line, y_line)
  
  point_overlap_line_cpp(x_point, y_point, distance, x_line, y_line)
}




npoints_overlap_line <- function(x_points, y_points, distance, x_line, y_line) {
  if (any(is.na(x_points)) || any(is.na(y_points))) {
    stop("[ERROR]: The `x_points`, `y_points` arguments cannot contain NA values!")
  }
  if (length(x_points) == 0 || length(y_points) == 0) {
    stop("[ERROR]: Looks like the `x_point` and `y_point` arguments were filled with empty vectors.")
  }
  check_line_coords__(x_line, y_line)
  
  npoints_overlap_line_cpp(x_points, y_points, distance, x_line, y_line)
}




npoints_in_polygon <- function(x_points, y_points, x_polygon, y_polygon) {
  max_x <- max(x_polygon)
  min_x <- min(x_polygon)
  max_y <- max(y_polygon)
  min_y <- min(y_polygon)
  npoints_in_polygon_cpp(
    x_points,
    y_points,
    x_polygon,
    y_polygon,
    max_x,
    min_x,
    max_y,
    min_y
  )
}


point_in_polygon <- function(x_point, y_point, x_polygon, y_polygon) {
  max_x <- max(x_polygon)
  min_x <- min(x_polygon)
  max_y <- max(y_polygon)
  min_y <- min(y_polygon)
  point_in_polygon_cpp(
    x_point,
    y_point,
    x_polygon,
    y_polygon,
    max_x,
    min_x,
    max_y,
    min_y
  )
}








check_line_coords__ <- function(x_line, y_line) {
  if (length(x_line) != 2) {
    stop("[ERROR]: The x coordinates of the line must be a numeric vector of length 2!")
  }
  if (length(y_line) != 2) {
    stop("[ERROR]: The y coordinates of the line must be a numeric vector of length 2!")
  }
}