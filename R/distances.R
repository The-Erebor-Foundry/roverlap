
distance_point_to_line <- function(x_point, y_point, x_line, y_line) {
  check_line_coords__(x_line, y_line)
  distance_point_to_line_impl(x_point, y_point, x_line, y_line)
}


point_overlapping_line <- function(x_point, y_point, distance, x_line, y_line) {
  check_line_coords__(x_line, y_line)
  point_overlap_line_impl(x_point, y_point, distance, x_line, y_line)
}

check_line_coords__ <- function(x_line, y_line) {
  if (length(x_line) != 2) {
    stop("[ERROR]: The x coordinates of the line must be a numeric vector of length 2!")
  }
  if (length(y_line) != 2) {
    stop("[ERROR]: The y coordinates of the line must be a numeric vector of length 2!")
  }
}