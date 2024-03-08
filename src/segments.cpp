#include <math.h>

#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"

// Overlapping functions for line segments ====================================

bool _impl_overlap_vertical_segment(double x_point,
                                    double y_point,
                                    double distance,
                                    NumericVector x_line,
                                    NumericVector y_line) {
  double x1_line = x_line[0];
  double x2_line = x_line[1];
  double y1_line = y_line[0];
  double y2_line = y_line[1];
  double d_x = abs(x1_line - x_point);
  if (d_x < distance) {
    double y_upper_line = dmax(y1_line, y2_line);
    double y_lower_line = dmin(y1_line, y2_line);
    if (y_point >= y_lower_line && y_point <= y_upper_line) {
      return true;
    }

    double d_upper = distance_(x_point, y_point, x1_line, y_upper_line);
    double d_lower = distance_(x_point, y_point, x1_line, y_lower_line);
    if (d_upper < distance || d_lower < distance) {
      return true;
    }

  }

  return false;
}

