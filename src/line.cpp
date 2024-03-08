#include <math.h>

#include <Rcpp.h>
using namespace Rcpp;

#include "line.h"
#include "utils.h"
#define DOUBLE_MAX 184467440737095516.0;


Line::Line(double x1, double y1, double x2, double y2) {
  if (x2 == x1) {
    is_vertical = true;
  }
  lower_y = dmin(y1, y2);
  upper_y = dmax(y1, y2);
  lower_x = dmin(x1, x2);
  upper_x = dmax(x1, x2);
  slope = line_slope_(x1, y1, x2, y2);
  y_intercept = y_intercept_(x1, y1, slope);
}

Line::Line(double x1, double x2, double y) {
  is_vertical = false;
  lower_y = y;
  upper_y = y;
  lower_x = dmin(x1, x2);
  upper_x = dmax(x1, x2);
  slope = 0;
  y_intercept = y;
}

bool Line::ray_intersect_(Line ray) {
  double eps = 0.00001;
  if (ray.lower_x == lower_x || ray.lower_x == upper_x) {
    ray.lower_x += eps;  
  }
  if (ray.lower_y == lower_y || ray.lower_y == upper_y) {
    ray.lower_y += eps;
  }

  if (is_vertical) {
    if ((ray.upper_y <= upper_y)
      && (ray.lower_y >= lower_y)
      && (ray.lower_x < lower_x)) {
      return true;
    }

    return false;
  }

  if (lower_y == upper_y) {
    if (ray.lower_y != lower_y) {
      return true;
    }

    return false;
  }

  if (slope == ray.slope) {
    // They are parallel
    return false;
  }
  // Point of intersection:
  double intersect_x = (ray.lower_y - y_intercept) / slope;
  double intersect_y = (slope * intersect_x + y_intercept);


  if ((intersect_x >= lower_x && intersect_x <= upper_x)
    && (intersect_y >= lower_y && intersect_y <= upper_y)) {

    return true;
  }

  return false;
}


double Line::line_slope_(double x1, double y1, double x2, double y2) {
  if (is_vertical) {
    return DOUBLE_MAX;
  }
  double dy = abs(upper_y - lower_y);
  double dx = abs(upper_x - lower_x);
  return dy / dx;
}

double Line::y_intercept_(double x, double y, double slope) {
  if (is_vertical) {
    return 0;
  }
  return y - x * slope;
}













// Overlapping functions for infinite lines ====================================

bool _impl_overlap_line(double x_point,
                        double y_point,
                        double distance,
                        NumericVector x_line,
                        NumericVector y_line) {

  double x1_line = x_line[0];
  double x2_line = x_line[1];
  double y1_line = y_line[0];
  double y2_line = y_line[1];
  double dist_to1 = distance_(x_point, y_point, x1_line, y1_line);
  if (dist_to1 < distance) {
    return true;
  }

  double dist_to2 = distance_(x_point, y_point, x2_line, y2_line);
  if (dist_to2 < distance) {
    return true;
  }

  double a = (y2_line - y1_line);
  double b = (x1_line - x2_line);
  double c = y1_line * (x2_line - x1_line) - (y2_line - y1_line) * x1_line;
  double d1 = abs(a * x_point + b * y_point + c);
  double d2 = sqrt(a * a + b * b);
  double d = d1 / d2;
  if (d < distance) {
    return true;
  }
  return false;
}


bool _impl_overlap_vertical_line(double x_point,
                                 double distance,
                                 NumericVector x_line) {

  double x1_line = x_line[0];
  double d_x = abs(x1_line - x_point);
  return wrap(d_x < distance);
}


bool _impl_overlap_horizontal_line(double y_point,
                                   double distance,
                                   NumericVector y_line) {
  double y1_line = y_line[0];
  double d_y = abs(y1_line - y_point);
  return wrap(d_y < distance);
}

