#include <math.h>

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP distance_point_to_line_impl(double x_point,
                                 double y_point,
                                 NumericVector x_line,
                                 NumericVector y_line) {

    double x1_line = x_line[0];
    double x2_line = x_line[1];
    double y1_line = y_line[0];
    double y2_line = y_line[1];
    
    double a = (y2_line - y1_line);
    double b = (x2_line - x1_line);
    double c = x1_line * y2_line - x2_line * y1_line;
    
    double d1 = abs(a * x_point + b * y_point + c);
    double d2 = sqrt(a * a + b * b);
    double distance = d1 / d2;
    return wrap(distance);
}


double distance_(double x1, double y1, double x2, double y2) {
  double s1 = pow(x2 - x1, 2.0);
  double s2 = pow(y2 - y1, 2.0);
  return sqrt(s1 + s2);
}

bool __point_overlap_line_impl(double x_point,
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
    double b = (x2_line - x1_line);
    double c = x1_line * y2_line - x2_line * y1_line;
    
    double d1 = abs(a * x_point + b * y_point + c);
    double d2 = sqrt(a * a + b * b);
    double d = d1 / d2;
    if (d < distance) {
        return true;
    }
    return false;
}


// [[Rcpp::export]]
SEXP point_overlap_line_impl(double x_point,
                             double y_point,
                             double distance,
                             NumericVector x_line,
                             NumericVector y_line) {
    bool result = __point_overlap_line_impl(
      x_point,
      y_point,
      distance,
      x_line,
      y_line
    );
    return wrap(result);
}


// [[Rcpp::export]]
SEXP npoints_overlapping_line_impl(NumericVector x_points,
                                   NumericVector y_points,
                                   double distance,
                                   NumericVector x_line,
                                   NumericVector y_line) {

  int n_points = x_points.length();
  LogicalVector results(n_points);
  for (int i = 0; i < n_points; i++) {
    results[i] = __point_overlap_line_impl(
      x_points[i],
      y_points[i],
      distance,
      x_line,
      y_line
    );
  }
  
  return wrap(results);
}
