#include <math.h>

#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"
#include "line.h"
#include "polygon.h"


// Exported functions to R ================================================================

// [[Rcpp::export]]
SEXP distance_point_to_line_cpp(double x_point,
                                double y_point,
                                NumericVector x_line,
                                NumericVector y_line) {

  double x1_line = x_line[0];
  double x2_line = x_line[1];
  double y1_line = y_line[0];
  double y2_line = y_line[1];

  double a = (y2_line - y1_line);
  double b = (x1_line - x2_line);
  double c = y1_line * (x2_line - x1_line) - (y2_line - y1_line) * x1_line;
  double d1 = abs(a * x_point + b * y_point + c);
  double d2 = sqrt(a * a + b * b);
  double d = d1 / d2;
  return wrap(d);
}


// [[Rcpp::export]]
SEXP point_in_polygon_cpp(double x_point,
                          double y_point,
                          NumericVector x_polygon,
                          NumericVector y_polygon,
                          double max_x_polygon,
                          double min_x_polygon,
                          double max_y_polygon,
                          double min_y_polygon) {
  bool result = _impl_ray_casting(
    x_point,
    y_point,
    x_polygon,
    y_polygon,
    max_x_polygon,
    min_x_polygon,
    max_y_polygon,
    min_y_polygon
  );

  return wrap(result);
}


// [[Rcpp::export]]
SEXP npoints_in_polygon_cpp(NumericVector x_points,
                            NumericVector y_points,
                            NumericVector x_polygon,
                            NumericVector y_polygon,
                            double max_x_polygon,
                            double min_x_polygon,
                            double max_y_polygon,
                            double min_y_polygon) {
  int n_points = x_points.length();
  LogicalVector results(n_points);
  for (int i = 0; i < n_points; i++) {
    results[i] = _impl_ray_casting(
      x_points[i],
      y_points[i],
      x_polygon,
      y_polygon,
      max_x_polygon,
      min_x_polygon,
      max_y_polygon,
      min_y_polygon
    );
  }

  return wrap(results);
}


// [[Rcpp::export]]
SEXP point_overlap_line_cpp(double x_point,
                            double y_point,
                            double distance,
                            NumericVector x_line,
                            NumericVector y_line) {

  double x1_line = x_line[0];
  double x2_line = x_line[1];
  double y1_line = y_line[0];
  double y2_line = y_line[1];
  if (x1_line == x2_line) {
    // Vertical line
    bool result = _impl_overlap_vertical_line(
      x_point,
      distance,
      x_line
    );
    return wrap(result);
  }
  if (y1_line == y2_line) {
    // Horizontal line
    bool result = _impl_overlap_horizontal_line(
      y_point,
      distance,
      y_line
    );
    return wrap(result);
  }


  bool result = _impl_overlap_line(
    x_point,
    y_point,
    distance,
    x_line,
    y_line
  );
  return wrap(result);
}


// [[Rcpp::export]]
SEXP npoints_overlap_line_cpp(NumericVector x_points,
                              NumericVector y_points,
                              double distance,
                              NumericVector x_line,
                              NumericVector y_line) {

  double x1_line = x_line[0];
  double x2_line = x_line[1];
  double y1_line = y_line[0];
  double y2_line = y_line[1];

  int n_points = x_points.length();
  LogicalVector results(n_points);
  for (int i = 0; i < n_points; i++) {

    if (x1_line == x2_line) {
      results[i] = _impl_overlap_vertical_line(
        x_points[i],
        distance,
        x_line
      );
      continue;
    }

    if (y1_line == y2_line) {
      results[i] = _impl_overlap_horizontal_line(
        y_points[i],
        distance,
        y_line
      );
      continue; 
    }

    results[i] = _impl_overlap_line(
      x_points[i],
      y_points[i],
      distance,
      x_line,
      y_line
    );

  }

  return wrap(results);
}








