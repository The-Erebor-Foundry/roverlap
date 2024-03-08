#pragma once
#include <Rcpp.h>

class Line {
public:
  double lower_x;
  double lower_y;
  double upper_x;
  double upper_y;
  bool is_vertical;
  double slope;
  double y_intercept;

public:
  Line(double x1, double y1, double x2, double y2);
  Line(double x1, double x2, double y);
  bool ray_intersect_(Line ray);
  double line_slope_(double x1, double y1, double x2, double y2);
  double y_intercept_(double x, double y, double slope);
};



bool _impl_overlap_line(double x_point,
                        double y_point,
                        double distance,
                        Rcpp::NumericVector x_line,
                        Rcpp::NumericVector y_line);

bool _impl_overlap_vertical_line(double x_point,
                                 double distance,
                                 Rcpp::NumericVector x_line);

bool _impl_overlap_horizontal_line(double y_point,
                                   double distance,
                                   Rcpp::NumericVector y_line);

