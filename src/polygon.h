#pragma once
#include <Rcpp.h>


bool _impl_ray_casting(double x_point,
                       double y_point,
                       Rcpp::NumericVector x_polygon,
                       Rcpp::NumericVector y_polygon,
                       double max_x_polygon,
                       double min_x_polygon,
                       double max_y_polygon,
                       double min_y_polygon);
