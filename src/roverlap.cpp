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
