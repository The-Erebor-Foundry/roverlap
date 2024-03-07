
#include <Rcpp.h>


// Min and Max functions between two double values.
double dmin(double v1, double v2);
double dmax(double v1, double v2);

// Distance between two points.
double distance_(double x1, double y1, double x2, double y2);
double line_slope_(double x1, double y1, double x2, double y2);
double y_intercept_(double x, double y, double slope);



bool _impl_overlap_vertical_line(double x_point,
                                 double distance,
                                 Rcpp::NumericVector x_line);
   
bool _impl_overlap_horizontal_line(double y_point,
                                 double distance,
                                 Rcpp::NumericVector y_line);

bool impl_ray_casting(double x_point,
                      double y_point,
                      Rcpp::NumericVector x_polygon,
                      Rcpp::NumericVector y_polygon,
                      double max_x_polygon,
                      double min_x_polygon,
                      double max_y_polygon,
                      double min_y_polygon);
                                        

                                        
                                        
bool _impl_overlap_line(double x_point,
                               double y_point,
                               double distance,
                               Rcpp::NumericVector x_line,
                               Rcpp::NumericVector y_line);