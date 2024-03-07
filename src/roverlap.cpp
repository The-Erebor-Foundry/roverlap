#include <math.h>

#include <Rcpp.h>
using namespace Rcpp;

#include "roverlap.h"


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
  bool result = impl_ray_casting(
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











// Raycasting algorithm

bool impl_ray_casting(double x_point,
                      double y_point,
                      NumericVector x_polygon,
                      NumericVector y_polygon,
                      double max_x_polygon,
                      double min_x_polygon,
                      double max_y_polygon,
                      double min_y_polygon) {
  
  double eps = 0.000001;
  int n_intersects = 0;
  int n_polygon_points = x_polygon.length();
  
  double x_ray_end = max_x_polygon + 1e6;
  double ray_slope = 0;
  double ray_y_intercept = y_point;
  
  for (int i = 0; i < n_polygon_points - 1; i++) {
    double ax = x_polygon[i];
    double ay = y_polygon[i];
    double bx = x_polygon[i + 1];
    double by = y_polygon[i + 1];
    
    if (y_point == ay || y_point == by) {
      y_point += eps;
    }
    if (x_point == ax || x_point == bx) {
      x_point += eps;
    }
    
    bool is_vertical_side = ax == bx;
    if (is_vertical_side) {
      if (x_point < ax) {
        n_intersects += 1;
      }
      continue;
    }
    
    double side_slope = line_slope_(ax, ay, bx, by);
    double side_y_intercept = y_intercept_(ax, ay, side_slope);
    if (side_slope == ray_slope) {
      // They are parallel
      continue;
    }
    // Point of intersection:
    double intersect_x = (side_y_intercept - ray_y_intercept) / (ray_slope - side_slope);
    double intersect_y = (side_slope * intersect_x + side_y_intercept);
    
    double lower_y = dmin(ay, by);
    double upper_y = dmax(ay, by);
    double lower_x = dmin(ax, bx);
    double upper_x = dmax(ax, bx);
    if ((intersect_x >= lower_x && intersect_x <= upper_x)
          && (intersect_y >= lower_y && intersect_y <= upper_y)) {
      
      n_intersects += 1;
      continue;
    }
  }
  
  return n_intersects % 2 != 0;
}




bool __impl_ray_casting(double x_point,
                      double y_point,
                      NumericVector x_polygon,
                      NumericVector y_polygon) {
  
  double eps = 0.000001;
  int n_intersects = 0;
  int n_polygon_points = x_polygon.length();
  
  double ax = x_polygon[0];
  double ay = y_polygon[0];
  double bx = x_polygon[n_polygon_points];
  double by = y_polygon[n_polygon_points];
  
  double lower_y = dmin(ay, by);
  double upper_y = dmax(ay, by);
  double lower_x = dmin(ax, bx);
  double upper_x = dmax(ax, bx);
  if (y_point == ay || y_point == by) {
    y_point += eps;
  }
  if (x_point == ax || x_point == bx) {
    x_point += eps;
  } 
  
  if (x_point > lower_x) {
    // The ray intersects with the edge
    n_intersects += 1;
  }
  
  double m1_edge = ((by - ay) == 0) ? eps : (by - ay);
  double m2_edge = ((bx - ax) == 0) ? eps : (bx - ax);
  double m_edge = m1_edge / m2_edge;
  
  double m_point = (y_point - ay) / (x_point - ax);
  if (m_point >= m_edge) {
    // The ray intersects with the edge
    n_intersects += 1;
  }
  
  
  for (int i = 0; i < n_polygon_points - 1; i++) {
    double ax = x_polygon[i];
    double ay = y_polygon[i];
    double bx = x_polygon[i + 1];
    double by = y_polygon[i + 1];
    
    double lower_y = dmin(ay, by);
    double upper_y = dmax(ay, by);
    double lower_x = dmin(ax, bx);
    double upper_x = dmax(ax, bx);
    if (y_point == ay || y_point == by) {
      y_point += eps;
    }
    if (x_point == ax || x_point == bx) {
      x_point += eps;
    }
    
    if (y_point > by || y_point < ay || x_point > upper_x) {
      // The horizontal ray does not intersect with the edge
      continue;
    }
    
    if (x_point > lower_x) {
      // The ray intersects with the edge
      n_intersects += 1;
      continue;
    }
    
    double m1_edge = ((by - ay) == 0) ? eps : (by - ay);
    double m2_edge = ((bx - ax) == 0) ? eps : (bx - ax);
    double m_edge = m1_edge / m2_edge;
    
    double m_point = (y_point - ay) / (x_point - ax);
    if (m_point >= m_edge) {
      // The ray intersets with the edge
      n_intersects += 1;
      continue;
    }
      
  }
  
  bool inside = n_intersects % 2 != 0;
  return inside;
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
















// Utility functions ==================================================

double line_slope_(double x1, double y1, double x2, double y2) {
  double dy = abs(y1 - y2);
  double dx = abs(x1 - x2);
  return dy / dx;
}

double y_intercept_(double x, double y, double slope) {
  return y - x * slope;
}

double distance_(double x1, double y1, double x2, double y2) {
  double s1 = pow(x2 - x1, 2.0);
  double s2 = pow(y2 - y1, 2.0);
  return sqrt(s1 + s2);
}


double dmax(double v1, double v2) {
  if (v1 > v2) {
    return v1;
  }
  return v2;
}

double dmin(double v1, double v2) {
  if (v1 < v2) {
    return v1;
  }
  return v2;
}