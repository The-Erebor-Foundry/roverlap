#include <math.h>

#include <Rcpp.h>
using namespace Rcpp;


#include "line.h"


// Point in polygon algorithm, using the ray-casting method.
bool _impl_ray_casting(double x_point,
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
  Line ray = Line(x_point, x_ray_end, y_point);
  Line side = Line(
    x_polygon[0],
    y_polygon[0],
    x_polygon[n_polygon_points - 1],
    y_polygon[n_polygon_points - 1]
  );
  if (side.ray_intersect_(ray)) {
    n_intersects += 1;
  }

  for (int i = 0; i < n_polygon_points - 1; i++) {
    double ax = x_polygon[i];
    double ay = y_polygon[i];
    double bx = x_polygon[i + 1];
    double by = y_polygon[i + 1];

    side = Line(ax, ay, bx, by);
    if (side.ray_intersect_(ray)) {
      n_intersects += 1;
    }
  }

  return n_intersects % 2 != 0;
}

