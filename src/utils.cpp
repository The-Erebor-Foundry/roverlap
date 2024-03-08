#include <math.h>

#define DOUBLE_MAX 184467440737095516.0;

// Distance between two points
double distance_(double x1, double y1, double x2, double y2) {
  double s1 = pow(x2 - x1, 2.0);
  double s2 = pow(y2 - y1, 2.0);
  return sqrt(s1 + s2);
}

// Get max double
double dmax(double v1, double v2) {
  if (v1 > v2) {
    return v1;
  }
  return v2;
}


// Get min double
double dmin(double v1, double v2) {
  if (v1 < v2) {
    return v1;
  }
  return v2;
}


double line_slope_(double x1, double y1, double x2, double y2) {
  if (x1 == x2) {
    // Line is vertical, so its slope is infinity.
    return DOUBLE_MAX;
  }
  double dy = abs(y1 - y2);
  double dx = abs(x1 - x2);
  return dy / dx;
}

double y_intercept_(double x, double y, double slope) {
  return y - x * slope;
}