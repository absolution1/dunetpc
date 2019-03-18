// ParabolicInterpolator.h
//
// David Adams
// December 2018
//
// Class to perform parabolic, i.e. fixed curvature, interpolation
// between two points.
//
//  y = a*(x-x1) + b*(x-x2) + c*(x-x1)*(x-x2)

#ifndef ParabolicInterpolator_H
#define ParabolicInterpolator_H

#include <string>

class TF1;

class ParabolicInterpolator {

public:

  struct Point {
    double x;
    double y;
    Point() : x(0.0), y(0.0) { };
    Point(double a_x, double a_y) : x(a_x), y(a_y) { }
  };

  using Name = std::string;

  // Ctors.
  ParabolicInterpolator(const Point xy1, const Point& xy2, double c);

  // Return the params.
  const Point& xy1() const { return m_xy1; }
  const Point& xy2() const { return m_xy2; }
  double a() const { return m_a; }
  double b() const { return m_b; }
  double c() const { return m_c; }

  // Return Root TF1.
  // If x2 > x1, then (x1, x2) is the range for the function.
  TF1* getTF1(Name ="parint", double x1=0.0, double x2=0.0) const;

private:

  Point m_xy1;
  Point m_xy2;
  double m_a;
  double m_b;
  double m_c;

};

#endif
