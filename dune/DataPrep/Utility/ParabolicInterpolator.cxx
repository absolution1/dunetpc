// ParabolicInterpolator.cxx

#include "ParabolicInterpolator.h"
#include "TF1.h"

//**********************************************************************

ParabolicInterpolator::
ParabolicInterpolator(const Point a_xy1, const Point& a_xy2, double c)
: m_xy1(a_xy1), m_xy2(a_xy2), m_a(0.0), m_b(0.0), m_c(c) {
  double dx = xy2().x - xy1().x;
  if ( fabs(dx) < 1.e-20 ) return;
  m_a = xy2().y/dx;
  m_b = xy1().y/dx;
}

//**********************************************************************

TF1* ParabolicInterpolator::getTF1(Name fnam, double x1, double x2) const {
  const Name sfun = "[2]*(x-[0]) + [3]*(x-[1]) + [4]*(x-[0])*(x-[1])";
  TF1* pf = new TF1(fnam.c_str(), sfun.c_str(), x1, x2);
  pf->SetParName(0, "x1");
  pf->SetParName(1, "x2");
  pf->SetParName(2, "a");
  pf->SetParName(3, "b");
  pf->SetParName(4, "c");
  pf->SetParameter(0, xy1().x);
  pf->SetParameter(1, xy2().x);
  pf->SetParameter(2, a());
  pf->SetParameter(3, b());
  pf->SetParameter(4, c());
  return pf;
}

//**********************************************************************
