// cePulserToAdc.cxx
//
// David Adams
// January 2019
//
// ADC response to CE pulser.ProtoDUNE CE ADC pulser charge. See the circuit diagram
//   http://internal.dunescience.org/people/dladams/docs/pulserCircuit.pdf

#ifndef cePulserToAdc_H
#define cePulserToAdc_H

#include "cePulser.h"

//**********************************************************************

// Function that returns the ADC height or area for a pulser gain setting.
//   ia - gain setting (minus inverts charge)
//   adcScale - d(ADC count)/dQin
//   rLines - FPGA IO resistances (high to low) [kOhm]
//   r6 - Fixed resistance to Vrail.
//   cinj - Injection capacitance (fF)
// Returns charge in fC.
// Use cinj = 1 to get voltage in mV.
inline
double cePulserToAdc(int ia, double adcScale, double pedestal, double negScale,
                     const double rfLines[6], double r6, double qvscale) {
  double dadc = adcScale*cePulser(ia, rfLines, r6, qvscale);
  if ( ia < 0.0 ) dadc *= negScale;
  return dadc + pedestal;
}

double cePulserToAdc(double* x, double* pars) {
  int ia = x[0];
  double adcScale = pars[0];
  double pedestal = pars[1];
  double negScale = pars[2];
  const double* rfLines = pars + 3;
  double r6 = pars[9];
  double qvscale = pars[10];
  return cePulserToAdc(ia, adcScale, pedestal, negScale, rfLines, r6, qvscale);
}

//**********************************************************************

// Generate a TF1 for the ADC response.
inline
TF1* cePulserToAdcTF1(std::string fname="adcpulser", double tol =0.01,
                      double adcScale = 3.0*14.0) {
  using Index = unsigned int;
  const Index npar = 11;
  TF1* pf = new TF1(fname.c_str(), cePulserToAdc, -64, 64, npar);
  initCePulserTF1(pf, 3, tol);
  pf->SetParName(0, "AdcScale");
  pf->SetParName(1, "Pedestal");
  pf->SetParName(2, "NegScale");
  pf->SetParameter(0, adcScale);
  pf->SetParameter(1, 0.0);
  pf->SetParameter(2, 1.0);
  pf->FixParameter(10, 1.8*183.0);
  return pf;
}
 
//**********************************************************************

#endif
