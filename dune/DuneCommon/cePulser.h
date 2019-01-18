// cePulser.cxx
//
// David Adams
// January 2019
//
// ProtoDUNE CE ADC pulser charge. See the circuit diagram
//   http://internal.dunescience.org/people/dladams/docs/pulserCircuit.pdf

#ifndef cePulser_H
#define cePulser_H

#include "TF1.h"
#include <string>

//**********************************************************************

// Function that returns the voltage on the FE injection capacitor.
//   ia - gain setting (minus inverts charge)
//   rLines - FPGA IO resistances (high to low) [kOhm]
//   r6 - Fixed resistance to Vrail.
//   cinj - Injection capacitance (fF)
// Returns charge in fC.
// Use cinj = 1 to get voltage in mV.
inline
double cePulserVoltage(int iaIn, double rfLines[6], double r6, double vfpga) {
  using Index = unsigned int;
  if ( r6 < 0.0 ) return 0.0;
  if ( iaIn < 0 ) return 0.0;
  if ( iaIn > 63 ) return 0.0;
  Index ia = iaIn;
  double con1 = 0.0;                      // Conductance to V=0
  double con2 = r6 > 0.0 ? 1.0/r6 : 0.0;  // Conductance to V=Vfpga
  for ( Index ilin=0; ilin<6; ++ilin ) {
    double rfLine = rfLines[ilin];
    if ( rfLine <= 0.0 ) continue;    // Means this line is disconnected.
    double conLine = 1.0/rfLine;
    // If this bit is set then the path is to Vfpga.
    if ( (ia >> ilin) & 0x1) con2 += conLine;
    // Otherwise, the path is to ground.
    else con1 += conLine;
  }
  double fac = 1.0;
  if ( con1 > 0.0 ) {
    double rtot = 1.0/con1 + 1.0/con2;
    fac = 1.0/con1/rtot;
  }
  const double vrail = 1.8;
  return fac*vfpga;
}

//**********************************************************************

// Function that returns the voltage or charge on the FE injection capacitor
// offset to be zero for pulser gain setting ia = 0;
//   ia - gain setting: 0 - 63 (minus inverts charge)
//   rLines - FPGA IO resistances (high to low) [kOhm]
//   r6 - Fixed resistance to Vrail.
//   qvscale - Vfpga*Cinj - Injection capacitance (fF)
inline
double cePulser(int ia, double rfLines[6], double r6, double qvscale) {
  if ( ia == 0 ) return 0.0;
  if ( ia > 63 ) return 0.0;
  if ( ia < 0 ) return cePulser(-ia, rfLines, r6, qvscale);
  return cePulserVoltage(ia, rfLines, r6, qvscale) -
         cePulserVoltage( 0, rfLines, r6, qvscale);
}

//**********************************************************************

// Pulser functions with Root syntax.
inline
double cePulserVoltage(double* x, double* pars) {
  int ia = x[0];
  double* rfLines = pars;
  double r6 = pars[6];
  double qvscale = pars[7];
  return cePulserVoltage(ia, pars, r6, qvscale);
}

//**********************************************************************

inline
double cePulser(double* x, double* pars) {
  int ia = x[0];
  double* rfLines = pars;
  double r6 = pars[6];
  double qvscale = pars[7];
  return cePulser(ia, pars, r6, qvscale);
}

//**********************************************************************

// Generate a TF1 for the pulser Voltage function.
inline
TF1* cePulserVoltageTF1(std::string fname="cepulserVoltage", double tol =0.01) {
  using Index = unsigned int;
  const Index npar = 8;
  TF1* pf = new TF1(fname.c_str(), cePulserVoltage, 0, 63, npar);
  pf->SetParName(0, "R0");
  pf->SetParName(1, "R1");
  pf->SetParName(2, "R2");
  pf->SetParName(3, "R3");
  pf->SetParName(4, "R4");
  pf->SetParName(5, "R5");
  pf->SetParName(6, "R6");
  pf->SetParName(7, "Vfpga");
  pf->SetParameter(0, 32.40);
  pf->SetParameter(1, 16.02);
  pf->SetParameter(2,  8.06);
  pf->SetParameter(3,  4.02);
  pf->SetParameter(4,  2.00);
  pf->SetParameter(5,  1.00);
  pf->SetParameter(6,  1.00);
  pf->SetParameter(7, 1.8);
  if ( tol > 0.0 && tol < 1.0 ) {
    double fmin = 1 - tol;
    double fmax = 1 + tol;
    for ( Index ipar=0; ipar<npar; ++ipar ) {
      double val = pf->GetParameter(ipar);
      pf->SetParLimits(ipar, fmin*val, fmax*val);
    }
  }
  return pf;
}

//**********************************************************************

// Generate a TF1 for the pulser function.
inline
TF1* cePulserTF1(std::string fname="cepulser", double tol =0.01) {
  using Index = unsigned int;
  const Index npar = 8;
  TF1* pf = new TF1(fname.c_str(), cePulser, -64, 64, npar);
  pf->SetParName(0, "R0");
  pf->SetParName(1, "R1");
  pf->SetParName(2, "R2");
  pf->SetParName(3, "R3");
  pf->SetParName(4, "R4");
  pf->SetParName(5, "R5");
  pf->SetParName(6, "R6");
  pf->SetParName(7, "QVscale");
  pf->SetParameter(0, 32.40);
  pf->SetParameter(1, 16.02);
  pf->SetParameter(2,  8.06);
  pf->SetParameter(3,  4.02);
  pf->SetParameter(4,  2.00);
  pf->SetParameter(5,  1.00);
  pf->SetParameter(6,  1.00);
  pf->SetParameter(7, 1.8*183.0);
  if ( tol > 0.0 && tol < 1.0 ) {
    double fmin = 1 - tol;
    double fmax = 1 + tol;
    for ( Index ipar=0; ipar<npar; ++ipar ) {
      double val = pf->GetParameter(ipar);
      pf->SetParLimits(ipar, fmin*val, fmax*val);
    }
  }
  return pf;
}
 
//**********************************************************************

#endif
