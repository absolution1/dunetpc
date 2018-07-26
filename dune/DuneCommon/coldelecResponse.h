// coldelecReponse.cxx
//
// David Adams
// May 2018
// -----
//  From Brett Viren May 2018
//  Cold Electronics response function.
//  How was this function derived?
//    1. Hucheng provided a transfer function of our electronics, which
//       is the Laplace transformation of the shaping function.
//    2. We did a anti-Laplace inverse of the shaping function
//    3. The one you saw in the code is basically the result of the inversion
//  - time is time in system of units
//  - gain in units of volts/charge gives peak value of the response to
//    a delta function of current integrating to unit charge.
//  - shaping is the shaping time in system of units
//  - the hard-coded numbers are the result of the inverting the
//    Lapalace transformation in Mathematica.
// -----

#ifndef coldelecReponse_H
#define coldelecReponse_H

#include "TF1.h"
#include <string>

inline
double coldelecResponse(double time, double gain, double shaping) {
    if (time <=0 || time >= 10*shaping) { // range of validity
	return 0.0;
    }

    const double reltime = time/shaping;

    // a scaling is needed to make the anti-Lapalace peak match the expected gain
    // fixme: this scaling is slightly dependent on shaping time.  See response.py
    gain *= 10*1.012;

    return 4.31054*exp(-2.94809*reltime)*gain
	-2.6202*exp(-2.82833*reltime)*cos(1.19361*reltime)*gain
	-2.6202*exp(-2.82833*reltime)*cos(1.19361*reltime)*cos(2.38722*reltime)*gain
	+0.464924*exp(-2.40318*reltime)*cos(2.5928*reltime)*gain
	+0.464924*exp(-2.40318*reltime)*cos(2.5928*reltime)*cos(5.18561*reltime)*gain
	+0.762456*exp(-2.82833*reltime)*sin(1.19361*reltime)*gain
	-0.762456*exp(-2.82833*reltime)*cos(2.38722*reltime)*sin(1.19361*reltime)*gain
	+0.762456*exp(-2.82833*reltime)*cos(1.19361*reltime)*sin(2.38722*reltime)*gain
 	-2.620200*exp(-2.82833*reltime)*sin(1.19361*reltime)*sin(2.38722*reltime)*gain 
	-0.327684*exp(-2.40318*reltime)*sin(2.5928*reltime)*gain + 
	+0.327684*exp(-2.40318*reltime)*cos(5.18561*reltime)*sin(2.5928*reltime)*gain
	-0.327684*exp(-2.40318*reltime)*cos(2.5928*reltime)*sin(5.18561*reltime)*gain
	+0.464924*exp(-2.40318*reltime)*sin(2.5928*reltime)*sin(5.18561*reltime)*gain;
}

inline
double coldelecResponseFunction(double* x, double* pars) {
  double gain = pars[0];
  double time = x[0];
  double shaping = pars[1];
  double offset = pars[2];
  return coldelecResponse(time-offset, gain, shaping);
}

inline
TF1* coldelecResponseTF1(double gainIn, double shapingIn, double t0, std::string fname ="ceresp") {
  double gain = gainIn != 0.0 ? gainIn : 10.0;
  double shaping = shapingIn > 0 ? shapingIn : 1.0;
  bool havePars = gainIn > 0 && shaping > 0.0;
  //TF1* pf = new TF1(fname.c_str(), coldelecResponseFunction, t0, t0+10.0, 3, 1, false);
  TF1* pf = new TF1(fname.c_str(), coldelecResponseFunction, t0, t0+10.0*shaping, 3);
  pf->SetParName(0, "Gain");
  pf->SetParName(1, "Shaping");
  pf->SetParName(2, "T0");
  pf->SetParameter(0, gain);
  pf->SetParameter(1, shaping);
  pf->SetParameter(2, t0);
  double limfac = havePars ? 2.0 : 10.0;
  if ( gain > 0.0 ) {
    pf->SetParLimits(0, gain/limfac, limfac*gain);
  } else {
    pf->SetParLimits(0, limfac*gain, gain/limfac);
  }
  pf->SetParLimits(1, shaping/limfac, limfac*shaping);
  return pf;
}
 
#endif
