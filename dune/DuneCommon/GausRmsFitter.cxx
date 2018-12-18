// GausRmsFitter.cxx

#include "GausRmsFitter.h"
#include "dune/DuneCommon/gausTF1.h"
#include "dune/DuneCommon/quietHistFit.h"
#include "TF1.h"
#include "TH1.h"
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;
using std::setw;

//**********************************************************************

GausRmsFitter::GausRmsFitter(Name fnam)
: m_sigma0(0.0), m_nsigma(0.0), m_fnam(fnam), m_LogLevel(0) { }

//**********************************************************************

GausRmsFitter::GausRmsFitter(double sigma0, double nsigma, Name fnam)
: m_sigma0(sigma0), m_nsigma(nsigma), m_fnam(fnam), m_LogLevel(1) { }

//**********************************************************************

int GausRmsFitter::fit(TH1* ph, double mean0) const {
  Name myname = "GausRmsFitter::fit: ";
  if ( ph->GetEntries() == 0 ) {
    if ( m_LogLevel >= 2 ) cout << "Histogram has no entries." << endl;
    return 101;  // E.g. empty histo.
  }
  if ( m_LogLevel >= 2 ) cout << myname << "Fitting with sigma " << m_sigma0
                              << " and nsigma " << m_nsigma << endl;
  double mean = mean0;
  double sigma = ph->GetRMS();
  double height = ph->GetMaximum();
  int ibin1 = ph->GetXaxis()->GetFirst();
  int ibin2 = ph->GetXaxis()->GetLast();
  double xmin = ph->GetXaxis()->GetBinLowEdge(ibin1);
  double xmax = ph->GetXaxis()->GetBinUpEdge(ibin2);
  double x1 = xmin;
  double x2 = xmax;
  // Evaluate mean and sigma.
  if ( m_sigma0 > 0.0 && m_nsigma > 0.0 ) {
    sigma = 0.0;
    double sigmaOld = m_sigma0;
    if ( m_LogLevel >= 2 ) {
      cout << myname << "Start loop with mean, sigma: " << mean << ", " << sigmaOld << endl;
    }
    const Index maxtry = 100;
    Index itry = 100;
    for ( Index itry=0; itry<maxtry; ++itry ) {
      double dx = m_nsigma*sigmaOld;
      double x1 = mean - dx;
      if ( x1 < xmin ) x1 = xmin;
      double x2 = mean + dx;
      if ( x1 > xmax ) x1 = xmax;
      ph->GetXaxis()->SetRangeUser(x1, x2);
      double sigmaNew = ph->GetRMS();
      mean = ph->GetMean();
      if ( m_LogLevel >= 2 ) {
        cout << myname << setw(5) << itry << ": " << x1 << ", " << x2 << "): "
             << mean << ", " << sigmaNew << endl;
      }
      if ( sigmaNew < 1.001*sigmaOld ) {
        sigma = sigmaNew;
        break;
      }
      sigmaOld = sigmaNew;
    }
    if ( itry >= 100 ) cout << myname << "WARNING: too many iterations." << endl;
  // Use configured sigma.
  } else if ( m_sigma0 > 0.0 ) {
    sigma = m_sigma0;
  } else {
    cout << myname << "Invalid configuration." << endl;
    return 102;
  }
  ph->GetXaxis()->SetRange(ibin1, ibin2);
  if ( sigma <= m_sigma0 && m_sigma0 > 0.0 ) {
    cout << myname << "Unexpected value for sigma: " << sigma << endl;
    return 103;
  }
  TF1* pffix = gausTF1(height, mean, sigma, m_fnam);
  pffix->FixParameter(1, mean);
  pffix->FixParameter(2, sigma);
  pffix->SetRange(x1, x2);
  string fopt = "S";
  if ( m_LogLevel <= 1 ) fopt += "Q";
  int fstat = quietHistFit(ph, pffix, fopt);
  if ( m_LogLevel >= 1 ) cout << myname << "  status " << fstat << endl;
  if ( fstat == 0 ) {
    if ( m_LogLevel >= 1 ) cout << myname << "Fit suceeded." << endl;
    return 0;
  }
  // Otherwise we may get a crash when we try to view saved copy of histo:
  ph->GetListOfFunctions()->Clear();
  delete pffix;
  if ( m_LogLevel >= 1 ) cout << myname << "Fit failed." << endl;
  return 1;
}

//**********************************************************************
