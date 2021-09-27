// AdcCodeMitigator_tool.cc

#include "AdcCodeMitigator.h"
#include <iostream>
#include <sstream>

using std::string;
using std::cout;
using std::endl;
using std::istringstream;
using Index = AdcCodeMitigator::Index;
using IndexVector = AdcCodeMitigator::IndexVector;
using Name = AdcCodeMitigator::Name;
using NameVector = std::vector<Name>;

//**********************************************************************
// Class methods.
//**********************************************************************

AdcCodeMitigator::AdcCodeMitigator(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_FixFlags(ps.get<IndexVector>("FixFlags")),
  m_InterpolateFlags(ps.get<IndexVector>("InterpolateFlags")),
  m_SkipFlags(ps.get<IndexVector>("SkipFlags")),
  m_FixedCurvThresh(ps.get<double>("FixedCurvThresh")) {
  const string myname = "AdcCodeMitigator::ctor: ";
  // Display configuration.
  if ( m_LogLevel ) {
    cout << myname << "Configuration: " << endl;
    cout << myname << "          LogLevel: " << m_LogLevel << endl;
    cout << myname << "          FixFlags: [";
    bool first = true;
    for ( Index iflg : m_FixFlags ) {
      if ( ! first ) cout << ", ";
      cout << iflg;
    }
    cout << "]" << endl;
    cout << myname << "  InterpolateFlags: [";
    first = true;
    for ( Index iflg : m_InterpolateFlags ) {
      if ( ! first ) cout << ", ";
      else first = false;
      cout << iflg;
    }
    cout << "]" << endl;
    cout << myname << "         SkipFlags: [";
    first = true;
    for ( Index iflg : m_SkipFlags ) {
      if ( ! first ) cout << ", ";
      else first = false;
      cout << iflg;
    }
    cout << "]" << endl;
    cout << myname << "   FixedCurvThresh: " << m_FixedCurvThresh << endl;
  }
  for ( Index iflg : m_FixFlags ) m_fixSet.insert(iflg);
  for ( Index iflg : m_InterpolateFlags ) {
    m_interpolateSet.insert(iflg);
    m_skipSet.insert(iflg);
  }
  for ( Index iflg : m_SkipFlags ) m_skipSet.insert(iflg);
}

//**********************************************************************

DataMap AdcCodeMitigator::update(AdcChannelData& acd) const {
  const string myname = "AdcCodeMitigator::view: ";
  DataMap ret;
  Index nsam = acd.samples.size();
  if ( acd.flags.size() != acd.samples.size() ) {
    cout << myname << "ERROR: Flag and sample sizes disagree: "
         << acd.flags.size() << " != " << nsam << endl;
    return ret.setStatus(1);
  }
  Index mitCount = 0;
  // Fixed mitigation.
  for ( Index isam=0; isam<nsam; ++isam ) {
    AdcFlag iflg = acd.flags[isam];
    if ( m_fixSet.find(iflg) == m_fixSet.end() ) continue;
    ++mitCount;
    acd.samples[isam] = 0.0;
    acd.flags[isam] = AdcSetFixed;
  }
  // Interpolation mitigation.
  //   Ordering of samples: isamLo2 isamLo1 isam isamHi1 isamHi2
  Index isamLo2 = 0;  // Second sample used for low-side interpolation.
  Index isamLo1 = 0;  // Firstample used for low-side interpolation.
  bool haveLo1 = false;
  bool haveLo2 = false;
  Index isamHi1 = 0;  // Sample used for high-side interpolation.
  Index isamHi2 = 0;  // 2nd Sample used for high-side interpolation.
  for ( Index isam=0; isam<nsam; ++isam ) {
    AdcFlag iflg = acd.flags[isam];
    // Record if this sample can be used for low-side interpolation.
    if ( m_skipSet.find(iflg) == m_skipSet.end() ) {
      isamLo2 = isamLo1;
      isamLo1 = isam;
      haveLo1 = true;
      haveLo2 = isamLo1 != isamLo2;
      continue;
    }
    // Skip sample if it is not to be interpolated.
    if ( m_interpolateSet.find(iflg) == m_interpolateSet.end() ) continue;
    // Check if we have a valid high-side sample.
    bool haveHi1 = isamHi1 > isam && isamHi1 < nsam;
    if ( ! haveHi1 && isamHi1 < nsam ) {
      isamHi1 = isam;
      while ( !haveHi1 && ++isamHi1 < nsam ) haveHi1 = m_skipSet.find(acd.flags[isamHi1]) == m_skipSet.end();
    }
    bool haveHi2 = false;
    if ( haveHi1 ) {
      haveHi2 = isamHi2 > isamHi1 && isamHi2 < nsam;
      if ( ! haveHi2 ) {
        isamHi2 = isamHi1;
        while ( !haveHi2 && ++isamHi2 < nsam ) haveHi2 = m_skipSet.find(acd.flags[isamHi2]) == m_skipSet.end();
        if ( ! haveHi2 ) isamHi2 = isamHi1;
      }
    }
    if ( m_LogLevel >= 5 ) {
      cout << myname << "Samples: " << isamLo2 << " " << isamLo1 << " | " << isam
           << " | " << isamHi1 << " " << isamHi2 << endl;
    }
    double y = 0.0;
    //  If threshold is set and samples are available, consider doing interpolation with
    // a varying slope but constant curvature.
    bool done = false;
    if ( m_FixedCurvThresh > 0.0 && haveLo1 && haveHi1 && ( haveLo2 || haveHi2 ) ) {
      double ylo1 = acd.samples[isamLo1];
      double ylo2 = acd.samples[isamLo2];
      double yhi1 = acd.samples[isamHi1];
      double yhi2 = acd.samples[isamHi2];
      //int nBigJump = 0;
      //if ( fabs(ylo2 - ylo1) > m_FixedCurvThresh ) ++nBigJump;
      //if ( fabs(yhi2 - yhi1) > m_FixedCurvThresh ) ++nBigJump;
      //if ( fabs(yhi1 - ylo1) > m_FixedCurvThresh ) ++nBigJump;
      //bool doFixedCurvature = nBigJump > 1;
      double adiflo = fabs(ylo2 - ylo1);
      double adifhi = fabs(yhi2 - yhi1);
      bool doFixedCurvature = adiflo > m_FixedCurvThresh ||
                              adifhi > m_FixedCurvThresh;
      if ( m_LogLevel >= 5 ) {
        cout << myname << "Lo, hi jumps for sample " << isam << ": " << adiflo << ", " << adifhi << endl;
      }
      if ( doFixedCurvature ) {
        // j = isam - isamLo1
        // k = isam - isamHi1
        // y = a*(i-isamLo1) +  b*(i-isamHi1) + c*(i-isamLo1)*(i-isamHi1)
        //   = a*j + b*k + c*j*k
        // There are three parameters: a, b, c, and four measurements: ylo1, yhi1, ylo2, yhi2
        // Require the  curve pass through the interpolation endpoints ylo1 and yhi1:
        //   a = yhi1/(isamHi1 - isamLo1)
        //   b = ylo1/(isamLo1 - isamHi1)
        // Fix the remaining param c so that
        //   wlo*[y(isamlo2) - ylo2] + whi*[y(isamhi2) - yhi2] = 0
        // E.g. if w1 == w2, the curve is equally close to ylo2 and yhi2
        int jlo2 = int(isamLo2) - int(isamLo1);
        int jhi2 = int(isamHi2) - int(isamLo1);
        int klo2 = int(isamLo2) - int(isamHi1);
        int khi2 = int(isamHi2) - int(isamHi1);
        double a = yhi1/(isamHi1 - isamLo1);
        double b = -ylo1/(isamHi1 - isamLo1);
        //int iden = jlo2*klo2 + jhi2*khi2;
        //double wlo = 1.0;
        //double whi = 1.0;
        double w0 = acd.pedestalRms;
        if ( w0 < 1.0 ) w0 = 1.0;
        double wlo = adifhi + w0;
        double whi = adiflo + w0;
        double den = wlo*jlo2*klo2 + whi*jhi2*khi2;
        if ( den != 0.0 ) {
          //double num = ylo2 + yhi2 - a*(jlo2+jhi2) - b*(klo2+khi2);
          //double c = num/double(iden);
          double num = wlo*ylo2 + whi*yhi2 - a*(wlo*jlo2+whi*jhi2) - b*(wlo*klo2+whi*khi2);
          double c = num/den;
          if ( m_LogLevel >= 4 ) {
            cout << myname << "Quadratic interpolation for sample " << isam << endl;
            cout << myname << "  isamLo2: " << isamLo2 << endl;
            cout << myname << "  isamLo1: " << isamLo1 << endl;
            cout << myname << "  isamHi1: " << isamHi1 << endl;
            cout << myname << "  isamHi2: " << isamHi2 << endl;
            cout << myname << "  ylo2: " << ylo2 << endl;
            cout << myname << "  ylo1: " << ylo1 << endl;
            cout << myname << "  yhi1: " << yhi1 << endl;
            cout << myname << "  yhi2: " << yhi2 << endl;
            cout << myname << "  a: " << a << endl;
            cout << myname << "  b: " << b << endl;
            cout << myname << "  c: " << c << endl;
          }
          int jsam = isam - isamLo1;
          int ksam = isam - isamHi1;
          y = a*jsam + b*ksam + c*jsam*ksam;
          done = true;
        } else {
          if ( m_LogLevel >= 4 ) {
            cout << myname << "Skipping quadratic interpolation for den=0 for sample " << isam << endl;
          }
        }
      }
    }
    // Two point linear interpolation.
    if ( done ) {
      if ( m_LogLevel >= 3 ) cout << myname << "Quadratic interpolation: isam=" << isam
                                  << ": [" << isamLo1 << "], [" << isamHi1 << "] ==> " << y << endl;
    } else if ( haveLo1 && haveHi1 ) {
      double ylo = acd.samples[isamLo1];
      double yhi = acd.samples[isamHi1];
      double x01 = isam - isamLo1;
      double x12 = isamHi1 - isamLo1;
      y = ylo + (yhi - ylo)*x01/x12;
      iflg = AdcInterpolated;
      if ( m_LogLevel >= 3 ) cout << myname << "Linear interpolating isam=" << isam
                                  << ": [" << isamLo1 << "], [" << isamHi1 << "] ==> " << y << endl;
    // One point lo extrapolation.
    } else if ( haveLo1 ) {
      y = acd.samples[isamLo1];
      iflg = AdcExtrapolated;
      if ( m_LogLevel >= 3 ) cout << myname << "Setting low extrapolation for isam=" << isam
                                  << ": [" << isamLo1 << "] = " << y << endl;
    // One point lo extrapolation.
    } else if ( haveHi1 ) {
      y = acd.samples[isamHi1];
      iflg = AdcExtrapolated;
      if ( m_LogLevel >= 3 ) cout << myname << "Setting high extrapolation for isam=" << isam
                                  << ": [" << isamHi1 << "] = " << y << endl;
    // Use fixed value.
    } else {
      iflg = AdcSetFixed;
    }
    acd.samples[isam] = y;
    acd.flags[isam] = iflg;
    ++mitCount;
  }
  ret.setInt("mitCount", mitCount);
  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcCodeMitigator)
