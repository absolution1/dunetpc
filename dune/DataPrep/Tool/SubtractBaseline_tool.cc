// SubtractBaseline_tool.cc

#include "SubtractBaseline.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/LArFFT.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"

using std::string;
using std::cout;
using std::endl;

//**********************************************************************
// Class methods.
//**********************************************************************

SubtractBaseline::SubtractBaseline(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel"))
, m_BaseSampleBins(ps.get<int>("BaseSampleBins"))
, m_BaseVarCut(ps.get<float>("BaseVarCut")) {
  const string myname = "SubtractBaseline::ctor: ";
  if ( m_LogLevel >= 1 ) {
    cout << myname << "Parameters:" << endl;
    cout << myname << "  LogLevel: " << m_LogLevel << endl;
    cout << myname << "  BaseSampleBins: "<< m_BaseSampleBins << endl;
    cout << myname << "  BaseVarCut: "<< m_BaseVarCut << endl;
  }
}

//**********************************************************************

DataMap SubtractBaseline::update(AdcChannelData& acd) const {
  const string myname = "SubtractBaseline::view: ";
  DataMap ret;
  AdcSignalVector& samples = acd.samples;

  // number of points to characterize the baseline
  unsigned short nBasePts = 1 + samples.size() / m_BaseSampleBins;
    // the baseline offset vector
    std::vector<float> base(nBasePts, 0.);
    // find the average value in each region, using values that are
    // similar
    float fbins = m_BaseSampleBins;
    unsigned short nfilld = 0;
    for(unsigned short ii = 0; ii < nBasePts; ++ii) {
      unsigned short loBin = ii * m_BaseSampleBins;
      unsigned short hiBin = loBin + m_BaseSampleBins;
      float ave = 0.;
      float sum = 0.;
      for(unsigned short bin = loBin; bin < hiBin; ++bin) {
        ave += samples[bin];
        sum += samples[bin] * samples[bin];
      } // bin
      ave = ave / fbins;
      float var = (sum - fbins * ave * ave) / (fbins - 1.);
      // Set the baseline for this region if the variance is small
      if(var < m_BaseVarCut) {
        base[ii] = ave;
        ++nfilld;
      }
    } // ii
    
    // fill in any missing points if there aren't too many missing
    if(nfilld < nBasePts && nfilld > nBasePts / 2) {
      bool baseOK = true;
      // check the first region
      if(base[0] == 0) {
        unsigned short ii1 = 0;
        for(unsigned short ii = 1; ii < nBasePts; ++ii) {
          if(base[ii] != 0) {
            ii1 = ii;
            break;
          } // base[ii] != 0
        } // ii
        unsigned short ii2 = 0;
        for(unsigned short ii = ii1 + 1; ii < nBasePts; ++ii) {
          if(base[ii] != 0) {
            ii2 = ii;
            break;
          } // base[ii] != 0
        } // ii
        // failure
        if(ii2 > 0) {
          float slp = (base[ii2] - base[ii1]) / (float)(ii2 - ii1);
          base[0] = base[ii1] - slp * ii1;
        } else {
          baseOK = false;
        }
      } // base[0] == 0
      // check the last region
      if(baseOK && base[nBasePts] == 0) {
        unsigned short ii2 = 0;
        for(unsigned short ii = nBasePts - 1; ii > 0; --ii) {
          if(base[ii] != 0) {
            ii2 = ii;
            break;            
          }
        } // ii
        baseOK = false; // assume the worst, hope for better
        unsigned short ii1 = 0;
        if (ii2 >= 1) {
          for(unsigned short ii = ii2 - 1; ii > 0; --ii) {
            if(base[ii] != 0) {
              ii1 = ii;
              baseOK = true;
              break;
            } // if base[ii]
          } // ii
        } // if ii2
        if (baseOK) {
          float slp = (base[ii2] - base[ii1]) / (float)(ii2 - ii1);
          base[nBasePts] = base[ii2] + slp * (nBasePts - ii2);
        }
      } // baseOK && base[nBasePts] == 0
      // now fill in any intermediate points
      for(unsigned short ii = 1; ii < nBasePts - 1; ++ii) {
        if(base[ii] == 0) {
          // find the next non-zero region
          for(unsigned short jj = ii + 1; jj < nBasePts; ++jj) {
            if(base[jj] != 0) {
              float slp = (base[jj] - base[ii - 1]) / (jj - ii + 1);
              base[ii] = base[ii - 1] + slp;
              break;
            }
          } // jj
        } // base[ii] == 0
      } // ii
    } // nfilld < nBasePts
    // interpolate and subtract
    float slp = (base[1] - base[0]) / (float)m_BaseSampleBins;
    // bin offset to the origin (the center of the region)
    short bof = m_BaseSampleBins / 2;
    short lastRegion = 0;
    for(unsigned short bin = 0; bin < samples.size(); ++bin) {
      // in a new region?
      short region = bin / m_BaseSampleBins;
      if(region > lastRegion) {
        // update the slope and offset
        slp = (base[region] - base[lastRegion]) / (float)m_BaseSampleBins;
        bof += m_BaseSampleBins;
        lastRegion = region;
      }
      samples[bin] -= base[region] + (bin - bof) * slp;
    } // bin

  return ret;
}

//**********************************************************************
