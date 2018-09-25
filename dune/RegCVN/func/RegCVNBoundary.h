////////////////////////////////////////////////////////////////////////
/// \file    RegCVNBoundary.h
/// \brief   RegCVNBoundary for RegCVN PixelMap modified from CVNBoundary.h
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef REGCVN_BOUNDARY_H
#define REGCVN_BOUNDARY_H

#include <ostream>
#include <vector>


namespace cvn
{


  /// RegCVNBoundary object intended for use with cvn::RegPixelMap.  
  /// Stores mean of wire and TPCs
  class RegCVNBoundary
  {

  public:
    /// Create new RegCVNBoundary object based on number of wires, number of tdcs,
    /// minumum wire and mean tdc in odd and even view.
    RegCVNBoundary(const int& nWire, const int& nTDC, const int& tRes,
             const int& WireMeanX,
             const int& WireMeanY,
             const int& WireMeanZ,
             const int& TDCMeanX,
             const int& TDCMeanY,
             const int& TDCMeanZ);

    RegCVNBoundary(){};

    bool IsWithin(const int& wire, const int& tdc, const unsigned int& view);

    int FirstWire(const unsigned int& view) const {return fFirstWire[view];};
    int LastWire(const unsigned int& view) const {return fLastWire[view];};
    int FirstTDC(const unsigned int& view) const {return fFirstTDC[view];};
    int LastTDC (const unsigned int& view) const {return fLastTDC[view];};



  private:
    int fFirstWire[3];  ///< Minimum wire, inclusive
    int fLastWire[3];   ///< Maximum wire, inclusive
    int fFirstTDC[3]; ///< Minimum TDC in each view, inclusive
    int fLastTDC[3];  ///< Maximum TDC in each view, inclusive


  };

  std::ostream& operator<<(std::ostream& os, const RegCVNBoundary& b);
}

#endif  // REGCVN_BOUNDARY_H
