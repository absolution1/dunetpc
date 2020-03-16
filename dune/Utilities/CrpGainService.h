////////////////////////////////////////////////////////////////////////
// Class:       CrpGainService
// Plugin Type: service (art v3_02_06)
// File:        CrpGainService.h
//
// Generated at Thu Nov  7 16:24:55 2019 by Vyacheslav Galymov using cetskelgen
// from cetlib version v3_07_02.
//
// This services calculates charge amplification factors for charge readout 
// planes (CRPs) of dual-phase LArTPC
//  - LEM gain is calculated based on 2D coordinate of the 
//    projected (drifted) charge.
//  - LEM borders and dead areas can be taken into account
// 
//
////////////////////////////////////////////////////////////////////////

#ifndef CrpGainService_h
#define CrpGainService_h

#include <vector>
#include <string>

#include "larcore/Geometry/Geometry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace fhicl {
  class ParameterSet;
}

namespace art {
  class ActivityRegistry;
}

namespace geo {
  class Geometry;
}

namespace sim {
  class SimChannel;
}

namespace util {
  class CrpGainService;
}

class FloatArrayTool;


namespace util
{

//
class CrpGainService {
public:
  using FloatArrayPtr = const FloatArrayTool*;
  
  explicit CrpGainService(fhicl::ParameterSet const& ps, art::ActivityRegistry& areg);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  
  // get charge collected on a view after amplification in CRP
  double viewCharge( const sim::SimChannel* psc, unsigned itck ) const;
  
  // calculate the gain based on position information
  double crpGain( geo::Point_t const &pos ) const;

  // default value of the effective gain
  double crpDefaultGain() const { return m_CrpDefGain; }
  
private:
  // methods
  bool   checkGeoConfig() const;
  int    getLemId( unsigned crp, int chx, int chy ) const;
  double getCrpGain( unsigned crp, int chx, int chy ) const;
  double getLemTransparency( int chx, int chy ) const;
  double getLemGain( int lemid ) const;

  void   dumpLemEffMap(int nlems) const;

  //
  int    m_LogLevel;
  
  // default value of CRP gain effective gain
  double m_CrpDefGain;
  bool   m_UseDefGain;

  //
  unsigned int m_CrpNLem;
  unsigned int m_CrpNLemPerSide;
  
  //
  unsigned int m_LemViewChans;
  unsigned int m_LemTotChans;

  // LEM transparency tool
  std::string m_LemEffTool;
  FloatArrayPtr m_plemeff;

  // dummy for now ...
  std::vector<float> m_lemgainmap;

  // detector geometry
  const geo::Geometry* m_geo;
};

}
DECLARE_ART_SERVICE(util::CrpGainService, LEGACY)
#endif /* CrpGainService_h */
