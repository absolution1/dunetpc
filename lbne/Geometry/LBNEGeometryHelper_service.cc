////////////////////////////////////////////////////////////////////////////////
/// \file LBNEGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

// Migration note:
// Geometry --> lbne/Geometry
#include "lbne/Geometry/LBNEGeometryHelper.h"

#include "Geometry/ChannelMapAlg.h"

// Migration note:
// Geometry --> lbne/Geometry for the two below
#include "lbne/Geometry/ChannelMap35Alg.h"
#include "lbne/Geometry/ChannelMap35OptAlg.h"
#include "lbne/Geometry/ChannelMapAPAAlg.h"

#include "TString.h"


namespace lbne
{

  LBNEGeometryHelper::LBNEGeometryHelper( fhicl::ParameterSet const & pset, art::ActivityRegistry & reg )
  :  fPset( pset ),
     fReg( reg ),
     fChannelMap()
  {}

  LBNEGeometryHelper::~LBNEGeometryHelper() throw()
  {}  
  
  void LBNEGeometryHelper::doConfigureChannelMapAlg( const TString & detectorName,
                                                     fhicl::ParameterSet const & sortingParam,
                                                     std::vector<geo::CryostatGeo*> & c )
  {
    fChannelMap = nullptr;
    
    if ( detectorName.Contains("lbne35t") ) 
    {
// Migration note:
// Change all geo::ChannelMap... to lbne::ChannelMap... throughout file

      TString detectorVersion = sortingParam.get< std::string >("DetectorVersion");

      if ( detectorVersion.Contains("v3") )
	fChannelMap = std::shared_ptr<geo::ChannelMapAlg>( new geo::ChannelMap35OptAlg( sortingParam ) );
      else
	fChannelMap = std::shared_ptr<geo::ChannelMapAlg>( new geo::ChannelMap35Alg( sortingParam ) );
    }
    else if ( detectorName.Contains("lbne10kt") ) 
    {
      fChannelMap = std::shared_ptr<geo::ChannelMapAlg>( new geo::ChannelMapAPAAlg( sortingParam ) );
    }
    else if ( detectorName.Contains("lbne34kt") )
    {
      fChannelMap = std::shared_ptr<geo::ChannelMapAlg>( new geo::ChannelMapAPAAlg( sortingParam ) );
    }
    else
    {
      fChannelMap = nullptr;
    }
    if ( fChannelMap )
    {
      fChannelMap->Initialize( c );
    }
  }
  
  std::shared_ptr<const geo::ChannelMapAlg> LBNEGeometryHelper::doGetChannelMapAlg() const
  {
    return fChannelMap;
  }

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(lbne::LBNEGeometryHelper, geo::ExptGeoHelperInterface)
