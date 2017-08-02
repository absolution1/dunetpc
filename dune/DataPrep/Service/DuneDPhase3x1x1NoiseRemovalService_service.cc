// DuneDPhase3x1x1NoiseRemovalService_service.cc

#include "DuneDPhase3x1x1NoiseRemovalService.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"

//**********************************************************************

DuneDPhase3x1x1NoiseRemovalService::
DuneDPhase3x1x1NoiseRemovalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&) :
    fGeometry( &*art::ServiceHandle<geo::Geometry>() )
{
}
//**********************************************************************

int DuneDPhase3x1x1NoiseRemovalService::update(AdcChannelDataMap& datamap) const {
  const std::string myname = "DuneDPhase3x1x1NoiseRemovalService:update: ";
  if ( datamap.size() == 0 ) {
    std::cout << myname << "WARNING: No channels found." << std::endl;
    return 1;
  }

  unsigned int nsig = 0;
  bool first = true;
  for (const AdcChannelDataMap::value_type& ent : datamap)
  {
    const AdcChannelData& data = ent.second;
    if (first) { nsig = data.samples.size(); first = false; }
    else if (data.samples.size() != nsig)
    {
      std::cout << myname << "WARNING: Channels have inconsistent sample counts." << std::endl;
      return 2;
    }
  }
  
  if (nsig == 0)
  {
    std::cout << myname << "WARNING: No ADC samples found." << std::endl;
    return 3;
  }

  std::cout << myname << "Processing..." << std::endl;

  return 0;
}
//**********************************************************************

std::ostream& DuneDPhase3x1x1NoiseRemovalService::print(std::ostream& out, std::string prefix) const
{
  out << prefix << "DuneDPhase3x1x1NoiseRemovalService:  ...info" << std::endl;
  return out;
}
//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DuneDPhase3x1x1NoiseRemovalService, AdcNoiseRemovalService)

//**********************************************************************
