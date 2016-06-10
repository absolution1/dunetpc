// ActivityTriggerService_service.cc
//
// Michael Baird
// May 2016
//
// Basic trigger algorithm that looks for some minimal activity within
// a single APA.
//

// ART and larsoft includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "larcore/SimpleTypesAndConstants/RawTypes.h"

// DUNETPC specific includes
#include "dune/DAQTriggerSim/Service/ActivityTriggerService.h"

// Framework includes
#include "fhiclcpp/ParameterSet.h"

// C++ includes
#include <iostream>



//.....................................................................
ActivityTriggerService::ActivityTriggerService() :
  fMinHits(0),
  fThreshold(0),
  fMinTicks(0),
  fNHits(0)
{ }



//.....................................................................
ActivityTriggerService::ActivityTriggerService(const fhicl::ParameterSet& pset, art::ActivityRegistry&) :
  fMinHits(0),
  fThreshold(0),
  fMinTicks(0),
  fNHits(0)
{
  fMinHits   = pset.get<unsigned int>("minHits");
  fThreshold = pset.get<short>       ("threshold");
  fMinTicks  = pset.get<unsigned int>("minTicks");
}



//.....................................................................
bool ActivityTriggerService::TPCTrigger(std::vector<raw::RawDigit> const & raw)
{

  this->Reset();
  bool decision = false;
  
  // make the geometry service handle
  art::ServiceHandle<geo::Geometry> geom;
  
  // loop over the channels:
  std::cout << "\n\n\n";
  for(unsigned int i = 0; i < raw.size(); ++i) {

    // Only use the collection wires for hit finding. Skip all other wire types.
    if(geom->SignalType(raw[i].Channel()) != geo::kCollection) continue;

    // loop over the vector of ADC values for this wire and compute
    // the average and sd ADC values
    std::vector<short> ADCvec = raw[i].ADCs();
    double ADCave = 0.0;
    double ADCsd  = 0.0;
    short  ADCmax = 0;
    for(unsigned int j = 0; j < ADCvec.size(); ++j) {
      ADCave += (double)ADCvec[j];
      ADCsd  += ((double)ADCvec[j])*((double)ADCvec[j]);
      if(ADCvec[j] > ADCmax) ADCmax = ADCvec[j];
    }
    if(ADCvec.size() > 0) {
      ADCave = ADCave/(double)ADCvec.size();
      ADCsd  = sqrt(ADCsd/(double)ADCvec.size() - ADCave*ADCave);
    }

    /*
    std::cout << "\nADCave     = " << ADCave
	      << "\nADCsd      = " << ADCsd
	      << "\nADCmax     = " << ADCmax
	      << "\nADC.size() = " << ADCvec.size()
	      << "\n";
    */
    
  } // end loop over channels (i)
  std::cout << "\n\n\n";



  return decision;
}



//.....................................................................
void ActivityTriggerService::Reset()
{
  fNHits = 0;
}


DEFINE_ART_SERVICE_INTERFACE_IMPL(ActivityTriggerService, BaseTriggerService)
