/////////////////////////////////////////////////////////
//
// DAQTriggerSim/TriggerDataProducts/classes.h
//
// Michael Baird
// March 2016
//
/////////////////////////////////////////////////////////

// Framework includes:
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"

// DUNETPC includes:
#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"



// make dummy variables of the PtrVector types so that we are
// sure to generate the dictionaries for them
namespace {
  art::PtrVector<triggersim::BasicTrigger> basictrig;
}



//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

template class std::vector<    triggersim::BasicTrigger >;
template class art::Ptr<       triggersim::BasicTrigger >;
template class art::PtrVector< triggersim::BasicTrigger >;

// NOTE: May need template classes for the following in the future:
//
//     std::pair< art::Ptr<triggersim::BasicTrigger>, art::Ptr< {something else...} >
//     art::Assns< art::Ptr<triggersim::BasicTrigger>, art::Ptr< {something else...}, void >
//
//     ...and the art::Wrappers for these things too...

template class art::Wrapper< std::vector<triggersim::BasicTrigger>                  >;
template class art::Wrapper< std::vector<art::PtrVector<triggersim::BasicTrigger> > >;
