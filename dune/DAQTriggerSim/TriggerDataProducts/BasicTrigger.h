// BasicTrigger.h
//
// Michael Baird
// March 2016
//
// Basic trigger data product class for the DAQ trigger service framework.
//

#ifndef BASICTRIGGER_H
#define BASICTRIGGER_H

// DUNETPC specific includes
#include "dune/DAQTriggerSim/TriggerDataProducts/TriggerTypes.h"

// Framework includes

// C++ includes

namespace triggersim {

  class BasicTrigger {
    
  public:
    
    // Constructor...
    BasicTrigger(unsigned int trigtype = kNullTrigger,
		 bool trigissued = false);
    
    // Destructor...
    ~BasicTrigger();

    // Set and Get functions for the trigger type
    void setTrigType(unsigned int trigtype);
    unsigned int TrigType();

    // Set and Get functions for the trigger issued variable
    void setTrigIssued(bool trigissued);
    bool TrigIssued();

  private:
    
    // Parameters.
    unsigned int fTrigType;
    bool         fTrigIssued;
    
  };
} // end namespace triggersim

#endif
