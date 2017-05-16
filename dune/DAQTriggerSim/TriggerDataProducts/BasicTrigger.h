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
#include "dune/DAQTriggerSim/TriggerDataProducts/HardwareElements.h"

// Framework includes

// C++ includes
#include <vector>
#include <iostream>

#ifndef __GCCXML__
#include <iosfwd> // std::ostream
#endif


namespace triggersim {

  class BasicTrigger {
    
  public:
    
    // Constructor...
    BasicTrigger(bool trigdecision = false,
		 unsigned int trigtype = kNullTrigger,
		 unsigned int trigsubtype = kNullTrigger,
                 Hardware::HardwareID hardwareid = Hardware::UnknownID);
    
    // Destructor...
    ~BasicTrigger();

    // Set and Get functions for the trigger type
    void setTrigType(unsigned int trigtype);
    unsigned int TrigType() const;

    // Set and Get functions for the trigger sub-type
    void setTrigSubType(unsigned int trigsubtype);
    unsigned int TrigSubType() const;

    // Set and Get functions for the trigger decision variable
    void setTrigDecision(bool trigdecision);
    bool TrigDecision() const;

    // Set and Get functions for the trigger algorithm performance metrics
    void setMetrics(std::vector<double> metrics);
    void setMetric(double metric);
    std::vector<double> Metrics() const;
    double Metric(unsigned int i) const;

    void setTrigHardwareID(Hardware::HardwareID hardwareid){ fTrigHardwareID = hardwareid; }
    Hardware::HardwareID const& TrigHardwareID() const{ return fTrigHardwareID; }


#ifndef __GCCXML__
    friend std::ostream& operator << (std::ostream& o, BasicTrigger const& bt);
    friend bool          operator <  (BasicTrigger const& a, BasicTrigger const& b);    
#endif //__GCCXML__

  private:
    

    // Parameters:
    bool         fTrigDecision;    ///< did the trigger algorithm decide to keep this event?
    unsigned int fTrigType;        ///< trigger type (types and subtypes are defined in TriggerTypes.h)
    unsigned int fTrigSubType;     ///< trigger subtype (types and subtypes are defined in TriggerTypes.h)
    Hardware::HardwareID fTrigHardwareID;

    std::vector<double> fMetrics;  ///< a vector to store trigger algorithm performance metrics
  };
} // end namespace triggersim

#endif
