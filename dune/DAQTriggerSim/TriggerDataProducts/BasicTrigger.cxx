// BasicTrigger.cxx
//
// Michael Baird
// March 2016
//
// Basic trigger data product class for the DAQ trigger service framework.
//

// DUNETPC specific includes
#include "dune/DAQTriggerSim/TriggerDataProducts/BasicTrigger.h"

// Framework includes

// C++ includes
#include <iomanip>

namespace triggersim {

  //----------------------------------------------------------------------
  BasicTrigger::BasicTrigger(unsigned int trigtype, bool trigdecision):
    fTrigType(trigtype),
    fTrigDecision(trigdecision)
  {}


  
  //----------------------------------------------------------------------
  BasicTrigger::~BasicTrigger()
  {}


  
  //----------------------------------------------------------------------
  void BasicTrigger::setTrigType(unsigned int trigtype)
  {
    fTrigType = trigtype;
  }



  //----------------------------------------------------------------------
  unsigned int BasicTrigger::TrigType() const
  {
    return fTrigType;
  }



  //----------------------------------------------------------------------
  void BasicTrigger::setTrigDecision(bool trigdecision)
  {
    fTrigDecision = trigdecision;
  }



  //----------------------------------------------------------------------
  bool BasicTrigger::TrigDecision() const
  {
    return fTrigDecision;
  }


  //----------------------------------------------------------------------
  // ostream operator.
  //
  std::ostream& operator << (std::ostream& o, BasicTrigger const& bt) {

    o << "Trigger Type     = " << bt.TrigType()     << std::endl;
    o << "Trigger Decision = " << bt.TrigDecision() << std::endl;
    return o;

  }



  //----------------------------------------------------------------------
  // < operator.
  //
  bool operator < (BasicTrigger const& a, BasicTrigger const& b) {

    return a.TrigType() < b.TrigType();

  }
  
} // end namespace triggersim
