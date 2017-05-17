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
  BasicTrigger::BasicTrigger(bool trigdecision,
			     unsigned int trigtype,
			     unsigned int trigsubtype,
                             Hardware::HardwareID hardwareid):
    fTrigDecision(trigdecision),
    fTrigType(trigtype),
    fTrigSubType(trigsubtype),
    fTrigHardwareID(hardwareid)
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
  void BasicTrigger::setTrigSubType(unsigned int trigsubtype)
  {
    fTrigSubType = trigsubtype;
  }



  //----------------------------------------------------------------------
  unsigned int BasicTrigger::TrigSubType() const
  {
    return fTrigSubType;
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
  void BasicTrigger::setMetrics(std::vector<double> metrics)
  {
    fMetrics = metrics;
  }



  //----------------------------------------------------------------------
  void BasicTrigger::setMetric(double metric)
  {
    fMetrics.push_back(metric);
  }



  //----------------------------------------------------------------------
  std::vector<double> BasicTrigger::Metrics() const
  {
    return fMetrics;
  }



  //----------------------------------------------------------------------
  double BasicTrigger::Metric(unsigned int i) const
  {
    if(fMetrics.size() > i) return fMetrics[i];
    else {
      std::cout << "\n\n\nWARNING!!! Requested metric index is larger than the size of fMetrics!!!\n\n\n";
      return 0;
    }
  }



  //----------------------------------------------------------------------
  // ostream operator.
  //
#ifndef __GCCXML__
  std::ostream& operator << (std::ostream& o, BasicTrigger const& bt) {

    o << "Trigger Type       = " << bt.TrigType()          << std::endl;
    o << "Trigger Sub-Type   = " << bt.TrigSubType()       << std::endl;
    o << "Trigger Decision   = " << bt.TrigDecision()      << std::endl;
    o << "Trigger HardwareID = " << bt.fTrigHardwareID     << std::endl;
    return o;

  }



  //----------------------------------------------------------------------
  // < operator.
  //
  bool operator < (BasicTrigger const& a, BasicTrigger const& b) {

    return a.TrigType() < b.TrigType();

  }
#endif //__GCCXML__  
} // end namespace triggersim
