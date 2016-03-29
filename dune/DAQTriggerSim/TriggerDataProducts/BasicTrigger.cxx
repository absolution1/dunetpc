// BasicTrigger.cxx
//
// Michael Baird
// March 2016
//
// Basic trigger data product class for the DAQ trigger service framework.
//

// DUNETPC specific includes

// Framework includes

// C++ includes

namespace triggersim {

  //----------------------------------------------------------------------
  BasicTrigger::BasicTrigger(unsigned int trigtype, bool trigissued):
    fTrigType(trigtype),
    fTrigIssued(trigissued)
  {}


  
  //----------------------------------------------------------------------
  ~BasicTrigger::BasicTrigger()
  {}


  
  //----------------------------------------------------------------------
  void BasicTrigger::setTrigType(unsigned int trigtype)
  {
    fTrigType = trigtype;
  }



  //----------------------------------------------------------------------
  unsigned int BasicTrigger::TrigType()
  {
    return fTrigType;
  }



  //----------------------------------------------------------------------
  void BasicTrigger::setTrigIssued(bool trigissued)
  {
    fTrigIssued = trigissued;
  }



  //----------------------------------------------------------------------
  bool BasicTrigger::TrigIssued()
  {
    return fTrigIssued;
  }
  
} // end namespace triggersim
