////////////////////////////////////////////////////////////////////////
// \file SlowControlsProtoDUNE.h
//
// \brief header of class for accessing slow controls data for ProtoDUNE
//
// \author jpaley@fnal.gov
// 
////////////////////////////////////////////////////////////////////////
#ifndef SLOW_CONTROLS_PROTODUNE_H
#define SLOW_CONTROLS_PROTODUNE_H


// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <vector>
#include <map>

// dunetpc includes
#include "dune/SlowControls/SlowControls.h"

namespace slowctrls {

  class SlowControlsProtoDUNE : public SlowControls {
    
  public:

    SlowControlsProtoDUNE();
    SlowControlsProtoDUNE(fhicl::ParameterSet const& pset);
    SlowControlsProtoDUNE(SlowControlsProtoDUNE const&) = delete;
    virtual ~SlowControlsProtoDUNE() = default;
      
    bool Configure(fhicl::ParameterSet const& pset);
    bool Update(float ts=0);
    
    virtual double GetValue(std::string& chan, float t) override;
    void SetTimeWindow(int dt);
    void SetSlowCtrlFileName(std::string fname) { fSlowCtrlFileName = fname; }
    void SetUseCondb(bool v) { fUseCondb = v; }
    void SetVerbosity(int v) { fVerbosity = v; }

    void AddChannel(std::string ch, int id) { 
      fNameToIdMap[ch] = id; 
      fIdToNameMap[id] = ch;
    }

  protected:
    bool Load();

  protected:
    bool fUseCondb;
    bool fIsLoaded;
    int  fTimeWindow;
    int  fVerbosity;
    float fCurrentTS;
    std::string fSlowCtrlFileName;

    std::map<std::string, int> fNameToIdMap;
    std::map<int, std::string> fIdToNameMap;
    std::map<int, std::vector<double> > fValues;
    std::map<int, std::vector<float> > fTimes;

  }; // class SlowControlsProtoDUNE
} //namespace slowctrls
#endif // SLOW_CONTROLS_PROTODUNE_H
