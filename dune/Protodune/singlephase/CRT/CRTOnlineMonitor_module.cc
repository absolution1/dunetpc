////////////////////////////////////////////////////////////////////////
// Class:       CRTOnlineMonitor
// Plugin Type: analyzer (art v2_11_03)
// File:        CRTOnlineMonitor_module.cc
// Brief:       Runs the OnlinePlotter algorithm to make online monitoring 
//              plots for the CRT.  Provides a wrapper to use OnlinePlotter 
//              within ART.   
//
// Generated at Tue Aug 21 03:22:48 2018 by Andrew Olivier using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

//CRT includes
#include "dunetpc/dune/Protodune/singlephase/CRT/alg/monitor/OnlinePlotter.cpp"

//c++ includes
#include <memory> //For std::unique_ptr

class CRTOnlineMonitor;

class CRTOnlineMonitor : public art::EDAnalyzer {
public:
  explicit CRTOnlineMonitor(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTOnlineMonitor(CRTOnlineMonitor const &) = delete;
  CRTOnlineMonitor(CRTOnlineMonitor &&) = delete;
  CRTOnlineMonitor & operator = (CRTOnlineMonitor const &) = delete;
  CRTOnlineMonitor & operator = (CRTOnlineMonitor &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run const & r) override;
  void endRun(art::Run const & r) override;

private:

  // Use the PIMPL idiom to weaken coupling between this module and the 
  // OnlinePlotter source code.  
  std::unique_ptr<CRT::OnlinePlotter<art::ServiceHandle<art::TFileService>>> fPlotter; //Algorithm that makes online monitoring plots 
                                                                                       //from CRT::Triggers

  art::InputTag fCRTLabel; //The name of the module that created the CRT::Triggers this module will read
};


CRTOnlineMonitor::CRTOnlineMonitor(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p), fPlotter(nullptr), fCRTLabel(p.get<art::InputTag>("CRTLabel"))
 // More initializers here.
{
  consumes<std::vector<CRT::Trigger>>(fCRTLabel);
}

void CRTOnlineMonitor::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  const auto& triggers = e.getValidHandle<std::vector<CRT::Trigger>>(fCRTLabel);
  fPlotter->AnalyzeEvent(*triggers);
}

void CRTOnlineMonitor::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fPlotter = std::make_unique<CRT::OnlinePlotter<art::ServiceHandle<art::TFileService>>>(tfs);
}

void CRTOnlineMonitor::beginRun(art::Run const & r)
{
  // Implementation of optional member function here.
  fPlotter->ReactBeginRun(""); //TODO: Remove std::string from interface for OnlinePlotter and friends
}

void CRTOnlineMonitor::endRun(art::Run const & r)
{
  // Implementation of optional member function here.
  fPlotter->ReactEndRun(""); //TODO: Remove std::string from interface for OnlinePlotter and friends
}

DEFINE_ART_MODULE(CRTOnlineMonitor)
