////////////////////////////////////////////////////////////////////////
// Class:       TimeOffset
// Plugin Type: analyzer (art v2_11_03)
// File:        TimeOffset_module.cc
// Brief:       A module to quantify time offsets between ProtoDUNE-SP 
//              CRT raw data and TPC raw data.  This module shall produce 
//              a plot of time difference between CRT timestamps and 
//              RDTimeStamps in an event.  I expect lots of "false matches", 
//              but a peak in this plot indicates a relative time offset 
//              between the CRT and the TPC raw data.  This module will 
//              also produce diagnostics for CRT timing.  
//
// Generated at Tue Oct 16 08:06:54 2018 by Andrew Olivier using cetskelgen
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

//lardataobj includes
#include "lardataobj/RawData/RDTimeStamp.h"

namespace CRT {
  class TimeOffset;
}

class CRT::TimeOffset : public art::EDAnalyzer {
public:
  explicit TimeOffset(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TimeOffset(TimeOffset const &) = delete;
  TimeOffset(TimeOffset &&) = delete;
  TimeOffset & operator = (TimeOffset const &) = delete;
  TimeOffset & operator = (TimeOffset &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Create all of the histograms filled by this module.  Callback 
  // for TFileService in case file switching is enabled. 
  void onFileClose(); 

  // Labels to identify instances of CRT raw data and 
  // RDTimeStamps to be compared.  Use lar -c eventdump.fcl 
  // on a prospective input file to find what labels to 
  // use.  An art::InputTag consists of the label of a 
  // module that put() a data product into the Event and 
  // an optional instance name. 
  const art::InputTag fCRTLabel; //Instance of CRT::Triggers to read
  const art::InputTag fTimestampLabel; //Instance of RDTimeStamps to read 

  // Plots that may be produced for each input file depending on 
  // how TFileService is configured.  All pointers are observer 
  // pointers to objects owned by a TDirectory that is managed 
  // by a TFileService.  
  TH1D* fTimestampMinusCRT; // Time difference between each RDTimeStamp 
                            // and CRT::Trigger
};

CRT::TimeOffset::TimeOffset(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fCRTLabel(p.get<art::InputTag>("CRTLabel", "crt")), 
  fTimestampLabel(p.get<art::InputTag>("TimestampLabel"))
{
  // Tell "scheduler" which data products this module needs as input
  consumes<std::vector<CRT::Trigger>>(fCRTLabel);
  consumes<std::vector<raw::RDTimeStamp>>(fTimestampLabel);

  // Set up file change callback
  art::ServiceHandle<art::TFileService> tfs;
  tfs->registerFileSwitchCallback(this, &CRT::TimeOffset::onFileClose);
  onFileClose(); // Setting up the callback above doesn't actually call 
                 // it for the first file.  So, call it now to create 
                 // initial histograms.  
}

void CRT::TimeOffset::onFileClose()
{
  std::stringstream ss;
  ss << fTimestampLabel;

  art::ServiceHandle<art::TFileService> tfs;
  fTimestampMinusCRT = tfs->make<TH1D>("TimestampMinusCRT", (std::string("Time Difference Between CRT and ")+ss.str()+";Pairs;Time [ticks]").c_str(), 
                                       1000, -125000, 125000); //125000 is width of the DAQ readout window 
}

void CRT::TimeOffset::analyze(art::Event const & e)
{
  // Get the CRT::Triggers and RDTimeStamps that we will compare.  
  // If we fail to get either data product, skip this event without 
  // killing the job.  Propagate an error message to the user about 
  // which data product wasn't found first.  
  try
  {
    const auto& crtHandle = e.getValidHandle<std::vector<CRT::Trigger>>(fCRTLabel);
    const auto& timeHandle = e.getValidHandle<std::vector<raw::RDTimeStamp>>(fTimestampLabel);

    for(const auto& trigger: *crtHandle)
    {
      //Fill plots for all combinations of a CRT::Trigger and an RDTimeStamp
      for(const auto& time: *timeHandle)
      {
        fTimestampMinusCRT->Fill(time.GetTimeStamp() - trigger.Timestamp());
      } //For each RDTimeStamp from fTimestampLabel

      //TODO: Fill plots for this CRT::Trigger
    } //For each CRT::Trigger from fCRTLabel
  }
  catch(const cet::exception& e)
  {
    mf::LogWarning("Product not found") << "Failed to find one of the data products that "
                                        << "CRT::TimeOffset needs to make timing plots.  "
                                        << "This module needs CRT::Triggers and raw::RDTimeStamps.\n"
                                        << e.what() << "\n";
  }
}

//TODO: Decide which plots to create and write file change callback

void CRT::TimeOffset::beginJob()
{
  //TODO: Set up file chance callback here?
  //TODO: Call file change callback here to create plots for first time.  
}

DEFINE_ART_MODULE(CRT::TimeOffset)
