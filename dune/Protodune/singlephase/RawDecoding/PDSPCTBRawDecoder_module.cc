////////////////////////////////////////////////////////////////////////
// Class:       PDSPCTBRawDecoder
// Plugin Type: producer (art v2_11_03)
// File:        PDSPCTBRawDecoder_module.cc
//
// Generated at Mon Aug 13 16:42:23 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

// artdaq and dune-raw-data includes

#include "dune-raw-data/Overlays/CTBFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
#include "dune-raw-data/Overlays/FragmentType.hh"

// dunetpc includes

#include "dune/Protodune/singlephase/CTB/data/pdspctb.h"

class PDSPCTBRawDecoder;


class PDSPCTBRawDecoder : public art::EDProducer {
public:
  explicit PDSPCTBRawDecoder(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPCTBRawDecoder(PDSPCTBRawDecoder const &) = delete;
  PDSPCTBRawDecoder(PDSPCTBRawDecoder &&) = delete;
  PDSPCTBRawDecoder & operator = (PDSPCTBRawDecoder const &) = delete;
  PDSPCTBRawDecoder & operator = (PDSPCTBRawDecoder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  // Declare member data here.

  std::string fInputLabel;
  std::string fInputContainerInstance;
  std::string fInputNonContainerInstance;
  std::string fOutputLabel;

  void _process_CTB_AUX(const artdaq::Fragment& frag);

  std::vector<raw::ctb::Trigger> fTrigs;
  std::vector<raw::ctb::ChStatus> fChStats;
  std::vector<raw::ctb::Feedback> fFeedbacks;
  std::vector<raw::ctb::Misc> fMiscs;
  std::vector<raw::ctb::WordIndex> fWordIndexes;

};


PDSPCTBRawDecoder::PDSPCTBRawDecoder(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  fInputLabel = p.get<std::string>("InputLabel");
  fInputContainerInstance = p.get<std::string>("InputContainerInstance");
  fInputNonContainerInstance = p.get<std::string>("InputNonContainerInstance");
  fOutputLabel = p.get<std::string>("OutputLabel");

  produces<std::vector<raw::ctb::pdspctb> >(fOutputLabel);
}

void PDSPCTBRawDecoder::produce(art::Event & evt)
{

  fTrigs.clear();
  fChStats.clear();
  fFeedbacks.clear();
  fMiscs.clear();
  fWordIndexes.clear();

  // look first for container fragments and then non-container fragments

  std::vector<raw::ctb::pdspctb> pdspctbs;

  art::Handle<artdaq::Fragments> cont_frags;
  evt.getByLabel(fInputLabel, fInputContainerInstance, cont_frags);  

  if(cont_frags.isValid())
    {
      for (auto const& cont : *cont_frags)
	{
	  artdaq::ContainerFragment cont_frag(cont);
	  for (size_t ii = 0; ii < cont_frag.block_count(); ++ii)
	    {
	      _process_CTB_AUX(*cont_frag[ii]);
	    }
	}
    }

  art::Handle<artdaq::Fragments> frags;
  evt.getByLabel(fInputLabel, fInputNonContainerInstance, frags); 

  if(frags.isValid())
    {
      for(auto const& frag: *frags)
	{
	  _process_CTB_AUX(frag);
	}
    }

  raw::ctb::pdspctb ctbdp(fTrigs,fChStats,fFeedbacks,fMiscs,fWordIndexes);
  pdspctbs.push_back(ctbdp);  // just one for now
  evt.put(std::make_unique<std::vector<raw::ctb::pdspctb>>(std::move(pdspctbs)),fOutputLabel);

}

void PDSPCTBRawDecoder::_process_CTB_AUX(const artdaq::Fragment& frag)
{
  dune::CTBFragment ctbfrag(frag);

  // use the same logic in dune-raw-data/Overlays/CTBFragment.cc:operator<<

  for (size_t iword = 0; iword < ctbfrag.NWords(); ++iword)
    {
      size_t ix=0;
      uint32_t wt=0;
      if (ctbfrag.Trigger(iword))
	{
	  raw::ctb::Trigger tstruct;
	  tstruct.word_type = ctbfrag.Trigger(iword)->word_type;
	  wt = tstruct.word_type;
	  tstruct.trigger_word = ctbfrag.Trigger(iword)->trigger_word;
	  tstruct.timestamp = ctbfrag.Trigger(iword)->timestamp;
	  ix = fTrigs.size();
	  fTrigs.push_back(tstruct);
	}
      else if (ctbfrag.ChStatus(iword))
	{
	  raw::ctb::ChStatus cstruct;
	  cstruct.word_type = ctbfrag.ChStatus(iword)->word_type;
	  wt = cstruct.word_type;
	  cstruct.pds = ctbfrag.ChStatus(iword)->pds;
	  cstruct.crt = ctbfrag.ChStatus(iword)->crt;
	  cstruct.beam_hi = ctbfrag.ChStatus(iword)->beam_hi;
	  cstruct.beam_lo = ctbfrag.ChStatus(iword)->beam_lo;
	  cstruct.timestamp = ctbfrag.ChStatus(iword)->timestamp;
	  ix = fChStats.size();
	  fChStats.push_back(cstruct);
	}
      else if (ctbfrag.Feedback(iword))
	{
	  raw::ctb::Feedback fstruct;
	  fstruct.word_type = ctbfrag.Feedback(iword)->word_type;
	  wt = fstruct.word_type;
	  fstruct.padding = ctbfrag.Feedback(iword)->padding;
	  fstruct.source = ctbfrag.Feedback(iword)->source;
	  fstruct.code = ctbfrag.Feedback(iword)->code;
	  fstruct.timestamp = ctbfrag.Feedback(iword)->timestamp;
	  ix = fFeedbacks.size();
	  fFeedbacks.push_back(fstruct);
	}
      else
	{
	  raw::ctb::Misc mstruct;
	  mstruct.word_type = ctbfrag.Word(iword)->word_type;
	  wt = mstruct.word_type;
	  mstruct.payload = ctbfrag.Word(iword)->payload;
	  mstruct.timestamp = ctbfrag.Word(iword)->timestamp;
	  ix = fMiscs.size();
	  fMiscs.push_back(mstruct);
	}

      raw::ctb::WordIndex wstruct;
      wstruct.word_type = wt;
      wstruct.index = ix;
      fWordIndexes.push_back(wstruct);

    }
}


DEFINE_ART_MODULE(PDSPCTBRawDecoder)
