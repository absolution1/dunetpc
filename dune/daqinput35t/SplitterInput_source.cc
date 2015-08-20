// art
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Root/rootNames.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Persistency/Provenance/EventID.h"
#include "art/Persistency/Provenance/MasterProductRegistry.h"
#include "art/Persistency/Provenance/RunID.h"
#include "art/Persistency/Provenance/SubRunID.h"
#include "art/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Provenance/EventAuxiliary.h"

//#include <iostream>

// artdaq
#include "artdaq-core/Data/Fragments.hh"

// lardata
#include "RawData/RawDigit.h"

// dune
#include "tpcFragmentToRawDigits.h"
#include "SSPFragmentToOpDetWaveform.h"

// C++
#include <functional>
#include <memory>
#include <regex>
#include <string>
#include <vector>

// ROOT
#include "TTree.h"
#include "TFile.h"

using raw::RawDigit;
using raw::OpDetWaveform;
using std::vector;
using std::string;

// TODO:
//  Split the SSP data too.
//  Change the index in the loadeddigits based on where we want to draw the split
//  Put in external trigger (Penn board) info when we get it
//  Deal with ZS data and MC -- currently this assumes non-ZS data
//  Discover if an event is not contiguous with the next event and do not stitch in that case -- dump the
//    loaded digits and the partially constructed event and start over
//  There's an assumption that the channels are the same event to event and they come in the same order --
//    would need a list of channels from one event to another to look it up properly if we are missing some channels

//==========================================================================

namespace {

  struct LoadedWaveforms {

    LoadedWaveforms() : waveforms() {}
    vector<OpDetWaveform> waveforms;

    void load( vector<OpDetWaveform> const & v ) {
      waveforms = v;
    }

    // not really used 
    bool empty() const 
    { 
      if (waveforms.size() == 0) return true;
      return false;
    }

    // retrieve bits of opdetwaveforms between tpc index beginnings and ends, and append
    // them to the currently active ones if need be, and make new ones if not.

    void append( vector<OpDetWaveform> &v, std::set<unsigned int> &curwforms, 
                 unsigned int tpcindexbeg, unsigned int tpcindexend )
    {
      for (auto wf : waveforms) // see if any waveforms have pieces inside this TPC boundary
	{
	  
	}
    }

  };


  struct LoadedDigits {

    LoadedDigits() : digits(), index(), firstTimestamp(0) {}

    vector<RawDigit> digits;
    size_t index;
    lbne::TpcNanoSlice::Header::nova_timestamp_t firstTimestamp; // timestamp of the first nanoslice in the digits vector

    void loadTimestamp(lbne::TpcNanoSlice::Header::nova_timestamp_t ts)
    {
      firstTimestamp = ts;
    }

    // copy rawdigits to LoadedDigits and uncompress if necessary.  

    void load( vector<RawDigit> const & v ) {
      if (v.size() == 0 || v.back().Compression() == raw::kNone)
	{ 
	  digits = v; 
	}
      else
	{
	  // make a new raw::RawDigit object for each compressed one
	  // to think about optimization -- two copies of the uncompressed raw digits here.
	  digits = std::vector<RawDigit>();
          for (auto idigit = v.begin(); idigit != v.end(); ++idigit)
	    {
	       std::vector<short> uncompressed;
	       raw::Uncompress(idigit->ADCs(),uncompressed,idigit->Compression());
	       raw::RawDigit digit(idigit->Channel(),
				   idigit->Samples(),
				   uncompressed,
				   raw::kNone);
	       digits.push_back(digit);
	    }
	}
      index = 0ul;
    }

    bool empty() const 
    { 
      //std::cout << "Checking if loaded digits is empty" << std::endl;
      //std::cout << "Number of rawdigits: " << digits.size() << std::endl;
      if (digits.size() == 0) return true;
      //std::cout << "Number of samples: " << digits[0].Samples() << std::endl;
      if (digits[0].Samples() == 0) return true;
      if (index >= digits[0].Samples()) return true; 
      return false;
    }

    RawDigit::ADCvector_t next() {
      ++index;
      RawDigit::ADCvector_t digitsatindex;
      for (size_t ichan=0; ichan<digits.size(); ichan++)
	{
	  digitsatindex.push_back(digits[ichan].ADC(index-1));
	}
      return digitsatindex;
    }

  };

  // Retrieves branch name (a la art convention) where object resides
  const char* getBranchName( art::InputTag const & tag, const string inputDataProduct )
  {
    std::ostringstream pat_s;
    pat_s << inputDataProduct << "s" 
          << '_'
          << tag.label()
          << '_'
          << tag.instance()
          << '_'
          << tag.process()
          << ".obj";

    //std::cout << inputDataProduct << "s" 
    //    << '_'
    //    << tag.label()
    //    << '_'
    //    << tag.instance()
    //      << '_'
    //      << tag.process()
    //      << ".obj" << std::endl;

    return pat_s.str().data();
  }

  artdaq::Fragments*
  getFragments( TBranch* br, unsigned entry )
  {
    br->GetEntry( entry );
    return reinterpret_cast<artdaq::Fragments*>( br->GetAddress() );
  }

  vector<raw::RawDigit>*
  getRawDigits( TBranch* br, unsigned entry )
  {
    br->GetEntry( entry );
    return reinterpret_cast<vector<raw::RawDigit>*>( br->GetAddress() );
  }

  vector<raw::OpDetWaveform>*
  getSSPWaveforms( TBranch* br, unsigned entry )
  {
    br->GetEntry( entry );
    return reinterpret_cast<vector<raw::OpDetWaveform>*>( br->GetAddress() );
  }


  // Assumed file format is no longer: we get run and subrun out of evAux now
  //
  //     "lbne_r[digits]_sr[digits]_other_stuff.root"
  //
  // This regex object is used to extract the run and subrun numbers.
  // The '()' groupings correspond to regex matches that can be
  // extracted using the std::regex_match facility.
  std::regex const filename_format( "lbne_r(\\d+)_sr(\\d+).*\\.root" );

}

//==========================================================================

namespace raw {

  // Enable 'pset.get<raw::Compress_t>("compression")'
  void decode (boost::any const & a, Compress_t & result ){
    unsigned tmp;
    fhicl::detail::decode(a,tmp);
    result = static_cast<Compress_t>(tmp);
  }
}

//==========================================================================

namespace DAQToOffline
{
  // The class Splitter is to be used as the template parameter for
  // art::Source<T>. It understands how to read art's ROOT data files,
  // to extract artdaq::Fragments from them, how to convert the
  // artdaq::Fragments into vectors of raw::RawDigit objects, and
  // finally how to re-combine those raw::RawDigit objects to present
  // the user of the art::Source<Splitter> with a sequence of
  // art::Events that have the desired event structure, different from
  // the event structure of the data file(s) being read.

  class Splitter
  {
  public:
    Splitter(fhicl::ParameterSet const& ps,
             art::ProductRegistryHelper& prh,
             art::SourceHelper& sh);

    // See art/Framework/IO/Sources/Source.h for a description of each
    // of the public member functions of Splitter.
    bool readFile(string const& filename, art::FileBlock*& fb);

    bool readNext(art::RunPrincipal* const& inR,
                  art::SubRunPrincipal* const& inSR,
                  art::RunPrincipal*& outR,
                  art::SubRunPrincipal*& outSR,
                  art::EventPrincipal*& outE);

    void closeCurrentFile();

  private:

    using rawDigits_t = vector<RawDigit>;
    using SSPWaveforms_t = vector<OpDetWaveform>;

    string                 sourceName_;
    string                 lastFileName_;
    std::unique_ptr<TFile> file_;
    bool                   doneWithFiles_;
    art::InputTag          TPCinputTag_;
    art::InputTag          SSPinputTag_;
    string                 TPCinputDataProduct_;
    string                 SSPinputDataProduct_;
    double                 fNOvAClockFrequency; // MHz
    art::SourceHelper      sh_;
    TBranch*               TPCinputBranch_;
    TBranch*               SSPinputBranch_;
    TBranch*               EventAuxBranch_;
    LoadedDigits           loadedDigits_;
    LoadedWaveforms        loadedWaveforms_;
    size_t                 nInputEvts_;
    size_t                 treeIndex_;
    art::RunNumber_t       runNumber_;
    art::SubRunNumber_t    subRunNumber_;
    art::EventNumber_t     eventNumber_;
    art::RunNumber_t       cachedRunNumber_;
    art::SubRunNumber_t    cachedSubRunNumber_;
    art::RunNumber_t       inputRunNumber_;
    art::SubRunNumber_t    inputSubRunNumber_;
    art::EventNumber_t     inputEventNumber_;
    size_t                 ticksPerEvent_;
    rawDigits_t            bufferedDigits_;
    std::vector<RawDigit::ADCvector_t>  dbuf_;
    unsigned short         fTicksAccumulated;
    std::vector<raw::OpDetWaveform> sbuf_;
    std::set<unsigned int> curwaveforms_; // list of waveform indices that are currently being built up

    std::function<rawDigits_t(artdaq::Fragments const&, lbne::TpcNanoSlice::Header::nova_timestamp_t& )> fragmentsToDigits_;

    bool eventIsFull_(rawDigits_t const & v);

    bool loadDigits_();

    void makeEventAndPutDigits_( art::EventPrincipal*& outE );

    art::EventAuxiliary    evAux_;
    art::EventAuxiliary*   pevaux_;
  };
}

//=======================================================================================
DAQToOffline::Splitter::Splitter(fhicl::ParameterSet const& ps,
                                 art::ProductRegistryHelper& prh,
                                 art::SourceHelper& sh) :
  sourceName_("SplitterInput"),
  lastFileName_(ps.get<vector<string>>("fileNames",{}).back()),
  file_(),
  doneWithFiles_(false),
  //  TPCinputTag_("daq:TPC:DAQ"), // "moduleLabel:instance:processName"
  TPCinputTag_(ps.get<string>("TPCInputTag")), // "moduleLabel:instance:processName"
  SSPinputTag_(ps.get<string>("SSPInputTag")), // "moduleLabel:instance:processName"
  TPCinputDataProduct_(ps.get<string>("TPCInputDataProduct")),
  SSPinputDataProduct_(ps.get<string>("SSPInputDataProduct")),
  fNOvAClockFrequency(ps.get<double>("NOvAClockFrequency",64.0)),
  sh_(sh),
  TPCinputBranch_(nullptr),
  SSPinputBranch_(nullptr),
  EventAuxBranch_(nullptr),
  nInputEvts_(),
  runNumber_(1),      // Defaults in case input filename does not
  subRunNumber_(0),   // follow assumed filename_format above.
  eventNumber_(),
  cachedRunNumber_(-1),
  cachedSubRunNumber_(-1),
  inputRunNumber_(1),      
  inputSubRunNumber_(1),   
  inputEventNumber_(1),

  ticksPerEvent_(ps.get<size_t>("ticksPerEvent")),
  bufferedDigits_(),
  dbuf_(),
  fTicksAccumulated(0),
  fragmentsToDigits_( std::bind( DAQToOffline::tpcFragmentToRawDigits,
                                 std::placeholders::_1, // artdaq::Fragments
                                 std::placeholders::_2, // lbne::TpcNanoSlice::Header::nova_timestamp_t& firstTimestamp
                                 ps.get<bool>("debug",false),
                                 ps.get<raw::Compress_t>("compression",raw::kNone),
                                 ps.get<unsigned>("zeroThreshold",0) ) )
{
  // Will use same instance names for the outgoing products as for the
  // incoming ones.
  prh.reconstitutes<rawDigits_t,art::InEvent>( sourceName_, TPCinputTag_.instance() );
  prh.reconstitutes<SSPWaveforms_t,art::InEvent>( sourceName_, SSPinputTag_.instance() );
}

//=======================================================================================
bool
DAQToOffline::Splitter::readFile(string const& filename, art::FileBlock*& fb)
{

  // Now use the run and subrun number input from the EventAuxiliary object instead of the file name
  // needed so we can split MC which has arbitrary filenames

  // Run numbers determined based on file name...see comment in
  // anon. namespace above.
  //std::smatch matches;
  //if ( std::regex_match( filename, matches, filename_format ) ) {
  //  runNumber_       = std::stoul( matches[1] );
  //  subRunNumber_    = std::stoul( matches[2] );
  //}

  // Get fragments branches
  file_.reset( new TFile(filename.data()) );
  TTree* evtree    = reinterpret_cast<TTree*>(file_->Get(art::rootNames::eventTreeName().c_str()));

  TPCinputBranch_ = evtree->GetBranch( getBranchName(TPCinputTag_, TPCinputDataProduct_ ) ); // get branch for TPC input tag
  SSPinputBranch_ = evtree->GetBranch( getBranchName(SSPinputTag_, SSPinputDataProduct_ ) ); // get branch for SSP input tag

  nInputEvts_      = static_cast<size_t>( TPCinputBranch_->GetEntries() );
  size_t nevt_ssp  = static_cast<size_t>( SSPinputBranch_->GetEntries() );
  if (nevt_ssp != nInputEvts_) throw cet::exception("35-ton SplitterInput: Different numbers of RCE and SSP input events in file");
  treeIndex_       = 0ul;

  EventAuxBranch_ = evtree->GetBranch( "EventAuxiliary" );
  pevaux_ = &evAux_;
  EventAuxBranch_->SetAddress(&pevaux_);

  // New fileblock
  fb = new art::FileBlock(art::FileFormatVersion(),filename);
  if ( fb == nullptr ) {
    throw art::Exception(art::errors::FileOpenError)
      << "Unable to open file " << filename << ".\n";
  }

  return true;
}

//=======================================================================================
bool
DAQToOffline::Splitter::readNext(art::RunPrincipal*    const& inR,
                                 art::SubRunPrincipal* const& inSR,
                                 art::RunPrincipal*    & outR,
                                 art::SubRunPrincipal* & outSR,
                                 art::EventPrincipal*  & outE)
{

  if ( doneWithFiles_ ) {
    return false;
  }


  while ( fTicksAccumulated < ticksPerEvent_ ) {


    while (loadedDigits_.empty())
      {
	bool rc = loadDigits_();
	if (!rc)
	  {
	    doneWithFiles_ = (file_->GetName() == lastFileName_);
	    return false;
	  }
      }

    std::vector<short> nextdigits = loadedDigits_.next(); 
    if (dbuf_.size() == 0)
      {
	RawDigit::ADCvector_t emptyvector; 
	for (size_t ichan=0;ichan<nextdigits.size();ichan++) dbuf_.push_back(emptyvector);
      }
    for (size_t ichan=0;ichan<nextdigits.size();ichan++)
      {
	dbuf_[ichan].push_back(nextdigits[ichan]);
      }
    fTicksAccumulated ++;
  }

  // copy all the info (channel number, digits, pedestal, sigma, compression) from the current loaded digits

  for (size_t ichan=0;ichan<dbuf_.size();ichan++)
    {
      RawDigit d(loadedDigits_.digits[ichan].Channel(),
                 fTicksAccumulated,dbuf_[ichan],
                 loadedDigits_.digits[ichan].Compression());
      d.SetPedestal(loadedDigits_.digits[ichan].GetPedestal(),
                    loadedDigits_.digits[ichan].GetSigma());
      bufferedDigits_.emplace_back(d);
    }

  runNumber_ = inputRunNumber_;
  subRunNumber_ = inputSubRunNumber_;

  art::Timestamp ts; // LBNE should decide how to initialize this
  if ( runNumber_ != cachedRunNumber_ ){
    outR = sh_.makeRunPrincipal(runNumber_,ts);
    cachedRunNumber_ = runNumber_;
    eventNumber_ = 0ul;
  }

  if ( subRunNumber_ != cachedSubRunNumber_ ) {
    outSR = sh_.makeSubRunPrincipal(runNumber_,subRunNumber_,ts);
    cachedSubRunNumber_ = subRunNumber_;
    eventNumber_ = 0ul;
  }

  makeEventAndPutDigits_( outE );
  return true;

}

//=======================================================================================
void
DAQToOffline::Splitter::closeCurrentFile()
{
  file_.reset(nullptr);
}

//=======================================================================================
bool
DAQToOffline::Splitter::eventIsFull_( vector<RawDigit> const & v )
{
  return v.size() == ticksPerEvent_;
}

//=======================================================================================
bool
DAQToOffline::Splitter::loadDigits_()
{
  if ( loadedDigits_.empty() && treeIndex_ != nInputEvts_ ) {

    EventAuxBranch_->GetEntry(treeIndex_);
    inputRunNumber_ = evAux_.run();
    inputSubRunNumber_ = evAux_.subRun();
    inputEventNumber_ = evAux_.event();

    if (TPCinputDataProduct_.find("Fragment") != std::string::npos)
      {
        lbne::TpcNanoSlice::Header::nova_timestamp_t firstTimestamp;
        auto* fragments = getFragments( TPCinputBranch_, treeIndex_ );
        rawDigits_t const digits = fragmentsToDigits_( *fragments, firstTimestamp );
        loadedDigits_.load( digits );
	loadedDigits_.loadTimestamp( firstTimestamp );
      }
    else
      {
	auto* digits = getRawDigits(TPCinputBranch_, treeIndex_ );
	loadedDigits_.load( *digits);
	loadedDigits_.loadTimestamp(0); // MC timestamp is zero (? assume?)

      }

    if (SSPinputDataProduct_.find("Fragment") != std::string::npos)
      {
        auto* SSPfragments = getFragments( SSPinputBranch_, treeIndex_ );
        std::vector<raw::OpDetWaveform> waveforms = DAQToOffline::SSPFragmentToOpDetWaveform(*SSPfragments, fNOvAClockFrequency);
        loadedWaveforms_.load( waveforms );
      }
    else
      {
	auto* waveforms = getSSPWaveforms(SSPinputBranch_, treeIndex_ );
        loadedWaveforms_.load( *waveforms );
      }

    treeIndex_++;

    return true;
  }
  else return false;
}

//=======================================================================================
void
DAQToOffline::Splitter::makeEventAndPutDigits_(art::EventPrincipal*& outE){

  // just keep incrementing the event number as we split along

  ++eventNumber_;

  outE = sh_.makeEventPrincipal( runNumber_, subRunNumber_, eventNumber_, art::Timestamp() );
  art::put_product_in_principal( std::make_unique<rawDigits_t>(bufferedDigits_),
                                 *outE,
                                 sourceName_,
                                 TPCinputTag_.instance() );
  mf::LogDebug("DigitsTest") << "Producing event: " << outE->id() << " with " << bufferedDigits_.size() << " digits";
  bufferedDigits_.clear();
  for (size_t ichan=0;ichan<dbuf_.size();ichan++) { dbuf_[ichan].clear(); }
  dbuf_.clear();
  fTicksAccumulated = 0;
}

//=======================================================================================
DEFINE_ART_INPUT_SOURCE(art::Source<DAQToOffline::Splitter>)
//=======================================================================================
