// art
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Root/rootNames.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Persistency/Provenance/EventID.h"
#include "art/Persistency/Provenance/RunID.h"
#include "art/Persistency/Provenance/SubRunID.h"
#include "art/Persistency/Provenance/MasterProductRegistry.h"
#include "art/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// artdaq
#include "artdaq-core/Data/Fragments.hh"

// lardata
#include "RawData/RawDigit.h"

// lbne
#include "lbne/daqinput35t/tpcFragmentToRawDigits.h"

// C++
#include <functional>
#include <memory>
#include <regex>
#include <vector>
#include <string>

// ROOT
#include "TTree.h"
#include "TFile.h"

using raw::RawDigit;
using std::vector;
using std::string;

//==========================================================================

namespace {


  struct LoadedDigits {

    LoadedDigits() : digits(), index(0) {}

    vector<RawDigit> digits;
    size_t index;

    void assign(vector<RawDigit> const & v ) {
      digits = v;
      index = 0ul;
    }

    bool empty() const { return index == digits.size(); }

    RawDigit next() {
      ++index;
      return digits.at(index-1);
    }

  };

  // Retrieves branch name (a la art convention) where object resides
  template <typename PROD>
  const char* getBranchName( art::InputTag const & tag )
  {
    std::ostringstream pat_s;
    pat_s << art::TypeID(typeid(PROD)).friendlyClassName()
          << '_'
          << tag.label()
          << '_'
          << tag.instance()
          << '_'
          << tag.process()
          << ".obj";
    return pat_s.str().data();
  }

  artdaq::Fragments*
  getFragments( TBranch* br, unsigned entry )
  {
    br->GetEntry( entry );
    return reinterpret_cast<artdaq::Fragments*>( br->GetAddress() );
  }

  // Assumed file format is
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

    string                 sourceName_;
    string                 lastFileName_;
    std::unique_ptr<TFile> file_;
    bool                   doneWithFiles_;
    art::InputTag          inputTag_;
    art::SourceHelper      sh_;
    TBranch*               fragmentsBranch_;
    LoadedDigits           loadedDigits_;
    size_t                 digitIndex_;
    size_t                 nInputEvts_;
    size_t                 treeIndex_;
    art::RunNumber_t       runNumber_;
    art::SubRunNumber_t    subRunNumber_;
    art::EventNumber_t     eventNumber_;
    art::RunNumber_t       cachedRunNumber_;
    art::SubRunNumber_t    cachedSubRunNumber_;
    size_t                 bufferLimit_;
    rawDigits_t            bufferedDigits_;

    std::function<rawDigits_t(artdaq::Fragments const&)> fragmentsToDigits_;

    bool eventIsFull_(rawDigits_t const & v);

    bool loadDigits_();

    void makeEventAndPutDigits_( art::EventPrincipal*& outE );
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
  inputTag_("daq:TPC:DAQ"), // "moduleLabel:instance:processName"
  sh_(sh),
  fragmentsBranch_(nullptr),
  nInputEvts_(),
  runNumber_(1),      // Defaults in case input filename does not
  subRunNumber_(0),   // follow assumed filename_format above.
  eventNumber_(),
  cachedRunNumber_(-1),
  cachedSubRunNumber_(-1),
  bufferLimit_(ps.get<size_t>("rawDigitsPerEvent")),
  bufferedDigits_(),
  fragmentsToDigits_( std::bind( DAQToOffline::tpcFragmentToRawDigits,
                                 std::placeholders::_1, // artdaq::Fragments
                                 ps.get<bool>("debug",false),
                                 ps.get<raw::Compress_t>("compression",raw::kNone),
                                 ps.get<unsigned>("zeroThreshold",0) ) )
{
  // Will use same instance name for the outgoing products as for the
  // incoming ones.
  prh.reconstitutes<rawDigits_t,art::InEvent>( sourceName_, inputTag_.instance() );
}

//=======================================================================================
bool
DAQToOffline::Splitter::readFile(string const& filename, art::FileBlock*& fb)
{

  // Run numbers determined based on file name...see comment in
  // anon. namespace above.
  std::smatch matches;
  if ( std::regex_match( filename, matches, filename_format ) ) {
    runNumber_       = std::stoul( matches[1] );
    subRunNumber_    = std::stoul( matches[2] );
  }

  // Get fragments branch
  file_.reset( new TFile(filename.data()) );
  TTree* evtree    = reinterpret_cast<TTree*>(file_->Get(art::rootNames::eventTreeName().c_str()));
  fragmentsBranch_ = evtree->GetBranch( getBranchName<artdaq::Fragments>( inputTag_ ) ); // get branch for specific input tag
  nInputEvts_      = static_cast<size_t>( fragmentsBranch_->GetEntries() );
  treeIndex_       = 0u;

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

  art::Timestamp ts; // LBNE should decide how to initialize this
  if ( runNumber_ != cachedRunNumber_ ){
    outR = sh_.makeRunPrincipal(runNumber_,ts);
    cachedRunNumber_ = runNumber_;
  }

  if ( subRunNumber_ != cachedSubRunNumber_ ) {
    outSR = sh_.makeSubRunPrincipal(runNumber_,subRunNumber_,ts);
    cachedSubRunNumber_ = subRunNumber_;
  }

  // eventIsFull_ is what LBNE should modify based on its needs
  while ( !eventIsFull_( bufferedDigits_ ) ) {

    if ( loadedDigits_.empty() && !loadDigits_() ) {
      if ( file_->GetName() != lastFileName_ )
        return false;
      else {
        doneWithFiles_ = true;
        break;
      }
    }

    bufferedDigits_.emplace_back( loadedDigits_.next() );

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
  return v.size() == bufferLimit_;
}

//=======================================================================================
bool
DAQToOffline::Splitter::loadDigits_()
{
  if ( loadedDigits_.empty() && treeIndex_ != nInputEvts_ ) {
    auto* fragments = getFragments( fragmentsBranch_, treeIndex_++ );
    rawDigits_t const digits = fragmentsToDigits_( *fragments );
    loadedDigits_.assign( digits );
    return true;
  }
  else return false;
}

//=======================================================================================
void
DAQToOffline::Splitter::makeEventAndPutDigits_(art::EventPrincipal*& outE){
  ++eventNumber_;
  outE = sh_.makeEventPrincipal( runNumber_, subRunNumber_, eventNumber_, art::Timestamp() );
  art::put_product_in_principal( std::make_unique<rawDigits_t>(bufferedDigits_),
                                 *outE,
                                 sourceName_,
                                 inputTag_.instance() );
  mf::LogDebug("DigitsTest") << "Producing event: " << outE->id() << " with " << bufferedDigits_.size() << " digits";
  bufferedDigits_.clear();
}

//=======================================================================================
DEFINE_ART_INPUT_SOURCE(art::Source<DAQToOffline::Splitter>)
//=======================================================================================
