#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "art/Framework/Core/ProductRegistryHelper.h"

// From lardata
#include "RawData/RawDigit.h"

#include <vector>
#include <string>

using std::vector;
using std::string;

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
    vector<string>        filenames_;
    vector<raw::RawDigit> bufferedDigits_;
    size_t                nextDigit_;
  };
}

DAQToOffline::Splitter::Splitter(fhicl::ParameterSet const& ps,
				 art::ProductRegistryHelper&,
				 art::SourceHelper& ) :
  filenames_(ps.get<vector<string>>("fileNames")),
  bufferedDigits_(),
  nextDigit_()
{
}

bool
DAQToOffline::Splitter::readFile(string const&, art::FileBlock*& fb)
{
  fb = nullptr;
  return true;
}

bool
DAQToOffline::Splitter::readNext(art::RunPrincipal* const&,
				 art::SubRunPrincipal* const&,
				 art::RunPrincipal*&,
				 art::SubRunPrincipal*&,
				 art::EventPrincipal*&)
{
  return false;
}

void
DAQToOffline::Splitter::closeCurrentFile()
{
}

// To use this source, use the module_type
// DAQToOffline::SplitterInput.
namespace DAQToOffline
{
  using SplitterInput = art::Source<Splitter>;
}

DEFINE_ART_INPUT_SOURCE(DAQToOffline::SplitterInput)
