#ifndef VDColdboxHDF5RawInput_h
#define VDColdboxHDF5RawInput_h
#include "art/Framework/Core/InputSourceMacros.h" 
#include "art/Framework/IO/Sources/Source.h" 
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "fhiclcpp/ParameterSet.h"


//Is raw a good namespace?
namespace raw {
//Forward declare the class
class VDColdboxHDF5RawInputDetail;
}

class raw::VDColdboxHDF5RawInputDetail {
 public:
  VDColdboxHDF5RawInputDetail(fhicl::ParameterSet const &,
                              art::ProductRegistryHelper &,
                              art::SourceHelper const &) {};
  void readFile(std::string const & filename, art::FileBlock*& fb) {};
  bool readNext(art::RunPrincipal const* const inR,
                art::SubRunPrincipal const* const inSR,
                art::RunPrincipal*& outR,
                art::SubRunPrincipal*& outSR,
                art::EventPrincipal*& outE) {return true;};
  void closeCurrentFile() {};

};
#endif

/////typedef for shorthand
///namespace raw {
///using VDColdboxHDF5RawInputSource = art::Source<VDColdboxHDF5RawInputDetail>;
///}
///
///DEFINE_ART_INPUT_SOURCE(raw::VDColdboxHDF5RawInput)
