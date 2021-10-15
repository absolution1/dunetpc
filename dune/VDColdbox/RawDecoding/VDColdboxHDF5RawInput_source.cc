#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "dune/VDColdbox/RawDecoding/VDColdboxHDF5RawInput.h"

// Use the standard service interfaces (CatalogInterface and
// FileTransfer) to obtain files.
//namespace art {
//template <raw::VDColdboxHDF5RawInputDetail>
//struct Source_wantFileServices {
//  static constexpr bool value = false
//};
//}

//typedef for shorthand
namespace raw {
using VDColdboxHDF5RawInputSource = art::Source<VDColdboxHDF5RawInputDetail>;
}

DEFINE_ART_INPUT_SOURCE(raw::VDColdboxHDF5RawInputSource)
