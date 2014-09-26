#ifndef lbne_artdaq_Overlays_FragmentType_hh
#define lbne_artdaq_Overlays_FragmentType_hh
#include "artdaq-core/Data/Fragment.hh"

namespace lbne {

  namespace detail {
    enum FragmentType : artdaq::Fragment::type_t
    { MISSED = artdaq::Fragment::FirstUserFragmentType,
        TPC,
        PHOTON,
        TRIGGER,
        TOY1,
        TOY2,
        INVALID // Should always be last.
        };

    // Safety check.
    static_assert(artdaq::Fragment::isUserFragmentType(FragmentType::INVALID - 1),
                  "Too many user-defined fragments!");
  }

  using detail::FragmentType;

  FragmentType toFragmentType(std::string t_string);
  std::string fragmentTypeToString(FragmentType val);
}
#endif /* lbne_artdaq_Overlays_FragmentType_hh */
