#include "dune/RegCVN/func/RegPixelMap.h"
#include "dune/RegCVN/func/RegCVNResult.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Wrapper.h"

template class std::vector<float>;
template class std::vector<cvn::RegPixelMap>;
template class art::Ptr<cvn::RegPixelMap>;
template class std::vector<art::Ptr<cvn::RegPixelMap> >;
template class art::Wrapper< std::vector<cvn::RegPixelMap> >;


template class std::vector<cvn::RegCVNResult>;
template class art::Ptr<cvn::RegCVNResult>;
template class std::vector<art::Ptr<cvn::RegCVNResult> >;

template class art::Wrapper< std::vector<cvn::RegCVNResult> >;


