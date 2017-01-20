// ArtServiceProvider.h

#ifndef ArtServiceProvider_H
#define ArtServiceProvider_H

// David Adams
// January 2017
//
// ArtServiceProvider provides access to an art service provider.
//
// Usage:
//   MyServiceProvider* pprv = ArtServiceProvider<MyService>()    
//
// Template specializations may be used enable service without
// direct use of ServiceHandle, i.e. inside root.

#include "dune/ArtSupport/ArtServicePointer.h"

// Access a provider via service and provider types.
template<class S, class P>
const P* ArtServiceProviderFromService() {
  S* psvc = ArtServicePointer<S>();
  if ( psvc == nullptr ) return nullptr;
  return psvc->GetProviderPtr();
}

// Access a provide via service type assuming the service
// knows its provider type.
template<class S>
const typename S::provider_type* DefaultArtServiceProvider() {
  return ArtServiceProviderFromService<S, typename S::provider_type>();
}

// Access a provider via provider type.
// This requires specialization.
template<class P>
const P* ArtServiceProvider() {
  std::cout << "ArtServiceProvider: No specialization found for " << typeid(P).name() << std::endl;
  return nullptr;
}

#endif
