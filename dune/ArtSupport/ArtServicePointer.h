// ArtServicePointer.h

#ifndef ArtServicePointer_H
#define ArtServicePointer_H

// David Adams
// January 2017
//
// ArtServicePointer provides access to an art service.
//
// Usage:
//   MyService* psvc = ArtServicePointer<MyService>()    
//
// Template specializations may be used enable service without
// direct use of ServiceHandle, i.e. inside root.

// This is the function that is specialized and used by clients.
template<class T>
T* ArtServicePointer();

#ifdef __CLING__

#include <iostream>
template<class T>
T* ArtServicePointer() {
  std::cout << "ArtServicePointer: No specialization found for " << typeid(T).name() << std::endl;
  return nullptr;
}

#else

#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"

// This should not be called directly.
template<class T>
T* GenericArtServicePointer() {
  try {
    art::ServiceRegistry& reg = art::ServiceRegistry::instance();
    if ( ! reg.isAvailable<T>() ) return nullptr;
    return &reg.get<T>();
  } catch(...) {
    std::cout << "WARNING: ArtServicePointer: ServiceRegistry not found. Instantiate with" << std::endl;
    std::cout << "           ArtServiceHelper::instance(\"myfcl.fcl\")" << std::endl;
  }
  return nullptr;
}

template<class T>
T* ArtServicePointer() {
  return GenericArtServicePointer<T>();
}

#endif

#endif
