// TriggerTypes.h
//
// Michael Baird
// March 2016
//
// Provides a simple set of unsigned ints defining the various trigger types in one file.

#ifndef TRIGGERTYPES_H
#define TRIGGERTYPES_H

namespace triggersim {

  // Trigger types:
  const unsigned int kNullTrigger = 0;   // a default value for uninitialized variables
  const unsigned int kActivity    = 1;   // a simple trigger based off of some minimal activity above threshold
  const unsigned int kSNeShort    = 2;   // supernova trigger generated from the quick neutrino burst
  const unsigned int kSNeLong     = 3;   // supernova trigger generated from the full time window

} // end namespace triggersim

#endif
