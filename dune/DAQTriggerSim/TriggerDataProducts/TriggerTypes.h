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
  const unsigned int kBeamNu      = 4;   // neutrino beam event
  const unsigned int kProtonDecay = 5;   // proton decay event
  const unsigned int kNnbar       = 6;   // n,n-bar oscillations
  const unsigned int kMonopole    = 7;   // magnetic monopole
  const unsigned int kAtmoNu      = 8;   // atmospheric neutrino
  const unsigned int kSolarNu     = 9;   // solar neutrino
  const unsigned int kCosmic      = 10;  // cosmic ray primary or secondary

} // end namespace triggersim

#endif
