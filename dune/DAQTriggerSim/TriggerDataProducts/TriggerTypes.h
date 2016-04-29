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
  const unsigned int kSNe         = 2;   // supernova trigger
  const unsigned int kNu          = 3;   // trigger from some kind of neutrino interaction
  const unsigned int kNucleon     = 4;   // a nucleon based event
  const unsigned int kMonopole    = 5;   // magnetic monopole trigger
  const unsigned int kCosmic      = 6;   // cosmic ray primary/secondary or some other cosmic induced event



  // Offset value for specifying trigger sub-types:
  const unsigned int kSubTypeOffset = 1000;



  // Activity trigger sub-types:
  // (specify activity trigger sub-types here...)



  // Supernova sub-types:
  const unsigned int kSNeShort = kSNe * kSubTypeOffset + 1;   // quick core-collapse burst
  const unsigned int kSNeLong  = kSNe * kSubTypeOffset + 2;   // trigger from a long (~10 seconds) time window



  // Neutrino sub-types:
  const unsigned int kNuBeam  = kNu * kSubTypeOffset + 1;   // beam neutrino
  const unsigned int kNuAtmo  = kNu * kSubTypeOffset + 2;   // atmospheric neutrino
  const unsigned int kNuSolar = kNu * kSubTypeOffset + 3;   // solar neutrino



  // Nucleon sub-types:
  const unsigned int kProtonDecay = kNucleon * kSubTypeOffset + 1;   // proton decay
  const unsigned int kNNbar       = kNucleon * kSubTypeOffset + 2;   // n -> n-bar oscillations



  // Cosmic sub-types:
  // (specify cosmic trigger sub-types here...)


} // end namespace triggersim

#endif
