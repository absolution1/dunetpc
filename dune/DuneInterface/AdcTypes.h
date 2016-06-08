// AdcTypes.h

// David Adams
// December 2015
//
// Integer and floating type for Adc signals.

#ifndef AdcTypes_H
#define AdcTypes_H

#include <vector>

typedef short AdcCount;
typedef std::vector<AdcCount> AdcCountVector;

typedef float AdcSignal;
typedef std::vector<AdcSignal> AdcSignalVector;
typedef std::vector<AdcSignalVector> AdcSignalVectorVector;

typedef float AdcPedestal;

typedef std::vector<bool> AdcFilterVector;

typedef short AdcFlag;
typedef std::vector<AdcFlag> AdcFlagVector;

const AdcFlag AdcGood      = 0; // ADC sample is fine
const AdcFlag AdcUnderflow = 1; // ADC sample is underflow
const AdcFlag AdcOverflow  = 2; // ADC sample is underflow
const AdcFlag AdcStuckOff  = 3; // ADC sample has bits stuck at 0
const AdcFlag AdcStuckOn   = 4; // ADC sample has bits stuck at 1

typedef unsigned int AdcChannel;

#endif
