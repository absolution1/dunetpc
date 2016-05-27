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

#endif
