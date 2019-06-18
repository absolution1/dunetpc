#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/Source.h"	
#include "dune/Protodune/dualphase/RawDecoding/PDDPRawInputDriver.h"

namespace lris 
{
  typedef art::Source<PDDPRawInputDriver> PDDPRawInput;	
}

DEFINE_ART_INPUT_SOURCE(lris::PDDPRawInput)
