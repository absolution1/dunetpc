#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/Source.h"	
#include "dune/3x1x1dp/DataImport/Services/RawData311InputDriver.h"

namespace lris 
{
  typedef art::Source<RawData311InputDriver> ImportFull311File;	
}

DEFINE_ART_INPUT_SOURCE(lris::ImportFull311File)
