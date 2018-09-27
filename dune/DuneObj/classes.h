// Build a dictionary.
//
// $Author:  $
// $Date: 2010/04/12 18:12:28 $
//
// Original author Rob Kutschke, modified by klg

#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"

#include "dune/DuneObj/CalibTreeRecord.h"
#include "dune/DuneObj/OpDetDivRec.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardataobj/RawData/OpDetWaveform.h"
//#include <bitset>
#include "dune/DuneObj/ProtoDUNEBeamEvent.h"

#include <vector>
#include <map>
//#include <array>

//template class art::wrapper<art::Assns<raw::OpDetWaveform,            sim::OpDetBacktrackerRecord,   void>  >
