#include "dune/DuneObj/ProtoDUNETimeStamp.h"

dune::ProtoDUNETimeStamp::ProtoDUNETimeStamp()
    :  m_cookie(0),
       m_triggerType((dune::ProtoDUNETimingCommand)0),
       m_reservedBits(0),
       m_timeStamp(0),
       m_eventCounter(0),
       m_checksumGood(0),
       m_lastRunStart(0),
       m_lastSpillStart(0),
       m_lastSpillEnd(0)
{}
