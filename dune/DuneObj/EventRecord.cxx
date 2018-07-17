//includes
#include "EventRecord.h"
namespace EventRecord {

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 
  void EventRecord::Clear() 
      {
        eves.clear();
        hits.clear();
        ophits.clear();
        run=0;
        subrun=0;
        event=0;
      }

}

