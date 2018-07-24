//includes
#include "CalibTreeRecord.h"
namespace CalibTreeRecord {

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= 
  void CalibTreeRecord::Clear() 
      {
        eves.clear();
        hits.clear();
        ophits.clear();
        run=0;
        subrun=0;
        event=0;
      }

}

