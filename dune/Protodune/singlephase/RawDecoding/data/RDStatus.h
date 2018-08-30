////////////////////////////////////////////////////////////////////////
//
// RDStatus.h
// Status flags for unpacking raw::RawDigits 
// Tom Junk, August 29, 2018
//
////////////////////////////////////////////////////////////////////////

#ifndef  RDStatus_H
#define  RDStatus_H

#include "RtypesCore.h"
#include <stdint.h>

namespace raw {

  class RDStatus
  {

  public:

  RDStatus() : fCorruptDataDropped(false), fCorruptDataKept(false), fStatWord(0) {}; // Default constructor

  RDStatus(bool corruptdatadropped,
	   bool corruptdatakept,
	   unsigned int statword) : 
    fCorruptDataDropped(corruptdatadropped), 
    fCorruptDataKept(corruptdatakept),
    fStatWord(statword) {};

    bool GetCorruptDataDroppedFlag() { return fCorruptDataDropped; };
    bool GetCorruptDataKeptFlag() { return fCorruptDataKept; };
    unsigned int GetStatWord() { return fStatWord; };

  private:

    bool fCorruptDataDropped;  ///< True if some data suspected to be corrupt were not decoded
    bool fCorruptDataKept;     ///< True if some data suspected to be corrupt were decoded
    unsigned int fStatWord; 

  };


} // namespace raw

#endif // RDStatus_H
