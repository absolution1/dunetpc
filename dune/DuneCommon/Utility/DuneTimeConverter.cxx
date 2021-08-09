// DuneTimeConverter.cxx

#include "DuneTimeConverter.h"
#include <sstream>
#include <iomanip>

using std::string;
using std::ostringstream;
using std::istringstream;
using std::setw;
using std::setfill;

//**********************************************************************

art::Timestamp DuneTimeConverter::fromNova(uint64_t novaTime) {
  const uint64_t nsecPerUsec = 1000;
  uint64_t secSinceNovaT0 = novaTime/novaTicksPerSec();
  uint64_t ticksRem = novaTime - novaTicksPerSec()*secSinceNovaT0;
  uint64_t sec = secSinceNovaT0 + novaT0Sec();
  uint64_t nsec = (ticksRem*nsecPerUsec)/novaTicksPerUsec();
  uint64_t tart = (sec << 32) + nsec;
  art::Timestamp ts(tart);
  return ts;
}

//**********************************************************************

uint64_t DuneTimeConverter::toNova(art::Timestamp tart) {
  uint64_t thi = tart.timeHigh();
  uint64_t tlo = tart.timeLow();
  uint64_t tnova = novaTicksPerSec()*(thi - novaT0Sec()) + (tlo*novaTicksPerUsec())/1000;
  return tnova;
}

//**********************************************************************

art::Timestamp DuneTimeConverter::makeTimestamp(uint32_t tsec, uint32_t tns) {
  uint64_t tthi = tsec;
  uint64_t thilo = (tthi << 32) + tns;
  return art::Timestamp(thilo);
}

//**********************************************************************

string DuneTimeConverter::toString(art::Timestamp tart) {
  ostringstream sstime;
  sstime << tart.timeHigh();
  sstime << ".";
  sstime << setw(9) << setfill('0') << tart.timeLow();
  return sstime.str();
}

//**********************************************************************

art::Timestamp DuneTimeConverter::fromString(string stime) {
  string::size_type ipos = stime.find(".");
  if ( ipos != string::npos ) ++ipos;
  while ( ipos < stime.size() && stime[ipos] == '0' ) ++ipos;
  bool hasns = ipos < stime.size();
  uint32_t tsec;
  istringstream sssec(stime.substr(0,ipos));
  sssec >> tsec;
  uint32_t tns = 0;
  if ( hasns ) {
    istringstream ssns(stime.substr(ipos));
    ssns >> tns;
  }
  return makeTimestamp(tsec, tns);
}

//**********************************************************************
