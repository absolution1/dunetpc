// DuneTimeConverter.cxx

#include "DuneTimeConverter.h"

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
