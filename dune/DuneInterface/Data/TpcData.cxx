// TpcData.cxx

#include "dune/DuneInterface/Data/TpcData.h"

//**********************************************************************

TpcData* TpcData::addTpcData(Name nam) {
  if ( nam == "" || nam == "." ) return nullptr;
  Name::size_type ipos = nam.rfind("/");
  if ( ipos == Name::npos ) return m_dat.count(nam) ? nullptr : &m_dat[nam];
  TpcData* pdat = getTpcData(nam.substr(0, ipos));
  if ( pdat == nullptr ) return nullptr;
  return pdat->addTpcData(nam.substr(ipos+1));
}

//**********************************************************************

TpcData::AdcData* TpcData::addAdcData(Name nam) {
  Name::size_type ipos = nam.rfind("/");
  if ( ipos == Name::npos ) return m_acm.count(nam) ? nullptr : &m_acm[nam];
  TpcData* pdat = getTpcData(nam.substr(0, ipos));
  if ( pdat == nullptr ) return nullptr;
  return pdat->addAdcData(nam.substr(ipos+1));
}

//**********************************************************************

TpcData::Tpc2dRoiVector* TpcData::addTpcRois(Name nam) {
  Name::size_type ipos = nam.rfind("/");
  if ( ipos == Name::npos ) return m_r2d.count(nam) ? nullptr : &m_r2d[nam];
  TpcData* pdat = getTpcData(nam.substr(0, ipos));
  if ( pdat == nullptr ) return nullptr;
  return pdat->addTpcRois(nam.substr(ipos+1));
}

//**********************************************************************

TpcData* TpcData::getTpcData(Name nam) {
  if ( nam == "" || nam == "." ) return this;
  Name::size_type ipos = nam.find("/");
  if ( ipos == Name::npos ) return m_dat.count(nam) ? &m_dat[nam] : nullptr;
  TpcData* pdat = getTpcData(nam.substr(0,ipos));
  if ( pdat == nullptr ) return nullptr;
  return pdat->getTpcData(nam.substr(ipos+1));
}

//**********************************************************************

TpcData::AdcData* TpcData::getAdcData(Name nam) {
  if ( nam == "" || nam == "." ) return nullptr;
  Name::size_type ipos = nam.rfind("/");
  if ( ipos == Name::npos ) return m_acm.count(nam) ? &m_acm[nam] : nullptr;
  TpcData* pdat = getTpcData(nam.substr(0,ipos));
  if ( pdat == nullptr ) return nullptr;
  return pdat->getAdcData(nam.substr(ipos+1));
}
