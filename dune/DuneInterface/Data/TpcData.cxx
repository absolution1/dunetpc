// TpcData.cxx

#include "dune/DuneInterface/Data/TpcData.h"

//**********************************************************************

TpcData::TpcData() : m_parent(nullptr) { }

//**********************************************************************

TpcData* TpcData::addTpcData(Name nam) {
  if ( nam == "" || nam == "." ) return nullptr;
  Name::size_type ipos = nam.rfind("/");
  if ( ipos != Name::npos ) {
    TpcData* pdat = getTpcData(nam.substr(0, ipos));
    if ( pdat == nullptr ) return nullptr;
    return pdat->addTpcData(nam.substr(ipos+1));
  }
  if ( m_dat.count(nam) ) return nullptr;
  m_dat[nam].m_parent = this;
  return &m_dat[nam];
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

const TpcData* TpcData::getTpcData(Name nam) const {
  if ( nam == "" || nam == "." ) return this;
  Name::size_type ipos = nam.find("/");
  if ( ipos == Name::npos ) {
    TpcDataMap::const_iterator idat = m_dat.find(nam);
    if ( idat == m_dat.end() ) return nullptr;
    return &idat->second;
  }
  const TpcData* pdat = getTpcData(nam.substr(0,ipos));
  if ( pdat == nullptr ) return nullptr;
  return pdat->getTpcData(nam.substr(ipos+1));
}

//**********************************************************************

TpcData::AdcData* TpcData::createAdcData() {
  if ( m_padcs ) return nullptr;
  m_padcs.reset(new AdcData);
  return getAdcData();
}

//**********************************************************************

void TpcData::clearAdcData() {
  m_padcs.reset(nullptr);
}

//**********************************************************************

TpcData::AdcData* TpcData::findAdcData() {
  AdcData* pacm = getAdcData();
  if ( pacm != nullptr ) return pacm;
  if ( getParent() == nullptr ) return nullptr;
  return getParent()->findAdcData();
}

//**********************************************************************

const TpcData::AdcData* TpcData::findAdcData() const {
  const AdcData* pacm = getAdcData();
  if ( pacm != nullptr ) return pacm;
  if ( getParent() == nullptr ) return nullptr;
  return getParent()->findAdcData();
}

