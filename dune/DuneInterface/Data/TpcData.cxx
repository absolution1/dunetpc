// TpcData.cxx

#include "dune/DuneInterface/Data/TpcData.h"

using std::cout;
using std::endl;

//**********************************************************************

TpcData::TpcData() : m_parent(nullptr) { }

//**********************************************************************

TpcData::TpcData(Index npla) : m_parent(nullptr), m_adcs(npla) {
  for ( Index ipla=0; ipla<npla; ++ipla ) m_adcs[ipla].reset(new AdcChannelDataMap);
}

//**********************************************************************

TpcData::TpcData(const AdcDataVector& adcs) : m_parent(nullptr), m_adcs(adcs) { }

//**********************************************************************

TpcData* TpcData::addTpcData(Name nam, bool copyAdcData) {
  if ( nam == "" || nam == "." ) return nullptr;
  Name::size_type ipos = nam.rfind("/");
  if ( ipos != Name::npos ) {
    TpcData* pdat = getTpcData(nam.substr(0, ipos));
    if ( pdat == nullptr ) return nullptr;
    return pdat->addTpcData(nam.substr(ipos+1), copyAdcData);
  }
  
  if ( m_dat.count(nam) ) return nullptr;
  TpcData& tpc = m_dat[nam];
  if ( copyAdcData ) tpc.m_adcs = m_adcs;
  tpc.m_parent = this;
  return &tpc;
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

TpcData::Index TpcData::getTpcData(Name nam, TpcDataVector& out) {
  Index nerr = 0;
  if ( nam == "" || nam == "." ) {
    out.push_back(this);
    return 0;
  }
  if ( nam == "*" ) {
    for ( TpcDataMap::value_type& idat : m_dat ) {
      out.push_back(&idat.second);
    }
    return 0;
  }
  // Split nam --> dnam/rnam
  Name::size_type ipos = nam.find("/");
  if ( ipos == Name::npos ) {
    TpcData* pdat = getTpcData(nam);
    if ( pdat == nullptr ) ++nerr;
    else out.push_back(pdat);
  } else {
    Name dnam = nam.substr(0,ipos);
    Name rnam = ipos == Name::npos ? "" : nam.substr(ipos+1);
    TpcDataVector dirs;
    nerr += getTpcData(dnam, dirs);
    for ( TpcData* pdat : dirs ) nerr += pdat->getTpcData(rnam, out);
  }
  return nerr;
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

TpcData::AdcDataPtr TpcData::createAdcData(bool updateParent) {
  return addAdcData(AdcDataPtr(new AdcChannelDataMap), updateParent);
}

//**********************************************************************

TpcData::AdcDataPtr TpcData::addAdcData(AdcDataPtr padc, bool updateParent) {
  m_adcs.push_back(padc);
  if ( updateParent && m_parent != nullptr ) m_parent->addAdcData(padc, true);
  return padc;
}

//**********************************************************************

void TpcData::clearAdcData() {
  m_adcs.clear();
}

//**********************************************************************

std::ostream& TpcData::print(Name prefix, Index depth) const {
  Index nmap = getAdcData().size();
  cout << prefix << "ADC Channel map count is " << nmap;
  if ( nmap ) cout << ":";
  cout << endl;
  for ( const AdcDataPtr& pacm : getAdcData() ) {
    if ( ! pacm ) cout << "NULL" << endl;
    else {
      cout << prefix << "  " << "Channel count " << pacm->size();
      if ( pacm->size() ) {
        cout << " in range [" << pacm->begin()->first << ", " << pacm->rbegin()->first << "]";
      }
      cout << endl;
    }
  }
  cout << prefix << "2D ROI count is " << get2dRois().size() << endl;
  cout << prefix << "TPC data count is " << getData().size();
  if ( getData().size() ) cout << ":";
  cout << endl;
  if ( depth > 0 ) {
    for ( const TpcDataMap::value_type& itpd : getData() ) {
      cout << prefix << itpd.first << ":" << endl;
      itpd.second.print(prefix + "  ", depth - 1);
    }
  }
  return cout;
}

//**********************************************************************
