// IcebergChannelGeo.cxx

#include "IcebergChannelGeo.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include <string>
#include <iostream>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

IcebergChannelGeo::IcebergChannelGeo(Index icha)
: IcebergChannelGeo(icha, &*art::ServiceHandle<geo::Geometry>()) { }

//**********************************************************************

IcebergChannelGeo::IcebergChannelGeo(Index icha, const geo::GeometryCore* pgeo)
: m_icha(icha), m_valid(false) {
  const string myname = "IcebergChannelGeo::ctor: ";
  if ( pgeo == nullptr ) return;
  if ( icha >= pgeo->Nchannels() ) return;
  m_valid = true;
  if ( pgeo->DetectorName() != "iceberg" ) {
    cout << myname << "WARNING: Iceberg mod is applied to non-Iceberg detector channels." << endl;
  }
  double ylo = 173.819 - 95.23;
  for ( geo::WireID wid : pgeo->ChannelToWire(icha) ) {
    EndPoints ends = pgeo->WireEndPoints<TVector3>(wid);
    if ( ends.first.y() > ends.second.y() ) {
      Point firstPoint = ends.first;
      ends.first = ends.second;
      ends.second = firstPoint;
    }
    Point& plo = ends.first;
    Point& phi = ends.second;
    if ( phi.y() < ylo ) continue;
    if ( plo.y() < ylo ) {
      // Set lower limit to ylo.
      double slop = (plo.z() - phi.z())/(plo.y() - phi.y());
      double zlo = phi.z() + slop*(ylo - phi.y());
      plo.SetY(ylo);
      plo.SetZ(zlo);
    }
    if ( m_ends.size() == 0 ) {
       m_bot = ends.first;
       m_top = ends.second;
    } else {
       if ( ends.first.y() < m_bot.y() ) m_bot = ends.first;
       if ( ends.second.y() > m_top.y() ) m_top = ends.second;
    }
    m_ends.push_back(ends);
  }
}

//**********************************************************************
