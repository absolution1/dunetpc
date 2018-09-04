// WireSelector.cxx

#include "WireSelector.h"
#include "larcore/Geometry/Geometry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include <string>
#include <iostream>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

WireSelector::WireSelector(Index icry) {
  selectCryostats({icry});
}

//**********************************************************************

const WireSelector::GeometryCore* WireSelector::geometry() {
  if ( m_pgeo == nullptr ) {
    art::ServiceHandle<geo::Geometry> hgeo;
    m_pgeo = hgeo.get();
  }
  return m_pgeo;
}

//**********************************************************************

void WireSelector::selectCryostats(const IndexVector& a_icrys) {
  m_icrys = a_icrys;
  selectPlanes();
}

//**********************************************************************

void WireSelector::selectView(View a_view) {
  m_view = a_view;
  selectPlanes();
}

//**********************************************************************

void WireSelector::selectWireAngle(double angle, double tol) {
  m_wireAngle = angle;
  m_wireAngleTolerance = tol;
  selectPlanes();
}

//**********************************************************************

void WireSelector::selectDrift(double dmin, double dmax) {
  m_driftMin = dmin;
  m_driftMax = dmax;
  selectPlanes();
}

//**********************************************************************

void WireSelector::selectTpcSet(Index itps) {
  m_tpcSets.clear();
  m_tpcSets.push_back(itps);
  selectPlanes();
}
//**********************************************************************

void WireSelector::selectTpcSets(const IndexVector& itpss) {
  m_tpcSets = itpss;
  selectPlanes();
}

//**********************************************************************

void WireSelector::selectPlanes() {
  const double pi = acos(-1.0);
  const double piOver2 = 0.5*pi;
  m_pids.clear();
  clearData();
  for ( Index icry : cryostats() ) {
    geo::CryostatID cid(icry);
    for ( geo::TPCID tid : geometry()->IterateTPCIDs(cid) ) {
      // Check TPC set.
      if ( m_tpcSets.size() ) {
        Index itps = geometry()->TPCtoTPCset(tid).TPCset;
        if ( find(m_tpcSets.begin(), m_tpcSets.end(), itps) == m_tpcSets.end() ) continue;
      }
      // Check drift.
      if ( driftMin() > 0.0 || driftMax() < 1.e20 ) {
        const geo::TPCGeo& gtpc = geometry()->TPC(tid);
        //double driftSize = gtpc.ActiveWidth();
        double driftSize = gtpc.DriftDistance(); // good for the last anode plane
        if ( driftSize < driftMin() ) continue;
        if ( driftSize >= driftMax() ) continue;
      }
      for ( geo::PlaneID pid : geometry()->IteratePlaneIDs(tid) ) {
        const geo::PlaneGeo& gpla = geometry()->Plane(pid);
        // Check view.
        if ( view() != geo::kUnknown && gpla.View() != view() ) continue;
        // Check wire angle.
        if ( wireAngleTolerance() < pi ) {
          double planeWireAngle = gpla.ThetaZ() - piOver2;
          double dang = planeWireAngle - wireAngle();
          while ( dang < -piOver2 ) dang += pi;
          while ( dang >  piOver2 ) dang -= pi;
          if ( fabs(dang) >= wireAngleTolerance() ) continue;
        }
        // Keep this plane.
        m_pids.push_back(pid);
      }
    }
  }
}

//**********************************************************************

const WireSelector::WireInfoVector& WireSelector::fillData() {
  const string myname = "WireSelector::fillData: ";
  const double piOver2 = 0.5*acos(-1.0);
  if ( haveData() ) return m_data;
  for ( geo::PlaneID pid : planeIDs() ) {
    const geo::PlaneGeo& gpla = geometry()->Plane(pid);
    const geo::TPCGeo& gtpc = geometry()->TPC(pid);
    const geo::PlaneGeo& gplaLast = gtpc.LastPlane();  // Plane furthest from TPC center.
    double xLastPlane = gplaLast.MiddleWire().GetCenter<TVector3>().x();
    double xThisPlane = gpla.MiddleWire().GetCenter<TVector3>().x();
    double driftOffset = fabs(xThisPlane - xLastPlane);
    for ( Index iwir=0; iwir<gpla.Nwires(); ++iwir ) {
      geo::WireID wid(pid, iwir);
      const geo::WireGeo* pgwir = gpla.WirePtr(wid);
      double wireAngle = gpla.ThetaZ() - piOver2;
      TVector3 xyzWire = pgwir->GetCenter<TVector3>();
      xyzWire.RotateX(wireAngle);
      double driftSign = gpla.GetNormalDirection().x();
      if ( fabs(fabs(driftSign) - 1.0) > 0.001 ) {
        cout << myname << "ERROR: Plane normal is not along x." << endl;
        continue;
      }
      double driftDist = gtpc.DriftDistance() - driftOffset;
      m_data.emplace_back(xyzWire.x(), xyzWire.y(), xyzWire.z(),
                          driftSign*driftDist,
                          pgwir->Length(),
                          gpla.WirePitch(),
                          geometry()->PlaneWireToChannel(wid));
    }
  }
  m_haveData = true;
  return m_data;
}

//**********************************************************************

const WireSelector::WireInfoMap& WireSelector::fillDataMap() {
  if ( haveData() && m_datamap.size() == m_data.size() ) return m_datamap;
  for ( const WireInfo& dat : fillData() ) {
    m_datamap.emplace(dat.channel, &dat);
  }
  return m_datamap;
}

//**********************************************************************

const WireSelector::WireSummary& WireSelector::fillWireSummary() {
  if ( haveData() && m_wireSummary.size() == m_data.size() ) return m_wireSummary;
  WireSummary& ws = m_wireSummary;
  Index nwir = data().size();
  ws.xmin = 1000.0;
  ws.ymin = 1000.0;
  ws.zmin = 1000.0;
  ws.xmax = -1000.0;
  ws.ymax = -1000.0;
  ws.zmax = -1000.0;
  ws.xWire.resize(nwir);
  ws.xCathode.resize(nwir);
  ws.zWire.resize(nwir);
  for ( Index iwir=0; iwir<nwir; ++iwir ) {
    const WireInfo& dat = data()[iwir];
    if ( dat.x1() < ws.xmin ) ws.xmin = dat.x1();
    if ( dat.y1() < ws.ymin ) ws.ymin = dat.y1();
    if ( dat.z1() < ws.zmin ) ws.zmin = dat.z1();
    if ( dat.x2() > ws.xmax ) ws.xmax = dat.x2();
    if ( dat.y2() > ws.ymax ) ws.ymax = dat.y2();
    if ( dat.z2() > ws.zmax ) ws.zmax = dat.z2();
    ws.xWire[iwir] = dat.x;
    ws.xCathode[iwir] = dat.x + dat.driftMax;
    ws.zWire[iwir] = dat.z;
  }
  return ws;
}

//**********************************************************************

void WireSelector::clearData() {
  m_data.clear();
  m_datamap.clear();
  m_wireSummary.clear();
  m_haveData = false;
}

//**********************************************************************
