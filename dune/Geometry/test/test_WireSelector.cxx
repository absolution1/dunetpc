// test_WireSelector.cxx

#ifndef test_WireSelector_CXX
#define test_WireSelector_CXX

// David Adams
// November 2016
//
// Test the DUNE geometry including channel mapping.
//
// This test demonstrates how to configure and use the LArSoft Geometry
// service outside the art framework. DUNE geometry and geometry helper
// service are used.
//
// Note the geometry service requires the experiment-specific geometry
// helper with the channel map also be loaded. 
//
// The functions that set the expected values must be defined externally.

#undef NDEBUG

#include "larcore/Geometry/Geometry.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "../WireSelector.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::istringstream;
using std::ostringstream;
using std::setw;
using std::vector;
using std::abs;
using geo::View_t;
using geo::SigType_t;
using geo::CryostatID;
using geo::TPCID;
using geo::PlaneID;
using geo::WireID;
typedef readout::TPCsetID APAID;
using readout::ROPID;
using geo::Geometry;
using geo::CryostatGeo;
using geo::TPCGeo;
using geo::WireGeo;
using Index = WireSelector::Index;

string toString(const TVector3& xyz, int w =9) {
  ostringstream sout;
  sout.precision(3);
  sout << "("
       << setw(w) << std::fixed << xyz.x()
       << ", "
       << setw(w) << std::fixed << xyz.y()
       << ", "
       << setw(w) << std::fixed << xyz.z()
       << ")";
  return sout.str();
}

string toString(const WireSelector::WireInfo& dat, int w =9) {
  ostringstream sout;
  sout.precision(3);
  sout << "("
       << setw(w) << std::fixed << dat.x
       << ", "
       << setw(w) << std::fixed << dat.y
       << ", "
       << setw(w) << std::fixed << dat.z
       << ") "
       << " ["
       << setw(5) << dat.channel
       << "] "
       << " ("
       << setw(w) << std::fixed << dat.driftMax
       << ", "
       << setw(w) << std::fixed << dat.length
       << ", "
       << setw(w) << std::fixed << dat.pitch
       << ")";
  return sout.str();
}
int test_WireSelector(string gname, double wireAngle, double minDrift, unsigned int  nShowWires) {
  const string myname = "test_WireSelector: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "   Geometry: " << gname << endl;

  WireSelector::View view = geo::kUnknown;
  if ( wireAngle == 100 ) view =  geo::kU;
  if ( wireAngle == 101 ) view =  geo::kV;
  if ( wireAngle == 102 ) view =  geo::kZ;
  bool useView = wireAngle >= 100;
  cout << myname << line << endl;
  cout << myname << "Create configuration." << endl;
  const char* ofname = "test_WireSelector.fcl";
  bool useExistingFcl = false;
  if ( ! useExistingFcl ) {
    ofstream fout(ofname);
    fout << "#include \"geometry_dune.fcl\"" << endl;
    fout << "services.Geometry:                   @local::" + gname << endl;
    fout << "services.ExptGeoHelperInterface:     @local::dune_geometry_helper" << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& myash = ArtServiceHelper::instance();
  myash.setLogLevel(3);

  cout << myname << line << endl;
  cout << myname << "Add services from " << ofname << endl;
  assert( myash.addServices(ofname, true) == 0 );

  cout << myname << line << endl;
  cout << myname << "Display services" << endl;
  myash.print();

  cout << myname << line << endl;
  cout << myname << "Load the services." << endl;
  assert( myash.loadServices() == 1 );
  myash.print();

  cout << myname << line << endl;
  cout << myname << "Get Geometry service." << endl;
  art::ServiceHandle<geo::Geometry> pgeo;

  cout << myname << line << endl;
  cout << myname << "Input arguments:" << endl;
  cout << myname << "      gname: " << gname << endl;
  cout << myname << "  wireAngle: " << wireAngle << endl;
  cout << myname << "   minDrift: " << minDrift << endl;
  cout << myname << " nShowWires: " << nShowWires << endl;

  cout << myname << line << endl;
  cout << myname << "Construct selector." << endl;
  WireSelector ws;
  if ( useView ) {
    if ( view != geo::kUnknown ) {
      cout << myname << "  Selecting view " << view << endl;
      ws.selectView(view);
    }
  } else {
    cout << myname << "  Selecting wire angle " << wireAngle << endl;
    ws.selectWireAngle(wireAngle);
  }
  if ( minDrift > 0.0 ) {
    cout << myname << "  Selecting min drift " << minDrift << " cm" << endl;
    ws.selectDrift(minDrift);
  }
  assert( ws.geometry() != nullptr );

  cout << myname << line << endl;
  cout << myname << "Selector properties:" << endl;
  cout << myname << "   Geometry name: " << pgeo->DetectorName() << endl;
  cout << myname << "            View: " << ws.view() << endl;
  cout << myname << "      Wire angle: " << ws.wireAngle() << endl;
  cout << myname << "  Wire angle tol: " << ws.wireAngleTolerance() << endl;
  cout << myname << "       Drift min: " << ws.driftMin() << endl;
  cout << myname << "        # planes: " << ws.planeIDs().size() << endl;
  Index ncha = pgeo->Nchannels();
  const double piOver2 = 0.5*acos(-1.0);
  Index nwirSel = 0;
  for ( WireSelector::PlaneID pid : ws.planeIDs() ) {
    const geo::PlaneGeo& gpla = ws.geometry()->Plane(pid);
    double thtx = gpla.ThetaZ() - piOver2;
    Index nwir = gpla.Nwires();
    const geo::TPCGeo& gtpc = pgeo->TPC(pid);
    double driftSize = gtpc.ActiveWidth();
    cout << myname << pid
         << "  view=" << gpla.View()
         << " thtx=" << std::fixed << std::setprecision(3) << thtx
         << " nwir=" << nwir
         << " drift distance =" << driftSize << " cm"
         << endl;
    TVector3 wdir = gpla.GetWireDirection();
    TVector3 ndir = gpla.GetNormalDirection();
    cout << myname << "    Wire direction: " << toString(wdir) << endl;
    cout << myname << "  Normal direction: " << toString(ndir) << endl;
    Index icha1 = ncha;
    Index icha2 = ncha;
    Index icha = ncha;
    string wireText;
    bool usePlane = useView ? gpla.View() == view : fabs(thtx - wireAngle)<0.001;
    if ( ! usePlane ) continue;
    nwirSel += nwir;
    for ( Index iwir=0; iwir<=nwir; ++iwir ) {
      bool showWire = iwir < nShowWires;
      bool endPlane = iwir == nwir;
      bool endBlock = endPlane;
      ostringstream wout;
      if ( ! endBlock ) {
        WireID wid(pid, iwir);
        const WireGeo* pgwir = pgeo->WirePtr(wid);
        TVector3 xyzWire = pgwir->GetCenter<TVector3>();
        TVector3 xyz1 = pgwir->GetStart<TVector3>();
        TVector3 xyz2 = pgwir->GetEnd<TVector3>();
        icha = pgeo->PlaneWireToChannel(wid);
        wout << myname << setw(12) << iwir << " [" << setw(4) << icha << "] " << toString(xyzWire);
        xyzWire.RotateX(thtx);
        wout << " --> " << toString(xyzWire) << "\n";
        wout << myname << setw(20) << "" << toString(xyz1);
        xyz1.RotateX(thtx);
        wout << " --> " << toString(xyz1) << "\n";
        wout << myname << setw(20) << "" << toString(xyz2);
        xyz2.RotateX(thtx);
        wout << " --> " << toString(xyz2) << "\n";
        assert( icha < ncha );
        if ( icha1 == ncha ) {
          icha1 = icha;
          icha2 = icha;
          wireText = showWire ? wout.str() : "";
        } else if ( icha == icha2 + 1 ) {
          icha2 = icha;
          if ( showWire ) wireText += wout.str();
        } else {
          endBlock = true;
        }
      }
      if ( endBlock ) {
        cout << myname << "      [" << icha1 << ", " << icha2 << "]" << endl;
        if ( endPlane ) {
          icha1 = ncha;
          icha2 = ncha;
        } else {
          icha1 = icha;
          icha2 = icha;
        }
        cout << wireText;
        wireText = showWire ? wout.str() : "";
      }
    }
    assert( icha1 == ncha );
    assert( icha2 == ncha );
  }


  cout << myname << line << endl;
  cout << myname << "Check geometry: " << gname << endl;
  assert( ws.geometry() == pgeo.get() );

  cout << myname << line << endl;
  cout << myname << "Check data" << endl;
  ws.fillData();
  cout << myname << "  Expected wire count: " << nwirSel << endl;
  cout << myname << "  Reported wire count: " << ws.data().size() << endl;
  Index iwir = 0;
  for ( const WireSelector::WireInfo& dat : ws.data() ) {
    if ( iwir++ < nShowWires ) cout << setw(6) << iwir << ": " << toString(dat) << endl;
  }
  assert( ws.data().size() == nwirSel );

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, const char* argv[]) {
  string gname = "protodune_geo";
  double thtx = 102;
  double minDrift = 0.0;
  unsigned int nShowWires = 0;
  if ( argc > 1 ) {
    string sarg = argv[1];
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [gname] [thtx] [minDrift] [nShowWires]" << endl;
      return 0;
    }
    gname = sarg;
    if ( argc > 2 ) {
      istringstream ssarg(argv[2]);
      ssarg >> thtx;
      if ( argc > 3 ) {
        istringstream ssarg(argv[3]);
        ssarg >> minDrift;
        if ( argc > 4 ) {
          istringstream ssarg(argv[4]);
          ssarg >> nShowWires;
        }
      }
    }
  }
  test_WireSelector(gname, thtx, minDrift, nShowWires);
  cout << "Tests concluded." << endl;
  ArtServiceHelper::close();
  return 0;
}

//**********************************************************************

#endif
