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

#include "TGraph.h"
#include "TH2F.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "dune/DuneCommon/LineColors.h"

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
int test_WireSelector(string gname, double wireAngle, double minDrift, unsigned int nShow, int sigopt) {
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
  cout << myname << "      nShow: " << nShow << endl;
  cout << myname << "     sigopt: " << sigopt << endl;

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

  cout << myname << line << endl;
  cout << myname << "Geometry properties:" << endl;
  Index ncha = pgeo->Nchannels();
  cout << myname << "        # channels: " << ncha << endl;
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
      bool showWire = iwir < nShow;
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
  ws.fillDataMap();
  cout << myname << "  Expected wire count: " << nwirSel << endl;
  cout << myname << "  Reported wire count: " << ws.data().size() << endl;
  cout << myname << "    Mapped wire count: " << ws.dataMap().size() << endl;
  const WireSelector::WireSummary& wsum = ws.fillWireSummary();
  assert( ws.data().size() == nwirSel );
  assert( ws.dataMap().size() == nwirSel );
  assert( wsum.size() == nwirSel );

  // Build discriminated adcdata as a vector of x-values for each channel.
  cout << myname << line << endl;
  cout << myname << "Build ADC data." << endl;
  vector<vector<float>> adcdata(ncha);
  if ( sigopt == 2 ) {
    cout << myname << "  Strumming selected wires." << endl;
    Index nwir = ws.data().size();
    float dx = (wsum.xmax-wsum.xmin)/nwir;
    float xval = wsum.xmin + 0.5*dx;
    for ( Index iwir=0; iwir<nwir; ++iwir ) {
      const WireSelector::WireInfo& dat = ws.data()[iwir];
      adcdata[dat.channel].push_back(xval);
      xval += dx;
    }
  } else if ( sigopt == 1 ) {
    cout << myname << "  Strumming channels." << endl;
    float dx = (wsum.xmax-wsum.xmin)/ncha;
    float xval = wsum.xmin + 0.5*dx;
    for ( Index icha=0; icha<ncha; ++icha ) {
      adcdata[icha].push_back(xval);
      xval += dx;
    }
  }

  // Build a graph of z vs. x using adcdata.
  cout << myname << line << endl;
  cout << myname << "Build signal data." << endl;
  vector<float> xsigs;
  vector<float> zsigs;
  for ( Index icha=0; icha<ncha; ++icha ) {
    // Loop over the wires read out by this channel.
    auto range = ws.dataMap().equal_range(icha);
    for ( auto ient=range.first; ient!=range.second; ++ient ) {
      const WireSelector::WireInfo* pdat = ient->second;
      float zsig = pdat->z;
      // Loop over the ticks hit for this channel.
      for ( float xsig : adcdata[icha] ) {
         // Create a space point for wire and tick.
         xsigs.push_back(xsig);
         zsigs.push_back(zsig);
      }
    }
  }
  Index nsig = xsigs.size();
  TGraph gsigxz(nsig, &zsigs[0], &xsigs[0]);
  TGraph gsigzx(nsig, &xsigs[0], &zsigs[0]);
  cout << myname << "  # signals: " << nsig << endl;
  Index nShowSig = nShow < nsig ? nShow : nsig;
  for ( Index isig=0; isig<nShowSig; ++isig ) {
    ostringstream sout;
    sout.precision(3);
    sout << myname << setw(10) << std::fixed << xsigs[isig] << "," << setw(9) << std::fixed << zsigs[isig];
    cout << sout.str() << endl;
  }
    
  bool drawWires = true;
  LineColors lc;
  double ticklen = 0.015;
  Index wpadx = 1200;
  Index wpady = 1000;
  if ( drawWires ) {
    double b = 20;
    double xpmin = wsum.xmin;
    double xpmax = wsum.xmin;
    if ( gname == "protodune_geo" ) {
      xpmin = -400.0;
      xpmax = 400.0;
    }
    ostringstream ssttl;
    ssttl << gname << " ";
    if ( useView ) {
      ssttl << "view=" << view;
    } else {
      ssttl << "wireAngle=" << std::setprecision(3) << std::fixed << wireAngle;
    }
    if ( minDrift > 0.0 ) {
      ssttl << ", TPCs with drift > " << std::setprecision(0) << std::fixed << minDrift << " cm";
    }
    string sttlxz = ssttl.str() + "; z [cm]; x [cm]";
    string sttlzx = ssttl.str() + "; x [cm]; z [cm]";
    TH2* phaxz = new TH2F("phaxz", sttlxz.c_str(), 1, wsum.zmin-b, wsum.zmax+b, 1, xpmin, xpmax);
    TH2* phazx = new TH2F("phazx", sttlzx.c_str(), 1, xpmin, xpmax, 1, wsum.zmin-b, wsum.zmax+b);
    phaxz->SetStats(0);
    phazx->SetStats(0);
    TGraph gxz(wsum.size(), &wsum.zWire[0], &wsum.xWire[0]);
    TGraph gzx(wsum.size(), &wsum.xWire[0], &wsum.zWire[0]);
    gxz.SetMarkerColor(lc.red());
    gzx.SetMarkerColor(lc.red());
    TGraph gxzd(wsum.size(), &wsum.zWire[0], &wsum.xCathode[0]);
    TGraph gzxd(wsum.size(), &wsum.xCathode[0], &wsum.zWire[0]);
    gxzd.SetMarkerColor(lc.green());
    gzxd.SetMarkerColor(lc.green());
    TPadManipulator padxz(wpadx, wpady);
    padxz.add(phaxz);
    padxz.add(&gxz, "P");
    padxz.add(&gxzd, "P");
    if ( sigopt ) padxz.add(&gsigxz, "P");
    padxz.addAxis();
    padxz.setTickLength(ticklen);
    padxz.print("test_WireSelector_xz.png");
    TPadManipulator padzx(wpadx, wpady);
    padzx.add(phazx);
    padzx.add(&gzx, "P");
    padzx.add(&gzxd, "P");
    if ( sigopt ) padzx.add(&gsigzx, "P");
    padzx.addAxis();
    padzx.setTickLength(ticklen);
    padzx.print("test_WireSelector_zx.png");
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, const char* argv[]) {
  string gname = "protodune_geo";
  double thtx = 102;
  double minDrift = 0.0;
  unsigned int nShow = 0;
  int sigopt = 0;
  if ( argc > 1 ) {
    string sarg = argv[1];
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [gname] [thtx] [minDrift] [nShow] [sigopt]" << endl;
      cout << "  gname: Geometry name, e.g. protodune_geo" << endl;
      cout << "  thtx: Wire angle, e.g. 0 for vertical" << endl;
      cout << "        >= 100 means use view thtx-100, e.g. 2 for kZ" << endl;
      cout << "  minDrift: select TPC with drift > midDrift cm" << endl;
      cout << "  nShow: # wires, signals to print" << endl;
      cout << "  sigopt: Option to add signals to the displays" << endl;
      cout << "          0 for no signal" << endl;
      cout << "          1 for strum of detector channels" << endl;
      cout << "          2 for strum of selected wires" << endl;
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
          ssarg >> nShow;
          if ( argc > 5 ) {
            istringstream ssarg(argv[5]);
            ssarg >> sigopt;
          }
        }
      }
    }
  }
  test_WireSelector(gname, thtx, minDrift, nShow, sigopt);
  cout << "Tests concluded." << endl;
  ArtServiceHelper::close();
  return 0;
}

//**********************************************************************

#endif
