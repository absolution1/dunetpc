// test_Geometry_Dune35t.cxx

// David Adams
// OCtober 2016
//
// Test the DUNE 35t geometry including channel mapping.
//
// This test demonstrates how to configure and use the LArSoft Geometry
// service outside the art framework. DUNE geometry and geometry helper
// service are used.
//
// Note the geometry service requires the experiment-specific geometry
// helper with the channel map also be loaded. 

#undef NDEBUG

#include "larcore/Geometry/Geometry.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::istringstream;
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

typedef unsigned int Index;
const Index InvalidIndex = std::numeric_limits<Index>::max();

template<class T>
using Vector = vector<T>;

template<class T>
using TwoVector = Vector<Vector<T>>;

template<class T>
using ThreeVector = Vector<TwoVector<T>>;

template<class T, class V>
void resize(Vector<T>& vec1, Index n1, const V& val) {
  vec1.resize(n1, val);
}

template<class T, class V>
void resize(TwoVector<T>& vec2, Index n1, Index n2, const V& val) {
  vec2.resize(n1);
  for ( Vector<T>& vec1 : vec2 ) {
    vec1.resize(n2, val);
  }
}

template<class T, class V>
void resize(ThreeVector<T>& vec3, Index n1, Index n2, Index n3, const V& val) {
  vec3.resize(n1);
  for ( TwoVector<T>& vec2 : vec3 ) {
    vec2.resize(n2);
    for ( Vector<T>& vec1 : vec2 ) {
      vec1.resize(n3, val);
    }
  }
}

struct SpacePoint {
  double x;
  double y;
  double z;
  SpacePoint(double ax, double ay, double az) : x(ax), y(ay), z(az) { }
};

//**********************************************************************

// Helper functions

namespace {

template<class T>
void check(string name, T val) {
  cout << name << ": " << val << endl;
}

template<class T, class V>
void check(string name, T val, V checkval) {
  check(name, val);
  T eval = checkval;
  if ( val != eval ) {
    cout << val << " != " << checkval << endl;
    assert(false);
  }
}

template<class T, class V>
void checkval(string name, T val, V chkval) {
  T eval = chkval;
  if ( val != eval ) {
    cout << name << ": " << val << " != " << chkval << endl;
    assert(false);
  }
}

void checkfloat(double x1, double x2, double tol =1.e-5) {
  double den = abs(x1) + abs(x2);
  double num = abs(x2 - x1);
  if ( num/den > tol ) {
    cout << "checkfloat: " << x1 << " != " << x2 << endl;
    assert(false);
  }
}

//**********************************************************************

// Helper classes.

struct ExpectedValues {
  typedef void (*FunPtr)(Geometry*);
  string fullname;
  Index ncry = 0;
  Index ntpc = 0;
  Index npla = 0;
  Index napa = 0;
  Index nrop = 0;
  Index nchaPerApa = 0;
  Index nchatot = 0;
  ThreeVector<View_t> view;
  ThreeVector<SigType_t> sigType;
  TwoVector<Index> nwirPerPlane;
  // For ROP tests.
  Vector<Index> nchaPerRop;
  Vector<Index> chacry;          // Expected cryostat for each channel;
  Vector<Index> chaapa;          // Expected APA for each channel;
  Vector<Index> charop;          // Expected ROP for each channel;
  ThreeVector<Index> firstchan;  // First channel(cry, apa, rop)
  // For space point tests.
  FunPtr pfun = nullptr;
  vector<double> xfs;            // x-fractions for space points
  vector<double> yfs;            // y-fractions for space points
  vector<double> zfs;            // z-fractions for space points
  Vector<SpacePoint> posXyz;     // Space points.
  vector<double> posTpc;         // TPCs for space points
  vector<double> posPla;         // TPC planes for space points
  vector<double> posWco;         // Wire coordinates for space points
};

//**********************************************************************

// 35t geometry.

void set35tSpacePoints(Geometry*);

void set35t(ExpectedValues& ev) {
  ev.fullname = "dune35t4apa_v6";
  // Geometry counts.
  ev.ncry = 1;
  ev.ntpc = 8;
  ev.npla = 3;
  ev.napa = 4;
  ev.nrop = 4;
  // Signal type and view for each TPC plane.
  resize(ev.view, ev.ncry, ev.ntpc, ev.npla, geo::kUnknown);
  resize(ev.sigType, ev.ncry, ev.ntpc, ev.npla, geo::kMysteryType);
  for ( Index icry=0; icry<ev.ncry; ++icry ) {
    for ( Index itpc=0; itpc<ev.ntpc; ++itpc ) {
      for ( Index ipla=0; ipla<ev.npla; ++ipla ) {
        ev.sigType[icry][itpc][ipla] = (ipla < 2) ? geo::kInduction : geo::kCollection;
      }
      ev.view[icry][itpc][0] = geo::kU;
      ev.view[icry][itpc][1] = geo::kV;
      ev.view[icry][itpc][2] = geo::kZ;
    }
  }
  // Wire count for each TPC plane.
  resize(ev.nwirPerPlane, ev.ntpc, ev.npla, 0);
  ev.nwirPerPlane[0][0] = 359;
  ev.nwirPerPlane[0][1] = 345;
  ev.nwirPerPlane[0][2] = 112;
  ev.nwirPerPlane[2][0] = 194;
  ev.nwirPerPlane[2][1] = 188;
  ev.nwirPerPlane[2][2] = 112;
  ev.nwirPerPlane[4][0] = 236;
  ev.nwirPerPlane[4][1] = 228;
  ev.nwirPerPlane[4][2] = 112;
  for ( Index ipla=0; ipla<ev.npla; ++ipla ) ev.nwirPerPlane[1][ipla] = ev.nwirPerPlane[0][ipla];
  for ( Index ipla=0; ipla<ev.npla; ++ipla ) ev.nwirPerPlane[3][ipla] = ev.nwirPerPlane[2][ipla];
  for ( Index ipla=0; ipla<ev.npla; ++ipla ) ev.nwirPerPlane[5][ipla] = ev.nwirPerPlane[4][ipla];
  for ( Index ipla=0; ipla<ev.npla; ++ipla ) ev.nwirPerPlane[6][ipla] = ev.nwirPerPlane[0][ipla];
  for ( Index ipla=0; ipla<ev.npla; ++ipla ) ev.nwirPerPlane[7][ipla] = ev.nwirPerPlane[0][ipla];
  resize(ev.nchaPerRop, ev.nrop, 0);
  ev.nchaPerRop[0] = 144;
  ev.nchaPerRop[1] = 144;
  ev.nchaPerRop[2] = 112;
  ev.nchaPerRop[3] = 112;
  for ( Index irop=0; irop<ev.nrop; ++irop ) ev.nchaPerApa += ev.nchaPerRop[irop];
  ev.nchatot = ev.napa*ev.nchaPerApa;
  resize(ev.chacry, ev.nchatot, InvalidIndex);
  resize(ev.chaapa, ev.nchatot, InvalidIndex);
  resize(ev.charop, ev.nchatot, InvalidIndex);
  Index icha = 0;
  for ( Index icry=0; icry<ev.ncry; ++icry ) {
    for ( Index iapa=0; iapa<ev.napa; ++iapa ) {
      for ( Index irop=0; irop<ev.nrop; ++irop ) {
        for ( Index kcha=0; kcha<ev.nchaPerRop[irop]; ++kcha ) {
          ev.chacry[icha] = icry;
          ev.chaapa[icha] = iapa;
          ev.charop[icha] = irop;
          ++icha;
        }
      }
    }
  }
  resize(ev.firstchan, ev.ncry, ev.napa, ev.nrop, raw::InvalidChannelID);
  Index chan = 0;
  for ( Index icry=0; icry<ev.ncry; ++icry ) {
    for ( Index iapa=0; iapa<ev.napa; ++iapa ) {
      for ( Index irop=0; irop<ev.nrop; ++irop ) {
        ev.firstchan[icry][iapa][irop] = chan;
        chan += ev.nchaPerRop[irop];
      }
    }
  }
  ev.pfun = set35tSpacePoints;
  #include "set35tSpacePoints.dat"
}

void set35tSpacePoints(Geometry* pgeo) {
  const string myname = "set35tSpacePoints: ";
  string ofname = "set35tSpacePoints.dat";
  cout << myname << "Writing " << ofname << endl;
  const CryostatGeo& crygeo = pgeo->Cryostat(0);
  double origin[3] = {0.0};
  double crypos[3] = {0.0};
  crygeo.LocalToWorld(origin, crypos);
  double cxlo = crypos[0] - crygeo.HalfWidth();
  double cxhi = crypos[0] + crygeo.HalfWidth();
  double cylo = crypos[1] - crygeo.HalfHeight();
  double cyhi = crypos[1] + crygeo.HalfHeight();
  double czlo = crypos[2] - 0.5*crygeo.Length();
  double czhi = crypos[2] + 0.5*crygeo.Length();
  cout << "Cryostat limits: "
       << "(" << cxlo << ", " << cylo << ", " << czlo << "), "
       << "(" << cxhi << ", " << cyhi << ", " << czhi << ")" << endl;
  vector<double> zfs = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };
  vector<double> yfs = {0.2, 0.5, 0.8};
  vector<double> xfs = {0.1, 0.3, 0.6 };
  ofstream fout(ofname.c_str());
  //fout << "void set35tSpacepoints(Geometry* pgeo) {" << endl;
  for ( double zf : zfs ) {
    double z = czlo + zf*(czhi-czlo);
    for ( double yf : yfs ) {
      double y = cylo + yf*(cyhi-cylo);
      for ( double xf : xfs ) {
        double x = cxlo + xf*(cxhi-cxlo);
        double xyz[3] = {x, y, z};
        TPCID tpcid = pgeo->FindTPCAtPosition(xyz);
        unsigned int itpc = tpcid.TPC;
        const TPCGeo& tpcgeo = pgeo->TPC(tpcid);
        unsigned int npla = tpcgeo.Nplanes();
        fout << "  ev.posXyz.push_back(SpacePoint(" << x << ", " << y << ", " << z << "));" << endl;
        for ( unsigned int ipla=0; ipla<npla; ++ipla ) {
          PlaneID plaid(tpcid, ipla);
          double xwire = pgeo->WireCoordinate(x, y, plaid);
          fout << "  ev.posTpc.push_back(" << itpc << ");" << endl;
          fout << "  ev.posPla.push_back(" << ipla << ");" << endl;
          fout << "  ev.posWco.push_back(" << xwire << ");" << endl;
        }  // end loop over planes
      }  // end loop over x
    }  // end loop over y
  }  // end loop over z
  //fout << "}" << endl;
}

//**********************************************************************

// 10kt workspace geometry.

void setWorkspaceSpacePoints(Geometry* pgeo);

void setWorkspace(ExpectedValues& ev) {
  ev.fullname = "dune10kt_v1_workspace";
  // Geometry counts.
  ev.ncry = 1;
  ev.ntpc = 8;
  ev.npla = 3;
  ev.napa = 4;
  ev.nrop = 4;
  // Signal type and view for each TPC plane.
  resize(ev.sigType, ev.ncry, ev.ntpc, ev.npla, geo::kMysteryType);
  resize(ev.view, 1, 8, 3, geo::kUnknown);
  for ( Index icry=0; icry<ev.ncry; ++icry ) {
    for ( Index itpc=0; itpc<ev.ntpc; ++itpc ) {
      for ( Index ipla=0; ipla<ev.npla; ++ipla ) {
        ev.sigType[icry][itpc][ipla] = (ipla < 2) ? geo::kInduction : geo::kCollection;
      }
      ev.view[icry][itpc][0] = geo::kU;
      ev.view[icry][itpc][1] = geo::kV;
      ev.view[icry][itpc][2] = geo::kZ;
    }
  }
  // Wire count for each TPC plane.
  resize(ev.nwirPerPlane, ev.ntpc, ev.npla, 0);
  ev.nwirPerPlane[0][0] = 1149;
  ev.nwirPerPlane[0][1] = 1148;
  ev.nwirPerPlane[0][2] =  480;
  for ( Index itpc=1; itpc<ev.ntpc; ++itpc ) {
    for ( Index ipla=0; ipla<ev.npla; ++ipla ) ev.nwirPerPlane[itpc][ipla] = ev.nwirPerPlane[0][ipla];
  }
  // Channel count per ROP.
  resize(ev.nchaPerRop, ev.nrop, 0);
  ev.nchaPerRop[0] = 1149;
  ev.nchaPerRop[1] = 144;
  ev.nchaPerRop[2] = 112;
  ev.nchaPerRop[3] = 112;
  ev.nchatot = 10240;
  // Space points.
  ev.pfun = setWorkspaceSpacePoints;
  #include "setWorkspaceSpacePoints.dat"
}

void setWorkspaceSpacePoints(Geometry* pgeo) {
  const string myname = "setWorkspacetSpacePoints: ";
  string ofname = "setWorkspaceSpacePoints.dat";
  cout << myname << "Writing " << ofname << endl;
  const CryostatGeo& crygeo = pgeo->Cryostat(0);
  double origin[3] = {0.0};
  double crypos[3] = {0.0};
  crygeo.LocalToWorld(origin, crypos);
  double cxlo = crypos[0] - crygeo.HalfWidth();
  double cxhi = crypos[0] + crygeo.HalfWidth();
  double cylo = crypos[1] - crygeo.HalfHeight();
  double cyhi = crypos[1] + crygeo.HalfHeight();
  double czlo = crypos[2] - 0.5*crygeo.Length();
  double czhi = crypos[2] + 0.5*crygeo.Length();
  cout << "Cryostat limits: "
       << "(" << cxlo << ", " << cylo << ", " << czlo << "), "
       << "(" << cxhi << ", " << cyhi << ", " << czhi << ")" << endl;
  vector<double> zfs = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };
  vector<double> yfs = {0.2, 0.5, 0.8};
  vector<double> xfs = {0.2, 0.3, 0.6 };
  ofstream fout(ofname.c_str());
  for ( double zf : zfs ) {
    double z = czlo + zf*(czhi-czlo);
    for ( double yf : yfs ) {
      double y = cylo + yf*(cyhi-cylo);
      for ( double xf : xfs ) {
        double x = cxlo + xf*(cxhi-cxlo);
        double xyz[3] = {x, y, z};
        cout << "  (" << x << ", " << y << ", " << z << ")" << endl;
        TPCID tpcid = pgeo->FindTPCAtPosition(xyz);
        unsigned int itpc = tpcid.TPC;
        if ( itpc == TPCID::InvalidID ) continue;
        const TPCGeo& tpcgeo = pgeo->TPC(tpcid);
        unsigned int npla = tpcgeo.Nplanes();
        fout << "  ev.posXyz.push_back(SpacePoint(" << x << ", " << y << ", " << z << "));" << endl;
        for ( unsigned int ipla=0; ipla<npla; ++ipla ) {
          PlaneID plaid(tpcid, ipla);
          double xwire = pgeo->WireCoordinate(x, y, plaid);
          fout << "  ev.posTpc.push_back(" << itpc << ");" << endl;
          fout << "  ev.posPla.push_back(" << ipla << ");" << endl;
          fout << "  ev.posWco.push_back(" << xwire << ");" << endl;
        }  // end loop over planes
      }  // end loop over x
    }  // end loop over y
  }  // end loop over z
}

//**********************************************************************

}  //end unnamed namespace

//**********************************************************************

// Declare test functions.

int test_Geometry_Dune35t(string gname ="dune35t_geo", string chanmap ="Dune35tChannelMapAlg", bool dorop =true,
                          Index maxchanprint =10);

int test_GeometryWithExpectedValues(const ExpectedValues& ev,
                                    string gname ="dune35t_geo",
                                    string chanmap ="Dune35tChannelMapAlg",
                                    bool dorop =true,
                                    Index maxchanprint =10);

//**********************************************************************

int test_Geometry_Dune35t(string gname, string chanmap, bool dorop, Index maxchanprint) {
  const string myname = "test_Geometry_Dune35t: ";
  ExpectedValues ev;
  if ( gname == "dune35t_geo" ) set35t(ev);
  else if ( gname == "dune10kt_workspace_geo" ) setWorkspace(ev);
  else {
    cout << myname << "Unknown geometry specifier: " << gname << endl;
    assert(false);
    return 1;
  }
  return test_GeometryWithExpectedValues(ev, gname, chanmap, dorop, maxchanprint);
}

//**********************************************************************

int test_GeometryWithExpectedValues(const ExpectedValues& ev, string gname, string chanmap,
                                    bool dorop, Index maxchanprint) {
  const string myname = "test_Geometry_Dune35t: ";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "   Geometry: " << gname << endl;
  cout << myname << "Channel map: " << chanmap << endl;
  cout << myname << "     Do ROP: " << dorop << endl;

  cout << myname << line << endl;
  cout << myname << "Create configuration." << endl;
  cout << myname << gname << endl;
  cout << myname << chanmap << endl;
  const char* ofname = "test_Geometry_Dune35t.fcl";
  {
    ofstream fout(ofname);
    fout << "#include \"geometry_dune.fcl\"" << endl;
    fout << "services.Geometry:                   @local::" + gname << endl;
    fout << "services.ExptGeoHelperInterface:     @local::dune_geometry_helper" << endl;
    if ( chanmap.size() ) {
      fout << "services.ExptGeoHelperInterface.ChannelMapClass: " << chanmap << endl;
    }
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
  check("Default wiggle", pgeo->DefaultWiggle());
  check("Geometry name", pgeo->DetectorName(), ev.fullname);
  cout << myname << "ROOT name: " << pgeo->ROOTFile() << endl;
  cout << myname << "GDML name: " << pgeo->GDMLFile() << endl;

  cout << myname << line << endl;
  double xlo, ylo, zlo;
  double xhi, yhi, zhi;
  cout << myname << "World box"  << endl;
  pgeo->WorldBox(&xlo, &ylo, &zlo, &xhi, &yhi, &zhi);
  check("  xlo", xlo);
  check("  ylo", ylo);
  check("  zlo", zlo);
  check("  xhi", xhi);
  check("  yhi", yhi);
  check("  zhi", zhi);

  cout << myname << line << endl;
  check("SurfaceY", pgeo->SurfaceY());

  cout << myname << line << endl;
  if ( 0  ) check("GetWorldVolumeName", pgeo->GetWorldVolumeName());
  check("TotalMass", pgeo->TotalMass());

  cout << myname << line << endl;
  check("CryostatHalfWidth", pgeo->CryostatHalfWidth());
  check("CryostatHalfHeight", pgeo->CryostatHalfHeight());
  check("CryostatLength", pgeo->CryostatLength());
  //check("", pgeo->());

  cout << myname << line << endl;
  Index ncry = pgeo->Ncryostats();
  check("Ncryostats", ncry, ev.ncry);
  check("MaxTPCs", pgeo->MaxTPCs(), ev.ntpc);
  check("MaxPlanes", pgeo->MaxPlanes(), ev.npla);
  check("TotalNTPC", pgeo->TotalNTPC(), ev.ntpc);
  check("Nviews", pgeo->Nviews(), ev.npla);
  check("Nchannels", pgeo->Nchannels(), ev.nchatot);

  cout << myname << line << endl;
  cout << "Check TPC wire plane counts." << endl;
  for ( Index icry=0; icry<ncry; ++icry ) {
    Index ntpc = pgeo->NTPC(icry);
    cout << "  Cryostat " << icry << " has " << ntpc << " TPCs" << endl;
    for ( Index itpc=0; itpc<ntpc; ++itpc ) {
      Index npla = pgeo->Nplanes(itpc, icry);
      cout << "    TPC " << itpc << " has " << npla << " planes" << endl;
      assert( npla == 3 );
      for ( Index ipla=0; ipla<npla; ++ipla ) {
        Index nwir = pgeo->Nwires(ipla, itpc, icry);
        cout << "      Plane " << ipla << " has " << nwir << " wires" << endl;
        assert( nwir == ev.nwirPerPlane[itpc][ipla] );
      }
    }
  }

  cout << myname << line << endl;
  cout << "Check wire planes." << endl;
  for ( PlaneID plaid : pgeo->IteratePlaneIDs() ) {
    cout << "  Plane " << plaid << endl;
    cout << "    Signal type: " << pgeo->SignalType(plaid) << endl;
    cout << "           View: " << pgeo->View(plaid) << endl;
    assert( pgeo->SignalType(plaid) == ev.sigType[plaid.Cryostat][plaid.TPC][plaid.Plane] );
    assert( pgeo->View(plaid) == ev.view[plaid.Cryostat][plaid.TPC][plaid.Plane] );
  }

  cout << myname << line << endl;
  cout << "Check channel-wire mapping." << endl;
  Index itpc1Last = TPCID::InvalidID;
  Index ipla1Last = WireID::InvalidID;
  Index nprint = 0;
  TwoVector<Index> lastwire;
  resize(lastwire, ev.ntpc, ev.npla, 0);
  for ( Index itpc=0; itpc<ev.ntpc; ++itpc ) 
    for ( Index ipla=0; ipla<ev.npla; ++ipla )
      lastwire[itpc][ipla] = 0;
  for ( Index icha=0; icha<ev.nchatot; ++icha ) {
    vector<WireID> wirids = pgeo->ChannelToWire(icha);
    assert( wirids.size() > 0 );
    WireID wirid1 = wirids[0];
    Index itpc1 = wirid1.TPC;
    Index iapa1 = itpc1/2;
    Index ipla1 = wirid1.Plane;
    if ( itpc1 != itpc1Last || ipla1 != ipla1Last ) nprint = 0;
    itpc1Last = itpc1;
    ipla1Last = ipla1;
    bool print = nprint < maxchanprint;
    if ( print ) ++nprint;
    if ( print ) cout << "  Channel " << setw(4) << icha << " has " << wirids.size() << " wires:";
    for ( WireID wirid : wirids ) {
      Index itpc = wirid.TPC;
      Index iapa = itpc/2;
      Index ipla = wirid.Plane;
      Index iwir = wirid.Wire;
      if ( print ) cout << " " << itpc << "-" << ipla << "-"<< iwir;
      if ( iwir > lastwire[itpc][ipla] ) lastwire[itpc][ipla] = iwir;
      assert( iapa == iapa1 );
      assert( ipla == ipla1 );
      checkval("\nPlaneWireToChannel", pgeo->PlaneWireToChannel(wirid), icha );
      assert( pgeo->PlaneWireToChannel(wirid) == icha );
      assert( pgeo->SignalType(icha) == ev.sigType[wirid.Cryostat][wirid.TPC][wirid.Plane] );
      assert( pgeo->View(icha) == ev.view[wirid.Cryostat][wirid.TPC][wirid.Plane] );
    }
    if ( print ) cout << endl;
  }
  for ( Index itpc=0; itpc<ev.ntpc; ++itpc ) {
    for ( Index ipla=0; ipla<ev.npla; ++ipla ) {
      Index nwir = lastwire[itpc][ipla] + 1;
      cout << "  TPC-plane " << itpc << "-" << ipla << " has " << setw(3) << nwir << " wires" << endl;
      assert( nwir == ev.nwirPerPlane[itpc][ipla] );
    }
  }

  if ( dorop ) {
    cout << myname << line << endl;
    cout << "Check ROP counts and channels." << endl;
    check("MaxROPs", pgeo->MaxROPs(), ev.nrop);
    Index icry = 0;
    for ( CryostatID cryid: pgeo->IterateCryostatIDs() ) {
      Index napa = pgeo->NTPCsets(cryid);
      cout << "  Cryostat " << icry << " has " << napa << " APAs" << endl;
      assert( napa == ev.napa );
      for ( Index iapa=0; iapa<napa; ++iapa ) {
        APAID apaid(cryid, iapa);
        Index nrop = pgeo->NROPs(apaid);
        cout << "    APA " << iapa << " has " << nrop << " ROPs" << endl;
        assert( nrop == ev.nrop );
        for ( Index irop=0; irop<nrop; ++irop ) {
          ROPID ropid(apaid, irop);
          Index ncha = pgeo->Nchannels(ropid);
          Index icha1 = pgeo->FirstChannelInROP(ropid);
          Index icha2 = icha1 + ncha - 1;
          cout << "      ROP " << irop << " has " << ncha << " channels: ["
               << icha1 << ", " << icha2 << "]" << endl;
          assert( ncha == ev.nchaPerRop[irop] );
          assert( icha1 == ev.firstchan[icry][iapa][irop] );
        }
      }
      ++icry;
    }
    assert( icry == ncry );
    cout << myname << line << endl;
    cout << "Check channel-ROP mapping." << endl;
    icry = 0;
    for ( Index icha=0; icha<ev.nchatot; ++icha ) {
      ROPID ropid = pgeo->ChannelToROP(icha);
      Index icry = ropid.Cryostat;
      Index iapa = ropid.TPCset;
      Index irop = ropid.ROP;
      assert( icry == ev.chacry[icha] );
      assert( iapa == ev.chaapa[icha] );
      assert( irop == ev.charop[icha] );
    }
      
  } else {
    cout << myname << line << endl;
    cout << "Skipped APA and ROP tests." << endl;
  }  // end dorop

  cout << myname << line << endl;
  Index nspt = ev.posXyz.size();
  if ( nspt == 0 ) {
    cout << myname << "No space points found." << endl;
    if ( ev.pfun ) {
      cout << myname << "Creating space points." << endl;
      (*ev.pfun)(&*pgeo);
    }
  } else {
    cout << myname << "Check " << nspt << " space points." << endl;
    assert( ev.posPla.size() == ev.posTpc.size() );
    assert( ev.posWco.size() == ev.posTpc.size() );
    Index ires = 0;
    int w = 9;
    for ( Index ispt=0; ispt<nspt; ++ispt ) {
      double x = ev.posXyz[ispt].x;
      double y = ev.posXyz[ispt].y;
      double z = ev.posXyz[ispt].z;
      double xyz[3] = {x, y, z};
      TPCID tpcid = pgeo->FindTPCAtPosition(xyz);
      cout << "  (" << setw(w) << x << "," << setw(w) << y << "," << setw(w) << z << "): " << tpcid << endl;
      assert( tpcid.Cryostat != CryostatID::InvalidID );
      unsigned int itpc = tpcid.TPC;
      assert( itpc != TPCID::InvalidID );
      const TPCGeo& tpcgeo = pgeo->TPC(tpcid);
      unsigned int npla = tpcgeo.Nplanes();
      assert( npla == ev.npla );
      assert( npla > 0 );
      for ( unsigned int ipla=0; ipla<npla; ++ipla ) {
        assert( ires < ev.posTpc.size() );
        PlaneID plaid(tpcid, ipla);
        WireID wirid = pgeo->NearestWireID(xyz, plaid);
        double xwire = pgeo->WireCoordinate(x, y, plaid);
        cout << "    TPC " << setw(2) << itpc << " plane " << ipla
             << " nearest wire is " << setw(3) << wirid.Wire
             << " and coordinate is " << xwire << endl;
        assert( itpc == ev.posTpc[ires] );
        assert( ipla == ev.posPla[ires] );
        checkfloat(xwire, ev.posWco[ires], 2.e-4);
        ++ires;
      }  // end loop over planes
    }  // end loop over space points
    cout << "  # checked planes: " << ires << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, const char* argv[]) {
  string gname = "dune35t_geo";
  string chanmap = "";
  bool dorop = false;
  Index maxchanprint = 10;
  if ( argc > 1 ) {
    string sarg = argv[1];
    if ( sarg == "-h" ) {
      cout << argv[0] << ": [ChannelMapClass] [dorop]";
      cout << argv[0] << ": [Geometry/ChannelMapClass] [dorop]";
      return 0;
    }
    string::size_type ipos = sarg.find("/");
    if ( ipos == string ::npos ) {
      chanmap = sarg;
    } else {
      gname = sarg.substr(0, ipos);
      chanmap = sarg.substr(ipos + 1);
    }
  }
  if ( argc > 2 ) {
    string sarg = argv[2];
    dorop = sarg == "1" || sarg == "true";
  }
  if ( argc > 3 ) {
    istringstream ssarg(argv[3]);
    ssarg >> maxchanprint;
  }
  ExpectedValues ev;
  test_Geometry_Dune35t(gname, chanmap, dorop, maxchanprint);
  cout << "Tests concluded." << endl;
  ArtServiceHelper::close();
  return 0;
}

//**********************************************************************
