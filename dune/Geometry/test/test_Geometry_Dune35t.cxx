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
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using geo::CryostatID;
typedef readout::TPCsetID APAID;
using readout::ROPID;

typedef unsigned int Index;

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

int test_Geometry_Dune35t(string chanmap ="Dune35tChannelMapAlg", bool dorop =true) {
  const string myname = "test_Geometry_Dune35t: ";
  string gname = "dune35t_geo";
  cout << myname << "Starting test" << endl;
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create configuration." << endl;
  const char* ofname = "test_Geometry_Dune35t.fcl";
  {
    ofstream fout(ofname);
    fout << "#include \"geometry_dune.fcl\"" << endl;;
    fout << "services.Geometry:                   @local::" + gname << endl;;
    fout << "services.ExptGeoHelperInterface:     @local::dune_geometry_helper" << endl;;
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
  check("Geometry name", pgeo->DetectorName(), "dune35t4apa_v6");
  cout << myname << "ROOT name: " << pgeo->ROOTFile() << endl;
  cout << myname << "GDML name: " << pgeo->GDMLFile() << endl;

  cout << myname << line << endl;
  if ( 1 ) {
    cout << myname << "World box"  << endl;
    double xlo, ylo, zlo;
    double xhi, yhi, zhi;
    pgeo->WorldBox(&xlo, &ylo, &zlo, &xhi, &yhi, &zhi);
    check("  xlo", xlo);
    check("  ylo", ylo);
    check("  zlo", zlo);
    check("  xhi", xhi);
    check("  yhi", yhi);
    check("  zhi", zhi);
  }

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
  Index encry = 1;
  Index entpc = 8;
  Index enpla = 3;
  Index enapa = 4;
  Index enchatot = 2048;
  Index ncry = pgeo->Ncryostats();
  check("Ncryostats", ncry, encry);
  check("MaxTPCs", pgeo->MaxTPCs(), entpc);
  check("MaxPlanes", pgeo->MaxPlanes(), enpla);
  check("TotalNTPC", pgeo->TotalNTPC(), entpc);
  check("Nviews", pgeo->Nviews(), enpla);
  check("Nchannels", pgeo->Nchannels(), enchatot);

  cout << myname << line << endl;
  cout << "TPC wire plane counts." << endl;
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
      }
    }
  }

  cout << myname << line << endl;
  //check("NROPs", pgeo->NROPs(), 16);
  Index enrop = 4;
  Index encha[enrop] = {144, 144, 112, 112};
  if ( dorop ) {
    cout << "APAs and ROPs." << endl;
    check("MaxROPs", pgeo->MaxROPs(), enrop);
    Index icry = 0;
    for ( CryostatID cryid: pgeo->IterateCryostatIDs() ) {
      Index napa = pgeo->NTPCsets(cryid);
      cout << "  Cryostat " << icry << " has " << napa << " APAs" << endl;
      assert( napa == enapa );
      for ( Index iapa=0; iapa<napa; ++iapa ) {
        APAID apaid(cryid, iapa);
        Index nrop = pgeo->NROPs(apaid);
        cout << "    APA " << iapa << " has " << nrop << " ROPs" << endl;
        assert( nrop == enrop );
        for ( Index irop=0; irop<nrop; ++irop ) {
          ROPID ropid(apaid, irop);
          Index ncha = pgeo->Nchannels(ropid);
          Index icha1 = pgeo->FirstChannelInROP(ropid);
          Index icha2 = icha1 + ncha - 1;
          cout << "      ROP " << irop << " has " << ncha << " channels: ["
               << icha1 << ", " << icha2 << "]" << endl;
          assert( ncha = encha[irop] );
        }
      }
      ++icry;
    }
    assert( icry == ncry );
  } else {
    cout << "Skipped APA and ROP tests." << endl;
  }  // end dorop

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

int main(int argc, const char* argv[]) {
  string chanmap = "";
  bool dorop = false;
  if ( argc > 1 ) {
    string sarg = argv[1];
    if ( sarg == "-h" ) {
      cout << argv[0] << ": ChannelMapClass []";
      return 0;
    }
    chanmap = sarg;
  }
  if ( argc > 2 ) {
    string sarg = argv[2];
    dorop = sarg == "1" || sarg == "true";
  }
  test_Geometry_Dune35t(chanmap, dorop);
  return 0;
}
