// test_GeometryProtoDune.cxx

// David Adams
// November 2016
//
// Test the protoDune geometry including channel mapping.

#include "test_GeometryDune.cxx"

//**********************************************************************

// Set expected values.

void setExpectedValues(ExpectedValues& ev) {
  ev.gname = "protodune_geo";
  ev.fullname = "protodune";
  // Geometry counts.
  ev.ncry = 1;
  ev.ntpc = 12;
  ev.npla = 3;
  ev.napa = ev.ntpc/2;
  ev.nrop = 4;
  // Signal type and view for each TPC plane.
  resize(ev.sigType, ev.ncry, ev.ntpc, ev.npla, geo::kMysteryType);
  resize(ev.view, ev.ncry, ev.ntpc, ev.npla, geo::kUnknown);
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
  ev.nwirPerPlane[0][0] = 1148;
  ev.nwirPerPlane[0][1] = 1148;
  ev.nwirPerPlane[0][2] =  480;
  for ( Index itpc=1; itpc<ev.ntpc; ++itpc ) {
    for ( Index ipla=0; ipla<ev.npla; ++ipla ) ev.nwirPerPlane[itpc][ipla] = ev.nwirPerPlane[0][ipla];
  }
  // Channel count per ROP.
  resize(ev.nchaPerRop, ev.nrop, 0);
  ev.nchaPerRop[0] = 800;
  ev.nchaPerRop[1] = 800;
  ev.nchaPerRop[2] = 480;
  ev.nchaPerRop[3] = 480;
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
  // Space points.
  #include "setProtoDuneSpacePoints.dat"
  // Optical detectors.
  // ... # detectors
  ev.nopdet = 10*ev.napa;
  // ... # channels in each detector
  ev.nopdetcha.resize(ev.nopdet, 4);
  ev.nopdetcha[31] = 12;
  ev.nopdetcha[46] = 12;
  // ... First channel in each detector
  vector<Index> firstChan = {
      0, 144,   4, 148,   8, 152,  12, 156,  16, 160,
     20, 164,  24, 168,  28, 172,  32, 176,  36, 180,
     48, 240,  52, 244,  56, 248,  60, 252,  64, 256,
     68, 264,  72, 260,  76, 276,  80, 280,  84, 284,
     96, 216, 100, 220, 104, 224, 132, 228, 108, 192,
    112, 232, 116, 196, 120, 200, 124, 236, 128, 204
  };
  ev.opdetcha.resize(ev.nopdet);
  Index ncha = 0;
  // ... Channels in each detector
  for ( Index iodt=0; iodt<ev.nopdet; ++iodt ) {
    ev.opdetcha[iodt].resize(ev.nopdetcha[iodt]);
    icha = firstChan[iodt];
    for ( Index ioch=0; ioch<ev.nopdetcha[iodt]; ++ioch ) {
      ev.opdetcha[iodt][ioch] = icha++;
      ++ncha;
    }
  }
  ev.nopcha = ncha;
}

//**********************************************************************

void setExpectedValuesSpacePoints(Geometry* pgeo) {
  const string myname = "setExpectedValuesSpacePoints: ";
  string ofname = "setExpectedValuesSpacePoints.dat";
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

