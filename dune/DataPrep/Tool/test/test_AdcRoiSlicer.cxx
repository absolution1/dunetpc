// test_AdcRoiSlicer.cxx
//
// David Adams
// September 2019
//
// Test AdcRoiSlicer.

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/ArtSupport/DuneToolManager.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::setw;
using fhicl::ParameterSet;

using Index = unsigned int;
using IndexVector = std::vector<Index>;
using Name = string;
using NameVector = std::vector<Name>;

bool checkint(int i1, int i2) {
  const Name myname = "checkint: ";
  if ( i1 == i2 ) return true;
  cout << myname << i1 << " != " << i2 << endl;
  return false;
}

//**********************************************************************

int test_AdcRoiSlicer(int opt, Index nsam, bool useExistingFcl =false) {
  const string myname = "test_AdcRoiSlicer: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcRoiSlicer.fcl";
  NameVector tnams = {"", "mytool1", "mytool2", "mytool3"};
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    for ( int iopt : {1, 2, 3} ) {
      fout << "  mytool" <<iopt << ": {" << endl;
      fout << "    tool_type: AdcRoiSlicer" << endl;
      fout << "    LogLevel: 1" << endl;
      fout << "    OutViewName: rois" << endl;
      fout << "    SliceOpt: " << iopt << endl;
      fout << "     CopyRaw: false" << endl;
      fout << "  }" << endl;
    }
    fout << "}" << endl;
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  cout << myname << line << endl;
  cout << myname << "Fetching tool manager." << endl;
  DuneToolManager* ptm = DuneToolManager::instance(fclfile);
  assert ( ptm != nullptr );
  DuneToolManager& tm = *ptm;
  tm.print();
  assert( tm.toolNames().size() == 3 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto psli = tm.getPrivate<TpcDataTool>(tnams[opt]);
  assert( psli != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data and call tool." << endl;
  AdcChannelData data;
  bool endInSig = false;
  Index isam1 = 20;
  Index isam2 = 40;
  Index isam3 = 80;
  Index isam4 = 90;
  for ( AdcIndex isam=0; isam<nsam; ++isam ) {
    float xadc = rand()%20 - 10.0;
    bool isSig = (isam >= isam1 && isam < isam2) || (isam >= isam3 && isam < isam4);
    if ( isSig ) xadc += 100;
    data.signal.push_back(isSig);
    data.samples.push_back(xadc);
    endInSig = isSig;
    cout << myname << setw(4) << isam << setw(3) << data.signal[isam]
         << " " << data.samples[isam] << endl;
  }
  data.roisFromSignal();
  Index nroi = data.rois.size();
  Index nnot = nroi;
  if ( ! endInSig ) ++nnot;
  cout << myname << " opt: " <<  opt << endl;
  cout << myname << "nroi: " << nroi << endl;
  cout << myname << "nnot: " << nnot << endl;
  Index nexp = 0;
  IndexVector viewSizes;
  IndexVector viewTick0s;
  if ( nsam < isam4 ) isam4 = nsam;
  if ( nsam < isam3 ) isam3 = nsam;
  if ( opt == 1 ) {
    nexp = nroi;
    viewTick0s.push_back(isam1);
    viewTick0s.push_back(isam3);
    viewSizes.push_back(isam2 - isam1);
    viewSizes.push_back(isam4 - isam3);
  } else if ( opt == 2 ) {
    nexp = nnot;
    viewTick0s.push_back(0);
    viewTick0s.push_back(isam2);
    viewTick0s.push_back(isam4);
    viewSizes.push_back(isam1);
    viewSizes.push_back(isam3 - isam2);
    viewSizes.push_back(nsam - isam4);
  } else if ( opt == 3 ) {
    nexp = nroi + nnot;
    viewTick0s.push_back(0);
    viewTick0s.push_back(isam1);
    viewTick0s.push_back(isam2);
    viewTick0s.push_back(isam3);
    viewTick0s.push_back(isam4);
    viewSizes.push_back(isam1);
    viewSizes.push_back(isam2 - isam1);
    viewSizes.push_back(isam3 - isam2);
    viewSizes.push_back(isam4 - isam3);
    viewSizes.push_back(nsam - isam4);
  } else {
    assert(false);
  }
  cout << myname << "nexp: " << nexp << endl;
  assert( data.signal.size() == nsam );
  assert( data.samples.size() == nsam );

  cout << myname << line << endl;
  cout << myname << "Running tool." << endl;
  DataMap res = psli->update(data);
  res.print();

  cout << myname << line << endl;
  cout << myname << "Checking results." << endl;
  assert( res == 0 );
  assert( Index(res.getInt("nRoiView")) == nexp );
  assert( data.viewSize() == 1 );
  assert( data.viewNames().size() == 1 );
  Name vnam = "rois";
  assert( data.viewNames()[0] == vnam );
  assert( data.view(vnam).size() == nexp );
  for ( Index ivie=0; ivie<nexp; ++ivie ) {
    cout << myname << "..entry " << ivie << endl;
    const AdcChannelData& acd = data.view(vnam)[ivie];
    assert( checkint(acd.tick0, int(viewTick0s[ivie])) );
    assert( acd.viewParent != nullptr );
    assert( acd.viewParent == &data );
    Index nvsam = viewSizes[ivie];
    assert( checkint(acd.samples.size(), nvsam) );
    assert( checkint(acd.signal.size(), nvsam) );
    for ( Index ivsam=0; ivsam<nvsam; ++ivsam ) {
      int qv = 1000*acd.samples[ivsam];
      int q0 = 1000*data.samples[ivsam + acd.tick0];
      int sv = acd.signal[ivsam];
      int s0 = data.signal[ivsam + acd.tick0];
      cout << myname << setw(3) << ivsam << ": " << sv << " "
           << setw(7) << qv << endl;
      assert( checkint(qv, q0) );
      assert( checkint(sv, s0) );
    }
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  int rs = 0;
  if ( (rs = test_AdcRoiSlicer(1, 100, useExistingFcl)) ) return rs;
  if ( (rs = test_AdcRoiSlicer(2, 100, useExistingFcl)) ) return rs;
  if ( (rs = test_AdcRoiSlicer(3, 100, useExistingFcl)) ) return rs;
  if ( (rs = test_AdcRoiSlicer(1,  85, useExistingFcl)) ) return rs;
  return 0;
}

//**********************************************************************
