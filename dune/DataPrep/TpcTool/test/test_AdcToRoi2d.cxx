// test_AdcToRoi2d.cxx
//
// David Adams
// February 2021
//
// Test AdcToRoi2d.

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TRandom.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::vector;
using std::setw;
using fhicl::ParameterSet;

using Index = unsigned int;

//**********************************************************************

int test_AdcToRoi2d(bool useExistingFcl =false) {
  const string myname = "test_AdcToRoi2d: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcToRoi2d.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"dataprep_tools.fcl\"" << endl;
    fout << "tools.mytool1: {" << endl;
    fout << "  tool_type: AdcToRoi2d" << endl;
    fout << "  LogLevel: 3" << endl;
    fout << "  Option: 1" << endl;
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
  assert( tm.toolNames().size() >= 2 );

  cout << myname << line << endl;
  cout << myname << "Create test data." << endl;
  // The data for each channel is offset by one tick using acd.tick0.
  TpcData tpd;
  TpcData::AdcDataPtr pacm = tpd.createAdcData();
  AdcChannelDataMap acds;
  Index ncha = 4;
  Index ntck = 5;
  Index ntckexp = ntck + ncha - 1;
  Index icha0 = 100;
  Index itck0 = 1000;
  Float2dData::IndexArray idxs;
  Tpc2dRoi exp(ncha, ntckexp, icha0, itck0);
  Index dtck = 0;
  DuneEventInfo* peviMutable = new DuneEventInfo(123, 1);
  peviMutable->triggerTick0 = itck0;
  AdcChannelData::EventInfoPtr pevi(peviMutable);
  for ( Index kcha=0; kcha<ncha; ++kcha, ++dtck ) {
    Index icha = icha0 + kcha;
    AdcChannelData& acd = (*pacm)[icha];
    acd.setEventInfo(pevi);
    acd.setChannelInfo(icha);
    acd.tick0 = dtck;
    acd.samples.resize(ntck, 0.0);
    float val = 10*(kcha + 1);
    idxs[0] = kcha;
    idxs[1] = dtck;
    for ( Index ktck=0; ktck<3; ++ktck, ++idxs[1] ) {
      exp.data().setValue(idxs, val);
      acd.samples[ktck] = val;
    }
  }
  assert( tpd.getAdcData().size() == 1 );
  assert( tpd.getAdcData()[0]->size() == ncha );
  for ( auto& iacd : *(tpd.getAdcData()[0]) ) {
    const AdcChannelData& acd = iacd.second;
    assert( acd.samples.size() == ntck );
  }
  assert( tpd.get2dRois().size() == 0 );
  
  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto ptoo = tm.getPrivate<TpcDataTool>("mytool1");
  assert( ptoo != nullptr );

  cout << myname << line << endl;
  cout << myname << "Call tool." << endl;
  DataMap res = ptoo->updateTpcData(tpd);
  res.print();
  assert( res.status() == 0 );

  cout << myname << line << endl;
  cout << myname << "Check results." << endl;
  assert( tpd.get2dRois().size() == 1 );
  const Tpc2dRoi& roi = tpd.get2dRois().at(0);
  assert( roi.channelSize() == ncha );
  assert( roi.channelOffset() == icha0 );
  assert( roi.sampleSize() == ntckexp );
  assert( roi.sampleOffset() == itck0 );
  Index& kcha = idxs[0];
  Index& ktck = idxs[1];
  for ( kcha=0; kcha<ncha; ++kcha ) {
    Index icha = icha0 + kcha;
    for ( ktck=0; ktck<ntckexp; ++ktck ) {
      Index itck = itck0 + ktck;
      float expval = exp.data().value(idxs);
      float val = roi.data().value(idxs);
      cout << myname << setw(5) << icha << setw(5) << itck << ": "
           << setw(12) << expval << " ?= " << setw(12) << val << endl;
      assert( val == expval );
    }
    cout << myname << line << endl;
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
  return test_AdcToRoi2d(useExistingFcl);
}

//**********************************************************************
