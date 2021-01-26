// test_TpcDataTool.cxx
//
// David Adams
// January 2021
//
// Test TpcDataTool.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/ArtSupport/DuneToolManager.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::ofstream;
using fhicl::ParameterSet;
using std::vector;
using Index = unsigned int;

//**********************************************************************

void fillAdcData(AdcChannelDataMap* pacm, Index ncha =10) {
  for ( Index icha=0; icha<10; ++icha ) {
    (*pacm)[icha].setChannelInfo(icha, icha%4);
  }
}

//**********************************************************************

// Tool that modfies data.
// Sets femb = 100+channel.

class TpcDataTool_update : public TpcDataTool {
public:
  DataMap update(AdcChannelData& acd) const override;
  bool viewWithUpdate() const override { return true; }
};
  
DataMap TpcDataTool_update::update(AdcChannelData& acd) const {
  cout << "TpcDataTool_update::update: Modifying channel " << acd.channel() << endl;
  acd.setChannelInfo(acd.channel(), 100 + acd.channel());
  int fembchan = 1000*acd.fembID() + acd.channel();
  std::vector<int> fembchans(1, fembchan);
  DataMap ret;
  ret.setInt("fembchan", fembchan);
  ret.setIntVector("fembchans", fembchans);
  return ret;
}

//**********************************************************************

// Test with all default methods.

int test_TpcDataTool_default() {
  const string myname = "test_TpcDataTool_default: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << endl;
  cout << myname << line << endl;
  cout << myname << "Instantiate tool." << endl;
  TpcDataTool act;

  cout << endl;
  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  TpcData tpd;
  AdcChannelDataMap* pacm = tpd.createAdcData();
  assert( pacm->size() == 0 );
  Index ncha = 10;
  fillAdcData(pacm, ncha);
  assert( pacm->size() == ncha );
  AdcChannelData acd = (*pacm)[1];

  DataMap ret(111);

  cout << myname << line << endl;
  cout << myname << "Call update." << endl;
  ret = act.update(acd);
  assert( ret == act.interfaceNotImplemented() );
  ret.print();

  cout << myname << line << endl;
  cout << myname << "Call view." << endl;
  ret = act.view(acd);
  assert( ret == act.interfaceNotImplemented() );
  ret.print();

  cout << myname << line << endl;
  cout << myname << "Call update TPC." << endl;
  ret = act.updateTpcData(tpd);
  assert( ret == act.interfaceNotImplemented() );
  ret.print();

  cout << myname << line << endl;
  cout << myname << "Call view TPC." << endl;
  ret = act.viewTpcData(tpd);
  assert( ret == act.interfaceNotImplemented() );
  ret.print();

  cout << myname << line << endl;
  cout << myname << "Test complete." << endl;

  return 0;
}

//**********************************************************************

// Test with update.

int test_TpcDataTool_update() {
  const string myname = "test_TpcDataTool_update: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << endl;
  cout << myname << line << endl;
  cout << myname << "Instantiate tool." << endl;
  TpcDataTool_update act;

  cout << endl;
  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  TpcData tpd;
  AdcChannelDataMap* pacm = tpd.createAdcData();
  assert( pacm->size() == 0 );
  Index ncha = 10;
  fillAdcData(pacm, ncha);
  assert( pacm->size() == ncha );
  AdcChannelData acd = (*pacm)[1];

  cout << myname << line << endl;
  cout << myname << "Call update TPC." << endl;
  DataMap ret(111);
  ret = act.updateTpcData(tpd);
  ret.print();
  assert( ! ret );
  for ( const auto& iacd : *tpd.getAdcData() ) {
    Index icha = iacd.first;
    const AdcChannelData& acd = iacd.second;
    cout << myname << "  icha, femb = " << acd.channel() << ", " << acd.fembID() << endl;
    assert(acd.channel() == icha);
    assert(acd.fembID() == 100 + acd.channel());
  }

  cout << myname << line << endl;
  cout << myname << "Call view TPC." << endl;
  for ( auto& iacd : *tpd.getAdcData() ) iacd.second.setChannelInfo(iacd.second.channel(), 200);
  ret = act.viewTpcData(tpd);
  ret.print();
  assert( ! ret );
  for ( const auto& iacd : *tpd.getAdcData() ) {
    Index icha = iacd.first;
    const AdcChannelData& acd = iacd.second;
    cout << myname << "  icha, femb = " << acd.channel() << ", " << acd.fembID() << endl;
    assert(acd.channel() == icha);
    assert(acd.fembID() == 200);
  }
  cout << myname << line << endl;
  cout << myname << "Test complete." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  test_TpcDataTool_default();
  test_TpcDataTool_update();
  return 0;
}

//**********************************************************************
