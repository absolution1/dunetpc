// test_AdcChannelTool.cxx
//
// David Adams
// April 2018
//
// Test AdcChannelTool.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
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

AdcChannelDataMap makeAdcData(Index ncha =10) {
  AdcChannelDataMap acds;
  for ( Index icha=0; icha<10; ++icha ) {
    acds[icha].channel = icha;
    acds[icha].fembID = icha%4;
  }
  return acds;
}

//**********************************************************************

// Tool that modfies data.
// Sets femb = 100+channel.

class AdcChannelTool_update : public AdcChannelTool {
public:
  DataMap update(AdcChannelData& acd) const override;
  bool viewWithUpdate() const override { return true; }
};
  
DataMap AdcChannelTool_update::update(AdcChannelData& acd) const {
  cout << "AdcChannelTool_update::update: Modifying channel " << acd.channel << endl;
  acd.fembID = 100 + acd.channel;
  int fembchan = 1000*acd.fembID + acd.channel;
  std::vector<int> fembchans(1, fembchan);
  DataMap ret;
  ret.setInt("fembchan", fembchan);
  ret.setIntVector("fembchans", fembchans);
  return ret;
}

//**********************************************************************

// Test with all default methods.

int test_AdcChannelTool_default() {
  const string myname = "test_AdcChannelTool_default: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << endl;
  cout << myname << line << endl;
  cout << myname << "Instantiate tool." << endl;
  AdcChannelTool act;
  Index ncha = 10;
  AdcChannelDataMap acds = makeAdcData(ncha);
  AdcChannelData acd = acds[1];

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
  cout << myname << "Call updateMap." << endl;
  ret = act.updateMap(acds);
  assert( ret == act.interfaceNotImplemented() );
  ret.print();

  cout << myname << line << endl;
  cout << myname << "Call viewMap." << endl;
  ret = act.viewMap(acds);
  assert( ret == act.interfaceNotImplemented() );
  ret.print();

  cout << myname << line << endl;
  cout << myname << "Test complete." << endl;

  return 0;
}

//**********************************************************************

// Test with update.

int test_AdcChannelTool_update() {
  const string myname = "test_AdcChannelTool_update: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << endl;
  cout << myname << line << endl;
  cout << myname << "Instantiate tool." << endl;
  AdcChannelTool_update act;
  Index ncha = 10;
  AdcChannelDataMap acds = makeAdcData(ncha);
  AdcChannelData acd1 = acds[1];

  DataMap ret(111);

  cout << myname << line << endl;
  cout << myname << "Call update." << endl;
  acd1.channel = 1;
  acd1.fembID = -1;
  ret = act.update(acd1);
  ret.print();
  assert( ret == 0 );
  assert( acd1.channel == 1 );
  assert( acd1.fembID == 101 );

  cout << myname << line << endl;
  cout << myname << "Call view." << endl;
  acd1.channel = 1;
  acd1.fembID = 200;
  ret = act.view(acd1);
  ret.print();
  assert( ! ret );
  assert( acd1.channel == 1 );
  assert( acd1.fembID == 200 );

  cout << myname << line << endl;
  cout << myname << "Call updateMap." << endl;
  ret = act.updateMap(acds);
  ret.print();
  assert( ! ret );
  for ( const auto& iacd : acds ) {
    Index icha = iacd.first;
    const AdcChannelData& acd = iacd.second;
    cout << myname << "  icha, femb = " << acd.channel << ", " << acd.fembID << endl;
    assert(acd.channel == icha);
    assert(acd.fembID == 100 + acd.channel);
  }

  cout << myname << line << endl;
  cout << myname << "Call viewMap." << endl;
  for ( auto& iacd : acds ) iacd.second.fembID = 200;
  ret = act.viewMap(acds);
  ret.print();
  assert( ! ret );
  for ( const auto& iacd : acds ) {
    Index icha = iacd.first;
    const AdcChannelData& acd = iacd.second;
    cout << myname << "  icha, femb = " << acd.channel << ", " << acd.fembID << endl;
    assert(acd.channel == icha);
    assert(acd.fembID == 200);
  }
  cout << myname << line << endl;
  cout << myname << "Test complete." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  test_AdcChannelTool_default();
  test_AdcChannelTool_update();
  return 0;
}

//**********************************************************************
