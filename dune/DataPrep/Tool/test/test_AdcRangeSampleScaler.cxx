// test_AdcRangeSampleScaler.cxx
//
// David Adams
// March 2019
//
// Test AdcRangeSampleScaler.

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "TRandom.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::istringstream;
using std::setw;
using fhicl::ParameterSet;

using Vector = std::vector<float>;
using Index = unsigned int;

//**********************************************************************

int test_AdcRangeSampleScaler(bool useExistingFcl, bool useMod) {
  const string myname = "test_AdcRangeSampleScaler: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcRangeSampleScaler.fcl";
  Index chaMod = useMod ? 30 : 0;
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: AdcRangeSampleScaler" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    RangeLimits: [10, 20]" << endl;
    fout << "    RangeModulus: " << chaMod << endl;
    fout << "    ScaleFactors: [2, 4, 8]" << endl;
    fout << "  }" << endl;
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
  assert( tm.toolNames().size() == 1 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto ptoo = tm.getPrivate<AdcChannelTool>("mytool");
  assert( ptoo != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  Index ncha = 50;
  Index nsam = 20;
  AdcChannelDataMap acds;
  for ( Index icha=0; icha<ncha; ++icha ) {
    AdcChannelData& acd = acds[icha];
    acd.setEventInfo(123, 456);
    acd.setChannelInfo(icha);
    acd.samples.resize(nsam, 10.0);
    acd.sampleUnit = "ADC";
  }

  cout << myname << line << endl;
  cout << myname << "Create expected data." << endl;
  Vector qexps(ncha, 10.);
  for ( Index icha=0; icha<ncha; ++icha ) {
    Index ichaMod = useMod ? icha%chaMod : icha;
    if ( ichaMod >= 20 )      qexps[icha] = 80;
    else if ( ichaMod >= 10 ) qexps[icha] = 40;
    else                      qexps[icha] = 20;
  }

  cout << myname << line << endl;
  cout << myname << "Call tool." << endl;
  DataMap ret = ptoo->updateMap(acds);
  assert( ! ret );

  cout << myname << line << endl;
  cout << myname << "Check samples." << endl;
  for ( const AdcChannelDataMap::value_type& iacd : acds ) {
    Index icha = iacd.first;
    const AdcChannelData& acd = iacd.second;
    assert( acd.channel() == icha );
    assert( acd.samples.size() == nsam );
    float qsam = acd.samples[0];
    float qexp = qexps[icha];
    float qdif = fabs(qsam - qexp);
    cout << myname << setw(4) << icha << ":" << setw(10) << qexp << setw(10) << qsam << endl;
    assert( qdif < 1.e-5 );
    for ( float qchk : acd.samples ) {
      assert( qchk == qsam );
    }
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  bool useMod = false;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG [qsig [ped [noise]]]]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  if ( argc > 2 ) {
    useMod = true;
  }
  return test_AdcRangeSampleScaler(useExistingFcl, useMod);
}

//**********************************************************************
