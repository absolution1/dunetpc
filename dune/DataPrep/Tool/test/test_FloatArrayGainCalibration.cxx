// test_FloatArrayGainCalibration.cxx
//
// David Adams
// November 2017
//
// Test FloatArrayGainCalibration.

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

//**********************************************************************

int test_FloatArrayGainCalibration(bool useExistingFcl =false) {
  const string myname = "test_FloatArrayGainCalibration: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_FloatArrayGainCalibration.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  calvals: {" << endl;
    fout << "       tool_type: FclFloatArray" << endl;
    fout << "        LogLevel: 2" << endl;
    fout << "    DefaultValue: -1" << endl;
    fout << "          Offset: 0" << endl;
    fout << "           Label: myvals" << endl;
    fout << "            Unit: ke" << endl;
    fout << "          Values: [ 0.11, 0.12, 0.13, 0.14, 0.15 ]" << endl;
    fout << "  }" << endl;

    fout << "  mytool: {" << endl;
    fout << "    tool_type: FloatArrayGainCalibration" << endl;
    fout << "    LogLevel: 2" << endl;
    fout << "    Unit: fC" << endl;
    fout << "    GainDefault: 1.0" << endl;
    fout << "    AdcUnderflowDefault: 0" << endl;
    fout << "    AdcOverflowDefault: 255" << endl;
    fout << "    GainTool: calvals" << endl;
    fout << "    ScaleFactor: \"[gain]/14.0\"" << endl;
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
  assert( tm.toolNames().size() == 2 );

  cout << myname << line << endl;
  cout << myname << "Fetching tool." << endl;
  auto ptoo = tm.getPrivate<TpcDataTool>("mytool");
  assert( ptoo != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create data." << endl;
  AdcChannelData acd0;
  // Specify the true signal (same for all channels).
  AdcSignalVector sigvals= {20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
  AdcIndex nsam = sigvals.size();
  assert( nsam == 10 );
  // Build data for five channels with different pedestal and gain for each.
  AdcSignal ped = 10.0;
  AdcSignal dped = 1.0;
  AdcChannelDataMap acds;
  for ( AdcChannel icha=0; icha<5; ++icha ) {
    float gain = 0.1 + 0.01*(icha+1);
    AdcChannelData& acd = acds[icha];
    acd.setChannelInfo(icha);
    acd.pedestal = ped;
    for ( AdcSignal sigval : sigvals ) {
      AdcIndex adc = sigval/gain + ped;
      acd.raw.push_back(adc);
    }
    ped += dped;
  }

  cout << myname << line << endl;
  cout << myname << "Calibrate and check data." << endl;
  //int w = 8;
  for ( auto& ient : acds ) {
    AdcChannel icha = ient.first;
    cout << myname << "-------------------- Channel " << icha << endl;
    AdcChannelData& acd = ient.second;
    DataMap res = ptoo->update(acd);
    res.print();
    assert( res.status() == 0 );
    assert( res.getInt("calibSampleCount") == int(nsam) );
    assert( acd.channel() == icha );
    assert( acd.raw.size() == nsam );
    assert( acd.samples.size() == nsam );
    //cout.precision(2);
    cout << "       ADC     signal  flag" << endl;
    for ( Index isam=0; isam<nsam; ++isam ) {
      cout << setw(3) << isam << ": " << setw(5) << acd.raw[isam]
           << setw(11) << acd.samples[isam]
           << setw(6) << acd.flags[isam] << endl;
      assert( fabs(acd.samples[isam] - sigvals[isam]) < 0.5 );
    }
  }

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
  return test_FloatArrayGainCalibration(useExistingFcl);
}

//**********************************************************************
