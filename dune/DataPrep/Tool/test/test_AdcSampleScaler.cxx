// test_AdcSampleScaler.cxx
//
// David Adams
// March 2019
//
// Test AdcSampleScaler.

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

int test_AdcSampleScaler(bool useExistingFcl, float qsig, float ped, float noise) {
  const string myname = "test_AdcSampleScaler: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_AdcSampleScaler.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: AdcSampleScaler" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    ScaleFactor: 100" << endl;
    fout << "    OutputUnit: \"e\"" << endl;
    fout << "    InputUnit: \"ADC\"" << endl;
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
  Index nsam = 50;
  Vector qsams(nsam, 0.0);
  Vector qexps(nsam, 99.);
  for ( Index isam=0; isam<nsam; ++isam ) {
    float qsam = 200 + isam + 0.1*(isam%10);
    qsams[isam] = qsam;
    qexps[isam] = 100*qsam;
  }
  AdcChannelData acd;
  acd.run = 123;
  acd.event = 456;
  acd.channel = 12345;
  acd.samples = qsams;
  acd.sampleUnit = "ADC";
  assert ( qsams.size() == nsam );

  cout << myname << line << endl;
  cout << myname << "Call tool." << endl;
  DataMap ret = ptoo->update(acd);
  assert( ! ret );
  assert( ret.getInt("ascUnitCheck") == 0 );

  cout << myname << line << endl;
  cout << myname << "Check samples." << endl;
  for ( Index isam=0; isam<nsam; ++isam ) {
    float qsam = acd.samples[isam];
    float qexp = qexps[isam];
    float qdif = fabs(qsam - qexp);
    assert( qdif < 1.e-5 );
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  float qsig = 500.0;
  float ped = 2.0;
  float noise = 0.0;
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
    string sarg(argv[2]);
    istringstream ssarg(sarg);
    ssarg >> qsig;
  }
  if ( argc > 3 ) {
    string sarg(argv[3]);
    istringstream ssarg(sarg);
    ssarg >> ped;
  }
  if ( argc > 4 ) {
    string sarg(argv[4]);
    istringstream ssarg(sarg);
    ssarg >> noise;
  }
  return test_AdcSampleScaler(useExistingFcl, qsig, ped, noise);
}

//**********************************************************************
