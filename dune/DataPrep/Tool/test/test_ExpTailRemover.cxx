// test_ExpTailRemover.cxx
//
// David Adams
// Marxh 2019
//
// Test ExpTailRemover.

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

int test_ExpTailRemover(bool useExistingFcl, float qsig, float ped, float noise) {
  const string myname = "test_ExpTailRemover: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_ExpTailRemover.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mytool: {" << endl;
    fout << "    tool_type: ExpTailRemover" << endl;
    fout << "    LogLevel: 1" << endl;
    fout << "    SignalUnit: \"ke\"" << endl;
    fout << "    DecayTime: 10.0" << endl;
    fout << "    SignalThreshold: 10.0" << endl;
    fout << "    CorrectFlag: []" << endl;
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
  cout << myname << "Create data and call tool." << endl;
  AdcChannelData acd;
  Index nsam = 50;
  Vector qsigs(nsam, 0.0);
  assert ( qsigs.size() == nsam );
  Index isig = 5;
  qsigs[isig] += qsig;
  Vector qdats = qsigs;
  // Add tail.
  double tdec =10.0;
  double alpha = 1.0/tdec;
  double beta = exp(-alpha);
  double qtai = -alpha*sqrt(beta)*qsig;   // Tail in 1st bin after signal.
  for ( Index isam=isig+1; isam<nsam; ++isam ) {
    qdats[isam] += qtai;
    qtai *= beta;
  }
  // Add pedestal.
  for ( Index isam=0; isam<nsam; ++isam ) {
    qdats[isam] += ped;
  }
  // Add noise.
  if ( noise != 0.0 ) {
    for ( Index isam=0; isam<nsam; ++isam ) {
      qdats[isam] += gRandom->Gaus(0.0, noise);
    }
  }
  // Display data.
  cout << myname << "Data:" << endl;
  for ( Index isam=0; isam<nsam; ++isam ) {
     cout << myname << setw(4) << isam << ": " << qdats[isam] << endl;
  }
  acd.run = 123;
  acd.event = 456;
  acd.channel = 12345;
  acd.sampleUnit = "ke";
  acd.samples = qdats;

  cout << myname << "Apply tool." << endl;
  DataMap ret = ptoo->update(acd);
  ret.print(myname);

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
  return test_ExpTailRemover(useExistingFcl, qsig, ped, noise);
}

//**********************************************************************
