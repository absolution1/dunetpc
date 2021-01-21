// test_Adc2dConvolute.cxx
//
// David Adams
// April 2019
//
// Test Adc2dConvolute.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "dune/DuneInterface/Tool/AdcChannelTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include <TRandom.h>
#include <TH1F.h>
#include <TCanvas.h>
#include "TLegend.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::istringstream;
using std::ofstream;
using fhicl::ParameterSet;
using std::vector;
using std::setw;
using std::fixed;

using Index = unsigned int;
using Name = std::string;

//**********************************************************************
namespace {

class Plot {
public:
  TPadManipulator topman;
  float ymax = 42.0;
  Plot(Index npad) : topman(400, 800) {
    topman.split(1, npad);
  }
  int addData(const AdcChannelDataMap& data, int col, int lwid, int lsty) {
    string myname = "Plot::addData: ";
    Index iman = 0;
    Name dopt = "HIST";
    for ( auto& ent : data ) {
      TPadManipulator* pman = topman.man(iman);
      if ( pman == nullptr ) {
        cout << myname << "ERROR: Unable to find subpad " << iman << endl;
        return 2;
      }
      ++iman;
      Index icha = ent.first;
      string scha = std::to_string(icha);
      const AdcSignalVector& sams = ent.second.samples;
      Index nsam = sams.size();
      //cout << myname << "Sample count is " << nsam << " for channel " << icha << endl;
      if ( nsam > 0 ) {
        Name hnam = "hch" + scha;
        Name httl = "Channel " + scha + ";Tick;Charge";
        TH1* ph = new TH1F(hnam.c_str(), httl.c_str(), nsam, 0, nsam);
        ph->SetStats(0);
        ph->SetLineWidth(lwid);
        ph->SetLineStyle(lsty);
        ph->SetLineColor(col);
        for ( Index isam=0; isam<nsam; ++isam ) {
          //cout << myname << "... " << setw(3) << isam << ": " << sams[isam] << endl;
          ph->Fill(isam+1, sams[isam]);
        }
        pman->add(ph, dopt);
        pman->setGridY();
        pman->setRangeY(0, ymax);
      }
    }
    //pcan->cd();
    //ph->Draw(dopt.c_str());
    return 0;
  }
  int print(Name fname) {
    string myname = "Plot::print: ";
    topman.draw();
    //cout << myname << "Printing " << fname << endl;
    topman.print(fname);
    return 0;
  }
};

bool checkEqual(float x1, float x2) {
  if ( x1 == x2 ) return true;
  if ( x1 > 0.0 && x2 > 0.0 ) {
    return fabs(x1-x2)/x1+x2 < 1.e-6;
  }
  return false;
}

}  // end unnamed namepsace
//**********************************************************************

int test_Adc2dConvolute(bool useExistingFcl) {
  const string mypre = "test_Adc2dConvolute";
  const string myname = mypre + ": ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_Adc2dConvolute.fcl";
  // Response has four bins/channels and is defined for central and first neightbor.
  // The summed response to a unit input bin charge is 1.0 in the central for any bin
  // and {0.5, 0.25, 0.125, 0.0625} in the neighbor from closest to furthest bin.
  //
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "tools: {" << endl;
    fout << "  mycon: {" << endl;
    fout << "           tool_type: Adc2dConvolute" << endl;
    fout << "            LogLevel: 3" << endl;
    fout << "     ResponseVectors: [" << endl;
    fout << "  [0.2, 0.4, 0.3, 0.1], [0.2, 0.4, 0.3, 0.1]," << endl;
    fout << "  [0.1, 0.2, 0.15, 0.05], [0.05, 0.10, 0.075, 0.025]," << endl;
    fout << "  [0.025, 0.050, 0.0375, 0.0125], [0.0125, 0.025, 0.01875, 0.00625]" << endl;
    fout << "                      ]" << endl;
    fout << "     ResponseCenter: 1" << endl;;
    fout << "        BinsPerWire: 4" << endl;;
    fout << "        BinsPerTick: 1" << endl;;
    fout << "         MinChannel: 101" << endl;;
    fout << "         MaxChannel: 104" << endl;;
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
  auto pcondir = tm.getPrivate<AdcChannelTool>("mycon");
  assert( pcondir != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create empty input data." << endl;
  Index nsam = 100;
  AdcChannelDataMap data;
  AdcSignalVector empty(nsam, 0.0);
  for ( Index icha=100; icha<105; ++icha ) {
    AdcChannelData& acd = data[icha];
    acd.setEventInfo(123, 456);
    acd.setChannelInfo(icha);
    acd.binSamples.resize(4, empty);
  }
  std::map<Index,float> expAreas;
  for ( Index icha=100; icha<106; ++icha ) expAreas[icha] = 0.0;

  Index nsig = 3;
  Index isig = 0;

  if ( isig < nsig ) {
    cout << myname << "Add signal " << ++isig << endl;
    data[103].binSamples[3][10] = 100.0;
    expAreas[102] +=  6.25;
    expAreas[103] += 100.0;
    expAreas[104] += 50.0;
  }

  if ( isig < nsig ) {
    cout << myname << "Add signal " << ++isig << endl;
    data[102].binSamples[1][20] = 100.0;
    expAreas[101] +=  25.0;
    expAreas[102] += 100.0;
    expAreas[103] +=  12.5;
  }

  cout << myname << line << endl;
  cout << myname << "Check signals." << endl;
  if ( isig < nsig ) {
    cout << myname << "Add signal " << ++isig << endl;
    data[101].binSamples[0][30] = 100.0;
    expAreas[101] += 100.0;
    expAreas[102] +=   6.25;
  }

  cout << myname << line << endl;
  cout << myname << "Check signals." << endl;
  cout << myname << "Signal count: " << nsig << endl;
  assert( isig == nsig );

  cout << myname << line << endl;
  cout << myname << "Call tool." << endl;
  DataMap ret = pcondir->updateMap(data);
  ret.print();
  assert( ret == 0 );
  assert( ret.getInt("responseVectorCount") == 6 );
  assert( ret.getInt("channelMin") == 101 );
  assert( ret.getInt("channelMax") == 104 );

  cout << myname << line << endl;
  cout << myname << "Check areas" << endl;
  cout << myname << "Chan  Nsam         Area" << endl;
  cout << myname << "----  ----  -----------" << endl;
  for ( const auto& kdat : data ) {
    Index icha = kdat.first;
    const AdcSignalVector& sams = kdat.second.samples;
    //assert( sams.size() == nsam );
    float area = 0.0;
    for ( float sam : sams ) area += sam;
    cout << myname << setw(4) << icha << setw(6) << nsam << setw(13) << area << endl;
    assert ( checkEqual(area, expAreas[icha]) );
  }

  cout << myname << line << endl;
  string fnam = "convhist.png";
  cout << myname << "Create plot " << fnam << endl;
  Plot plt(data.size());
  plt.addData(data, 1, 2, 1);
  plt.print(fnam);

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
      cout << "Usage: " << argv[0] << " [KEEPFCL [NOISE [SIGMAFIL [NSAM [SETSEED]]]]]" << endl;
      cout << "  KEEPFCL = if 1 (or true), existing FCL file is used [false]" << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_Adc2dConvolute(useExistingFcl);
}

//**********************************************************************
