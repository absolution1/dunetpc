// test_Tpc2dDeconvolute.cxx
//
// David Adams
// April 2019
//
// Test Tpc2dDeconvolute.

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "dune/DuneInterface/Tool/TpcDataTool.h"
#include "dune/ArtSupport/DuneToolManager.h"
#include "dune/DuneCommon/Utility/LineColors.h"
#include "dune/DuneCommon/Utility/TPadManipulator.h"
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
  float ymin =  0.0;
  float ymax = 42.0;
  Name dopt = "HIST";
  Plot(Index npad) : topman(400, 800) {
    topman.split(1, npad);
  }
  void setRange(float a_ymin, float a_ymax) {
    ymin = a_ymin;
    ymax = a_ymax;
  }
  int addData(const TpcData& tpd, Name roinam, int col, int lwid, int lsty) {
    string myname = "Plot::addData: ";
    // Log messages
    //   0 - none
    //   1 - sample count for each channel
    //   2 - non-zero samples
    //   3 - all samples
    int dbg = 1;
    Index iman = 0;
    const AdcChannelDataMap& data = *tpd.getAdcData()[0];
    for ( auto& ent : data ) {
      const AdcSignalVector& sams = ent.second.samples;
      Index nsam = sams.size();
      Index icha = ent.first;
      string scha = std::to_string(icha);
      string slab = roinam.size() ? roinam : "ADC";
      if ( dbg > 0 ) cout << myname << slab << " sample count is " << nsam << " for channel "
                          << icha << endl;
      TPadManipulator* pman = topman.man(iman);
      if ( pman == nullptr ) {
        cout << myname << "ERROR: Unable to find subpad " << iman << endl;
        return 2;
      }
      ++iman;
      Name hnam = "hch" + scha;
      Name httl = "Channel " + scha + ";Tick;Charge";
      TH1* ph = new TH1F(hnam.c_str(), httl.c_str(), nsam, 0, nsam);
      ph->SetDirectory(nullptr);
      ph->SetStats(0);
      ph->SetLineColor(col);
      ph->SetLineWidth(lwid);
      ph->SetLineStyle(lsty);
      if ( roinam.size() ) {
        assert(tpd.getTpcData(roinam) != nullptr );
        const TpcData::Tpc2dRoiVector& rois = tpd.getTpcData(roinam)->get2dRois();
        assert( rois.size() == 1 );
        const Tpc2dRoi& roi = rois[0];
        assert( roi.sampleSize() == nsam );
        for ( Index isam=0; isam<nsam; ++isam ) {
          float val = roi.value(icha, isam);
          if ( dbg > 2 || (dbg > 1 && sams[isam]) ) {
            cout << myname << "... " << setw(3) << isam << ": " << val << endl;
          }
          ph->Fill(isam+1, val);
        }
      } else {
        for ( Index isam=0; isam<nsam; ++isam ) {
          if ( dbg > 2 || (dbg > 1 && sams[isam]) ) {
            cout << myname << "... " << setw(3) << isam << ": " << sams[isam] << endl;
          }
          ph->Fill(isam+1, sams[isam]);
        }
      }
      if ( nsam > 0 ) {
        pman->add(ph, dopt);
        pman->setGridY();
        pman->setRangeY(ymin, ymax);
      }
    }
    dopt = "HIST SAME";
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

int test_Tpc2dDeconvolute(bool useExistingFcl) {
  const string mypre = "test_Tpc2dDeconvolute";
  const string myname = mypre + ": ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_Tpc2dDeconvolute.fcl";
  // Response has four bins/channels and is defined for central and first neightbor.
  // The summed response to a unit input bin charge is 1.0 in the central for any bin
  // and {0.5, 0.25, 0.125, 0.0625} in the neighbor from closest to furthest bin.
  //
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "data.response: [" << endl;
    fout << "  [0.2, 0.4, 0.3, 0.1]," << endl;
    fout << "  [-0.04, -0.08, -0.06, -0.02]," << endl;
    fout << "  [0.04, 0.08, 0.06, 0.02]" << endl;
    fout << "]" << endl;
    fout << "tools: {" << endl;
    fout << "  mycon: {" << endl;
    fout << "           tool_type: Adc2dConvolute" << endl;
    fout << "            LogLevel: 1" << endl;
    fout << "     ResponseVectors: @local::data.response" << endl;
    fout << "     ResponseCenter: 1" << endl;;
    fout << "        BinsPerWire: 1" << endl;;
    fout << "        BinsPerTick: 0" << endl;;
    fout << "         MinChannel: 100" << endl;;
    fout << "         MaxChannel: 104" << endl;;
    fout << "  }" << endl;
    fout << "  myroi: {" << endl;
    fout << "           tool_type: AdcToRoi2d" << endl;
    fout << "            LogLevel: 1" << endl;
    fout << "            Option: 1" << endl;
    fout << "  }" << endl;
    fout << "  mydco: {" << endl;
    fout << "           tool_type: Tpc2dDeconvolute" << endl;
    fout << "            LogLevel: 1" << endl;
    fout << "             FftSize: 1000" << endl;
    fout << "     ResponseVectors: @local::data.response" << endl;
    fout << "     ResponseCenter: 1" << endl;
    fout << "             InPath: \"conv\"" << endl;
    fout << "            OutPath: \"dcon\"" << endl;
    fout << "        SampleSigma: 2.0" << endl;
    fout << "       ChannelSigma: 0.0" << endl;
    fout << "     LowFilterPower: 0.0" << endl;
    fout << "     LowFilterWidth: 0.0" << endl;
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
  assert( tm.toolNames().size() == 3 );

  cout << myname << line << endl;
  cout << myname << "Fetching tools." << endl;
  auto pcon = tm.getPrivate<TpcDataTool>("mycon");
  assert( pcon != nullptr );
  auto pdco = tm.getPrivate<TpcDataTool>("mydco");
  assert( pdco != nullptr );

  cout << myname << line << endl;
  cout << myname << "Create empty input data." << endl;
  Index nsam = 60;
  Index ncha = 5;
  Index chan0 = 100;
  TpcData tpd(1);
  AdcChannelDataMap& data = *tpd.getAdcData()[0];
  for ( Index icha=100; icha<105; ++icha ) {
    AdcChannelData& acd = data[icha];
    acd.setEventInfo(123, 456);
    acd.setChannelInfo(icha);
    acd.samples.resize(nsam, 0.0);
  }
  std::map<Index,float> expAreas;
  for ( Index icha=100; icha<106; ++icha ) expAreas[icha] = 0.0;

  Index nsig = 3;
  Index isig = 0;

  if ( isig < nsig ) {
    cout << myname << "Add signal " << ++isig << endl;
    float sfac = 1.0;
    data[102].samples[20] = 100.0*sfac;
    expAreas[100] += +20.0*sfac;
    expAreas[101] += -20.0*sfac;
    expAreas[102] += 100.0*sfac;
    expAreas[103] += -20.0*sfac;
    expAreas[104] += +20.0*sfac;
  }

  if ( isig < nsig ) {
    cout << myname << "Add signal " << ++isig << endl;
    data[103].samples[10] = 100.0;
    expAreas[101] +=  20.0;
    expAreas[102] += -20.0;
    expAreas[103] += 100.0;
    expAreas[104] += -20.0;
  }

  cout << myname << line << endl;
  cout << myname << "Check signals." << endl;
  if ( isig < nsig ) {
    cout << myname << "Add signal " << ++isig << endl;
    data[101].samples[30] = 100.0;
    expAreas[100] += -20.0;
    expAreas[101] += 100.0;
    expAreas[102] += -20.0;
    expAreas[103] +=  20.0;
  }

  cout << myname << "Check signals." << endl;
  cout << myname << "Signal count: " << nsig << endl;
  assert( isig == nsig );

  cout << myname << line << endl;
  string fnam = "noconvhist.png";
  cout << myname << "Create plot " << fnam << endl;
  Plot plt1(data.size());
  LineColors lc;
  plt1.addData(tpd, "", lc.green(), 2, 1);
  plt1.setRange(0, 105);
  plt1.print(fnam);

  cout << myname << line << endl;
  cout << myname << "Call convolution tool." << endl;
  DataMap ret = pcon->updateMap(data);
  ret.print();
  assert( ret == 0 );
  assert( ret.getInt("responseVectorCount") == 3 );
  assert( ret.getInt("channelMin") == 100 );
  assert( ret.getInt("channelMax") == 104 );

  cout << myname << line << endl;
  cout << myname << "Copy samples to ROI." << endl;
  tpd.addTpcData("conv");
  assert( tpd.getTpcData("conv") != nullptr );
  assert( tpd.getTpcData("dcon") == nullptr );
  assert( tpd.getTpcData("conv")->get2dRois().size() == 0 );
  tpd.getTpcData("conv")->get2dRois().emplace_back(5, nsam, 100, 0);
  assert( tpd.getTpcData("conv")->get2dRois().size() == 1 );
  Tpc2dRoi& roi = tpd.getTpcData("conv")->get2dRois().back();
  Tpc2dRoi::IndexArray idxs;
  Index& kcha = idxs[0];
  Index& ksam = idxs[1];
  for ( kcha=0; kcha<ncha; ++kcha ) {
    for ( ksam=0; ksam<nsam; ++ksam ) {
      roi.data().setValue(idxs, data[chan0+kcha].samples[ksam]);
    }
  }
  assert( tpd.getTpcData("conv")->get2dRois().size() == 1 );

  cout << myname << line << endl;
  fnam = "convhist.png";
  cout << myname << "Create plot " << fnam << endl;
  Plot plt2(data.size());
  plt2.setRange(-10, 45);
  plt2.addData(tpd,     "", lc.blue(), 2, 1);
  plt2.addData(tpd, "conv",  lc.red(), 2, 2);
  plt2.print(fnam);

  cout << myname << line << endl;
  cout << myname << "Check convoluted areas" << endl;
  cout << myname << "Chan  Nsam         Area" << endl;
  cout << myname << "----  ----  -----------" << endl;
  for ( const auto& kdat : data ) {
    Index icha = kdat.first;
    const AdcSignalVector& sams = kdat.second.samples;
    //assert( sams.size() == nsam );
    float area = 0.0;
    for ( float sam : sams ) area += sam;
    cout << myname << setw(4) << icha << setw(6) << nsam << setw(13) << area
         << " (" << expAreas[icha] << ")" << endl;
    assert ( checkEqual(area, expAreas[icha]) );
  }

  cout << myname << line << endl;
  cout << myname << "Call deconvolution tool." << endl;
  ret = pdco->updateTpcData(tpd);
  ret.print();
  assert( ret == 0 );
  assert( ret.getInt("dcoNroiIn") == 1 );
  assert( ret.getInt("dcoNroiOut") == 1 );
  assert( tpd.getTpcData("conv")->get2dRois().size() == 1 );
  assert( tpd.getTpcData("dcon")->get2dRois().size() == 1 );
  cout << myname << "Input normalization: "
       << tpd.getTpcData("conv")->get2dRois()[0].dft()->normalization().globalName() << ", "
       << tpd.getTpcData("conv")->get2dRois()[0].dft()->normalization().termName() << endl;
  cout << myname << "Output normaization: "
       << tpd.getTpcData("dcon")->get2dRois()[0].dft()->normalization().globalName() << ", "
       << tpd.getTpcData("dcon")->get2dRois()[0].dft()->normalization().termName() << endl;
  assert( tpd.getTpcData("conv")->get2dRois()[0].dft()->normalization() ==
          tpd.getTpcData("dcon")->get2dRois()[0].dft()->normalization() );

  cout << myname << line << endl;
  fnam = "dcovhist.png";
  cout << myname << "Create plot " << fnam << endl;
  Plot plt3(data.size());
  plt3.setRange(-10, 45);
  plt3.addData(tpd,     "",    lc.blue(), 2, 1);
  plt2.addData(tpd, "conv",     lc.red(), 2, 2);
  plt3.addData(tpd, "dcon",  lc.orange(), 2, 3);
  plt3.print(fnam);

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
  return test_Tpc2dDeconvolute(useExistingFcl);
}

//**********************************************************************
