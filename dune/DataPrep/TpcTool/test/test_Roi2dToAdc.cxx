// test_Roi2dToAdc.cxx
//
// David Adams
// February 2021
//
// Test Roi2dToAdc.

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
using std::vector;
using std::setw;
using fhicl::ParameterSet;

using Index = unsigned int;

int showAdcMap(const AdcChannelDataMap& acm, string myname, const AdcChannelDataMap* pacmchk =0) {
  Index nerr = 0;
  // Display samples.
  cout << myname << "  chan  toff  samples" << endl;
  for ( auto& iacd : acm ) {
    const AdcChannelData& acd = iacd.second;
    Index icha = acd.channel();
    bool dochk = pacmchk;
    if ( dochk && pacmchk->count(icha) == 0 ) {
      ++nerr;
      dochk = false;
    }
    const AdcChannelData* pacdchk = dochk ? &pacmchk->find(icha)->second : nullptr;
    if ( dochk ) {
      if ( pacdchk->samples.size() != acd.samples.size() ) ++nerr;
    }
    cout << myname << setw(6) << acd.channel() << setw(6) << acd.tickOffset() << ":";
    bool first = true;
    Index isam = 0;
    for ( float sam : acd.samples ) {
      if ( first ) first = false;
      else cout << ",";
      if ( dochk ) {
        if ( pacdchk->samples.size() <= isam ||
             fabs(pacdchk->samples[isam] - sam ) > 1.e-5 ) ++nerr;
      }
      cout << setw(4) << sam;
      ++isam;
    }
    cout << " (";
    for ( bool sig : acd.signal ) cout << (sig ? "X" : ".");
    cout << ")" << endl;
  }
  // Display ROIs.
  cout << myname << "  chan  nROI  ROIs" << endl;
  for ( auto& iacd : acm ) {
    const AdcChannelData& acd = iacd.second;
    Index icha = acd.channel();
    cout << myname << setw(6) << icha << setw(6) << acd.rois.size() << ":";
    bool first = true;
    for ( const AdcRoi& roi : acd.rois ) {
      if ( first ) first = false;
      else cout << ",";
      cout << " " << roi.first << "-" << roi.second;
    }
    cout << endl;
  }
  return nerr;
}
      
int showRoi(const Tpc2dRoi& roi, string myname) {
  cout << myname << "ROI channel offset: " << roi.channelOffset() << endl;
  cout << myname << "ROI sample offset: " << roi.sampleOffset() << endl;
  for ( Index kcha=0; kcha<roi.channelSize(); ++kcha ) {
    Index icha = roi.channelOffset() + kcha;
    cout << myname << setw(6) << icha << ":";
    bool first = true;
    for ( Index ksam=0; ksam<roi.sampleSize(); ++ksam ) {
      if ( first ) first = false;
      else cout << ",";
      Index isam = roi.sampleOffset() + ksam;
      cout << setw(6) << roi.value(icha, isam, 8888);
    }
    cout << endl;
  }
  return 0;
}
      
//**********************************************************************

int test_Roi2dToAdc(bool useExistingFcl =false) {
  const string myname = "test_Roi2dToAdc: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_Roi2dToAdc.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "#include \"dataprep_tools.fcl\"" << endl;
    fout << "tools.mytool1: {" << endl;
    fout << "  tool_type: Roi2dToAdc" << endl;
    fout << "  LogLevel: 1" << endl;
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
  AdcChannelDataMap acmchk;
  Index ncha = 4;
  Index ntck = 10;
  Index ntckroi = 6;
  Index icha0 = 100;
  Index itck0 = 1000;
  Index tckoff = 2;   // Offset from tick0 to signal
  Float2dData::IndexArray idxs;
  assert( tpd.get2dRois().size() == 0 );
  tpd.get2dRois().emplace_back(ncha, ntckroi, icha0, itck0+tckoff);
  assert( tpd.get2dRois().size() == 1 );
  Tpc2dRoi& roi = tpd.get2dRois().back();
  Index dtck = 0;
  DuneEventInfo* peviMutable = new DuneEventInfo(123, 1);
  peviMutable->triggerTick0 = itck0;
  AdcChannelData::EventInfoPtr pevi(peviMutable);
  Index nsig = 0;
  for ( Index kcha=0; kcha<ncha; ++kcha, ++dtck ) {
    Index icha = icha0 + kcha;
    AdcChannelData& acd = (*pacm)[icha];
    acd.setEventInfo(pevi);
    acd.setChannelInfo(icha);
    acd.tick0 = dtck;
    acd.samples.resize(ntck, 999);
    float val = 10*(kcha + 1);
    idxs[0] = kcha;
    idxs[1] = dtck;
    AdcChannelData& acdchk = acmchk[icha];
    acdchk.setEventInfo(pevi);
    acdchk.setChannelInfo(icha);
    acdchk.tick0 = dtck;
    acdchk.samples.resize(ntck, 0.0);
    for ( Index ktck=tckoff; ktck<tckoff+3; ++ktck, ++idxs[1] ) {
      roi.data().setValue(idxs, val);
      acdchk.samples[ktck] = val;
      if ( acdchk.rois.size() == 0 ) {
        acdchk.rois.push_back({ktck, ktck});
      } else {
        acdchk.rois.back().second = ktck;
      }
    }
    acdchk.signal.resize(ntck, false);
    Index isig1 = tckoff > kcha ? tckoff - kcha : 0;
    Index isig2 = std::min(8 - kcha, ntck);
    for ( Index ktck=isig1; ktck<isig2; ++ktck ) {
      acdchk.signal[ktck] = true;
      ++nsig;
    }
  }
  assert( tpd.getAdcData().size() == 1 );
  assert( tpd.getAdcData()[0]->size() == ncha );
  for ( auto& iacd : *(tpd.getAdcData()[0]) ) {
    const AdcChannelData& acd = iacd.second;
    assert( acd.samples.size() == ntck );
  }
  assert( roi.channelSize() == ncha );
  assert( roi.channelOffset() == icha0 );
  assert( roi.sampleSize() == ntckroi );
  assert( roi.sampleOffset() == itck0 + tckoff );
  
  cout << myname << line << endl;
  cout << myname << "Check ROI data." << endl;
  assert( tpd.get2dRois().size() == 1 );
  assert( &tpd.get2dRois().front() == &roi );
  showRoi(roi, myname);
  cout << myname << line << endl;

  cout << myname << line << endl;
  cout << myname << "Check ADC data before tool:" << endl;
  assert( tpd.getAdcData().size() == 1 );
  assert( tpd.getAdcData().front() == pacm );
  AdcChannelDataMap& acm = *pacm;
  showAdcMap(acm, myname);
  
  cout << myname << line << endl;
  cout << myname << "Expected data after tool:" << endl;
  showAdcMap(acmchk, myname);

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
  cout << myname << "Check ADC data after tool." << endl;
  assert( tpd.getAdcData().size() == 1 );
  assert( tpd.getAdcData().front() == pacm );
  assert( showAdcMap(acm, myname, &acmchk) == 0 );
  assert( res.getInt("r2a_nchaZeroed") == int(ncha) );
  assert( res.getInt("r2a_nchaFilled") == int(ncha) );
  assert( res.getInt("r2a_nsamZeroed") == int(ncha*ntck) );
  assert( res.getInt("r2a_nsamFilled") == int(nsig) );
  
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
  return test_Roi2dToAdc(useExistingFcl);
}

//**********************************************************************
