// AdcToRoi2d_tool.cc

#include "AdcToRoi2d.h"

using std::string;
using std::cout;
using std::endl;
using std::vector;
using fhicl::ParameterSet;

using Index = AdcToRoi2d::Index;
using LongIndex = unsigned long int;

namespace {

string optionString(Index iopt) {
  if ( iopt == 0 ) return "No copy";
  if ( iopt == 1 ) return "One ROI/map";
  if ( iopt == 2 ) return "One ROI for all maps";
  return "Invalid option";
}

}  // end unnamed namespace

//**********************************************************************
// Class methods.
//**********************************************************************

AdcToRoi2d::AdcToRoi2d(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")),
  m_Option(ps.get<Index>("Option")),
  m_InputAdcMaps(ps.get<IndexVector>("InputAdcMaps")),
  m_OutputNames(ps.get<NameVector>("OutputNames")) {
  const string myname = "AdcToRoi2d::ctor: ";
  if ( m_Option > 2 ) {
    cout << myname << "WARNING: Invalid option: " << m_Option
         << " will be treated as 0." << endl;
  }
  // Display the configuration.
  if ( m_LogLevel>= 1 ) {
    cout << myname << "     LogLevel: " << m_LogLevel << endl;
    cout << myname << "       Option: " << m_Option
         << " (" << optionString(m_Option) << ")" << endl;
    cout << myname << "  OutputNames: [";
    bool first = true;
    for ( Name nam : m_OutputNames ) {
      if ( first ) first = false;
      else cout << ", ";
      cout << nam;
    }
    cout << "]" << endl;
  }
}

//**********************************************************************

AdcToRoi2d::~AdcToRoi2d() { }

//**********************************************************************

DataMap AdcToRoi2d::updateTpcData(TpcData& tpd) const {
  const string myname = "AdcToRoi2d::updateTpcData: ";
  DataMap ret;
  Index imax = -1;
  LongIndex lmax = -1;
  // Find the channel and tick ranges for each ADC map.
  Index nacm = tpd.getAdcData().size();
  vector<LongIndex> itck1s(nacm, lmax);
  vector<LongIndex> itck2s(nacm, 0);
  vector<Index> icha1s(nacm, imax);
  vector<Index> icha2s(nacm, 0);
  vector<Index> ndats(nacm, 0);
  Index nerr = 0;
  if ( m_LogLevel >= 3 ) cout << myname << "Channel map count: " << nacm << endl;
  IndexVector mapIndices = m_InputAdcMaps;
  if ( mapIndices.size() == 0 ) {
    for ( Index iacm=0; iacm<nacm; ++iacm ) mapIndices.push_back(iacm);
  }
  for ( Index iacm : mapIndices ) {
    if ( iacm >= tpd.getAdcData().size() ) {
      cout << myname << "WARNING: Skipping invalid ADC map index " << iacm << endl;
      continue;
    }
    LongIndex ndat = 0;
    const TpcData::AdcDataPtr pacm = tpd.getAdcData().at(iacm);
    if ( ! pacm ) {
      if ( m_LogLevel >=2 ) cout << myname << "Skipping missing ADC channel map." << endl;
      ++nerr;
      continue;
    }
    if ( pacm->size() == 0 ) {
      if ( m_LogLevel >=2 ) cout << myname << "Skipping empty ADC channel map." << endl;
      ++nerr;
      continue;
    }
    if ( m_LogLevel >= 3 ) cout << myname << "  Channel count: " << pacm->size() << endl;
    for ( const auto& iacd : *pacm ) {
      const AdcChannelData& acd = iacd.second;
      Index nsam = acd.samples.size();
      Index icha1 = acd.channel();
      Index icha2 = icha1 + 1;
      Index itck1 = acd.tickOffset();
      Index itck2 = itck1 + nsam;
      if ( icha1 < icha1s[iacm] ) icha1s[iacm] = icha1;
      if ( icha2 > icha2s[iacm] ) icha2s[iacm] = icha2;
      if ( itck1 < itck1s[iacm] ) itck1s[iacm] = itck1;
      if ( itck2 > itck2s[iacm] ) itck2s[iacm] = itck2;
      if ( m_LogLevel >= 3 ) cout << myname << "    Channel " << acd.channel()
                                  << " sample count: " << nsam << endl;
      ndat += nsam;
    }
    if ( icha1s[iacm] >= icha2s[iacm] ) {
      cout << myname << "ERROR: Channels out of range. Skipping ADC map " << iacm << "." << endl;
      ++nerr;
      continue;
    }
    if ( ndat == 0 ) {
      cout << myname << "ERROR: No samples found. Skipping ADC map " << iacm << "." << endl;
      ++nerr;
      continue;
    }
    if ( itck1s[iacm] >= itck2s[iacm] ) {
      cout << myname << "ERROR: Ticks out of range. Skipping ADC map " << iacm << "." << endl;
      ++nerr;
      continue;
    }
    ndats[iacm] = ndat;
    LongIndex ncha = icha2s[iacm] - icha1s[iacm];
    LongIndex ntck = itck2s[iacm] - itck1s[iacm];
    LongIndex ncht = ncha*ntck;
    double frac = double(ndat)/double(ncht);
    if ( m_LogLevel >= 2 ) {
      cout << myname << "  ADC map is " << ncha << " x " << ntck
           << " with fill fraction " << frac << endl;
      cout << myname << "  Channel range: [ " << icha1s[iacm] << ", " << icha2s[iacm] << ")" << endl;
      cout << myname << "  Tick range: [ " << itck1s[iacm] << ", " << itck2s[iacm] << ")" << endl;
    }
  }
  Index nroi = 0;
  vector<int> roiNdats;
  if ( m_Option == 1 ) {
    for ( Index iacm=0; iacm<nacm; ++iacm ) {
      Index icha1 = icha1s[iacm];
      Index icha2 = icha2s[iacm];
      LongIndex itck1 = itck1s[iacm];
      LongIndex itck2 = itck2s[iacm];
      if ( icha1 >= icha2 ) continue;
      if ( itck1 >= itck2 ) continue;
      Index ncha = icha2 - icha1;
      Index ntck = Index(itck2 - itck1);
      // Create and fill ROI.
      tpd.get2dRois().emplace_back(ncha, ntck, icha1, itck1);
      ++nroi;
      TpcData* ptpdo = &tpd;
      if ( m_OutputNames.size() ) {
        if ( iacm >= m_OutputNames.size() ) {
          cout << myname << "WARNING: No outpu name supplied for ADC map " << iacm << endl;
        } else {
          Name nam = m_OutputNames[iacm];
          ptpdo = tpd.getTpcData(nam);
          if ( ptpdo == nullptr ) {
            ptpdo = tpd.addTpcData(nam);
            if ( ptpdo == nullptr ) {
              cout << myname << "ERROR: Unable to add TpcData directory " << nam << endl;
              continue;
            } else if ( m_LogLevel >= 2 ) {
              cout << myname << "Added TpcData directory " << nam << endl;
            }
          } else if ( m_LogLevel >= 2 ) {
            cout << myname << "Using existing TpcData directory " << nam << endl;
          }
        }
      }
      Tpc2dRoi& roi = ptpdo->get2dRois().back();
      const TpcData::AdcDataPtr pacm = tpd.getAdcData().at(iacm);
      Tpc2dRoi::DataArray::IndexArray idxs;
      Index ndat = 0;
      for ( const auto& iacd : *pacm ) {
        const AdcChannelData& acd = iacd.second;
        idxs[0] = acd.channel() - icha1;
        idxs[1] = acd.tickOffset() - itck1;
        for ( Index isam=0; isam<acd.samples.size(); ++isam, ++idxs[1], ++ndat ) {
          roi.data().setValue(idxs, acd.samples[isam]);
        }
      }
      roiNdats.push_back(ndat);
      if ( ndat != ndats[iacm] ) {
        cout << myname << "ERROR: Unexpected fill count: " << ndat << " != " << ndats[iacm] << endl;
        ++nerr;
      }
    }
  } else if ( m_Option == 2 ) {
    cout << myname << "ERROR: Merging of ADC data is not yet supported." << endl;
  }
  ret.setStatus(nerr);
  ret.setInt("a2r_nroi", nroi);
  ret.setIntVector("a2r_nsams", roiNdats);
  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(AdcToRoi2d)
