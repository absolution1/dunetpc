// Roi2dToAdc_tool.cc

#include "Roi2dToAdc.h"

using std::string;
using std::cout;
using std::endl;
using std::vector;
using fhicl::ParameterSet;

using Index = Tpc2dRoi::Index;
using LongIndex = Tpc2dRoi::LongIndex;
using IndexArray = Tpc2dRoi::IndexArray;

//**********************************************************************
// Class methods.
//**********************************************************************

Roi2dToAdc::Roi2dToAdc(fhicl::ParameterSet const& ps)
: m_LogLevel(ps.get<int>("LogLevel")) {
  const string myname = "Roi2dToAdc::ctor: ";
  // Display the configuration.
  if ( m_LogLevel>= 1 ) {
    cout << myname << "             LogLevel: " << m_LogLevel << endl;
  }
}

//**********************************************************************

Roi2dToAdc::~Roi2dToAdc() { }

//**********************************************************************

DataMap Roi2dToAdc::updateTpcData(TpcData& tpd) const {
  const string myname = "Roi2dToAdc::updateTpcData: ";
  DataMap ret;
  const float zero = 0.0;
  Index nerr = 0;
  Index nchaZeroed = 0;
  Index nsamZeroed = 0;
  if ( m_LogLevel >= 3 ) cout << myname << "ADC map count " << tpd.getAdcData().size() << endl;
  // Zero the ADC samples, signal, ROIs and DFTs.
  for ( TpcData::AdcDataPtr pacm : tpd.getAdcData() ) {
    if ( ! pacm ) {
      if ( m_LogLevel >= 4 ) cout << myname << "  Skipping null map." << endl;
      continue;
    }
    if ( m_LogLevel >= 4 ) cout << myname << "  Processing map with channel count "
                                << pacm->size() << endl;
    for ( auto& iacd : *(pacm) ) {
      AdcChannelData& acd = iacd.second;
      Index nsam = acd.samples.size();
      if ( m_LogLevel >= 5 ) cout << myname << "    Zeroing channel " << acd.channel()
                                  << " with sample count " << nsam << endl;
      std::fill(acd.samples.begin(), acd.samples.end(), zero);
      nsamZeroed += nsam;
      acd.signal.resize(nsam);
      std::fill(acd.signal.begin(), acd.signal.end(), false);
      acd.rois.clear();
      acd.dftmags.clear();
      acd.dftphases.clear();
      ++nchaZeroed;
    }
  }
  // Copy ROI samples to ADC.
  if ( m_LogLevel >= 3 ) cout << myname << "ROI count " << tpd.get2dRois().size() << endl;
  Index nchaFilled = 0;
  Index nsamFilled = 0;
  Index iroi = 0;
  for ( const Tpc2dRoi& roi : tpd.get2dRois() ) {
    if ( m_LogLevel >= 4 ) cout << myname << "  Processing ROI " << iroi << " with "
                                << roi.channelSize() << " channels and "
                                << roi.sampleSize() << " samples." << endl;
    Index nsamRoi = roi.sampleSize();
    LongIndex ioffRoi = roi.sampleOffset();
    IndexArray idxs;
    Index& kcha = idxs[0];
    Index& ksamRoi = idxs[1];
    for ( kcha=0; kcha<roi.channelSize(); ++kcha ) {
      Index icha = kcha + roi.channelOffset();
      Index foundChannel = 0;
      for ( TpcData::AdcDataPtr pacm : tpd.getAdcData() ) {
        if ( pacm && pacm->count(icha) ) {
          ++foundChannel = true;
          AdcChannelData& acd = (*pacm)[icha];
          Index nsamAdc = acd.samples.size();
          if ( m_LogLevel >= 5 ) cout << myname << "    Copying channel " << acd.channel()
                                      << " from ROI sample count " << nsamRoi
                                      << " to ADC sample count " << nsamAdc << endl;
          LongIndex ioffAdc = acd.tickOffset();
          bool agr = ioffAdc > ioffRoi;
          Index isamAdc = agr ? 0 : Index(ioffRoi - ioffAdc);
          ksamRoi = agr ? Index(ioffAdc - ioffRoi) : 0;
          while ( ksamRoi < nsamRoi && isamAdc < nsamAdc ) {
            acd.samples[isamAdc] = roi.data().value(idxs);
            acd.signal[isamAdc] = true;
            ++ksamRoi;
            ++isamAdc;
            ++nsamFilled;
          }  // End loop over copied samples
        }
      }  // End loop over ADC channel maps
      if ( foundChannel ) ++nchaFilled;
    }  // end loop over channels in the ROI
    ++iroi;
  }  // end loop over ROIs
  // Add ROI info.
  for ( TpcData::AdcDataPtr pacm : tpd.getAdcData() ) {
    if ( ! pacm ) continue;
    for ( auto& iacd : *(pacm) ) {
      AdcChannelData& acd = iacd.second;
      if ( m_LogLevel >= 5 ) cout << myname << "    Recording ROIs for channel " << acd.channel() << endl;
      acd.roisFromSignal();
    }
  }
  // Fill result.
  ret.setStatus(nerr);
  ret.setInt("r2a_nchaZeroed", nchaZeroed);
  ret.setInt("r2a_nsamZeroed", nsamZeroed);
  ret.setInt("r2a_nchaFilled", nchaFilled);
  ret.setInt("r2a_nsamFilled", nsamFilled);
  return ret;
}

//**********************************************************************

DEFINE_ART_CLASS_TOOL(Roi2dToAdc)
