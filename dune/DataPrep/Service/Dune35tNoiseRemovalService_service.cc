// Dune35tNoiseRemovalService_service.cc

#include "Dune35tNoiseRemovalService.h"
#include <iostream>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "lbne-raw-data/Services/ChannelMap/ChannelMapService.h"

using std::vector;
using std::string;
using std::ostream;
using std::cout;
using std::endl;
using art::ServiceHandle;

//**********************************************************************

Dune35tNoiseRemovalService::
Dune35tNoiseRemovalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1) {
  const string myname = "Dune35tNoiseRemovalService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  // Get services.
  m_pGeometry = &*ServiceHandle<geo:Geometry>;
  m_pChannelMap &*ServiceHandle<lbne::ChannelMapService>;
  print(cout, myname);
}

//**********************************************************************

int Dune35tNoiseRemovalService::update(AdcChannelDataMap& datamap) const {
  const string myname = "Dune35tNoiseRemovalService:update: ";
  if ( datamap.size() == 0 ) {
    cout << myname << "WARNING: No channels found." << endl;
    return 1;
  }
  unsigned int nsig = 0;
  bool first = true;
  for ( const AdcChannelDataMap::value_type& ent : datamap ) {
    const AdcChannelData& data = ent.second;
    if ( first ) nsig = data.samples.size();
    if ( data.samples.size() != nsig ) {
      cout << myname << "WARNING: Channels have inconsistent sample counts." << endl;
      return 2;
    }
  }
  if ( nsig == 0 ) {
    cout << myname << "WARNING: No ADC samples found." << endl;
    return 3;
  }
  if ( m_LogLevel >= 3 ) {
    cout << myname << "Entering..." << endl;
    cout << myname << "  # Channels: " << datamap.size() << endl;
    cout << myname << "     # Ticks: " << nsig << endl;
  }
  const unsigned int nchan = fGeom->Nchannels();
  // Build the channel groupings.
  vector<vector<AdcChannelVector>>> grpchans;  //grpchans[plane id][group id][channel id]
  if ( m_GroupingFlag == 0 ) { //group by regulators
    grpchans.resize(3);
    for (size_t i = 0; i<3; ++i){
      grpchans[i].resize(32);
    }
    for (unsigned int i = 0; i<nchan; ++i){//online channels
      unsigned int plane     = m_pChannelMap->PlaneFromOnlineChannel(i);
      unsigned int rce       = m_pChannelMap->RCEFromOnlineChannel(i);
      unsigned int regulator = m_pChannelMap->RegulatorFromOnlineChannel(i);
      grpchans[plane][rce*2+regulator].push_back(i);
    }
  } else if ( m_GroupingFlag == 1 ) {
    grpchans.resize(3);
    for (size_t i = 0; i<3; ++i){
      grpchans[i].resize(128);
    }
    for (unsigned int i = 0; i<nchan; ++i){//online channels
      unsigned int plane     = m_pChannelMap->PlaneFromOnlineChannel(i);
      unsigned int rce       = m_pChannelMap->RCEFromOnlineChannel(i);
      unsigned int asic      = m_pChannelMap->ASICFromOnlineChannel(i);
      //cout<<i<<" "<<plane<<" "<<rce<<" "<<asic<<endl;
      grpchans[plane][rce*8+asic].push_back(i);
    }
  } else {
    cout << myname << "ERROR: Invalid GroupingFlag: " << m_GroupingFlag << endl;
    return 4;
  }
  // derive correction factors - require raw adc waveform and pedestal for each channel
  vector<double> corrVals;
  // loop through time slices
  for ( unsigned int isig=0; isig<nsig; ++isig) {
    // Loop over wire orientations.
    for ( size_t iori=0; iori<grpchans.size(); ++iori ){
      // Loop over channel groups.
      for ( size_t igrp=0; igrp<grpchans[i].size(); ++igrp ){
        const AdcChannelVector& chans = grpchans[iori][igrp];
        if ( m_logLevel >= 3 && isig==0 ) {
          cout << myname << setw(3) << iori << setw(4) << igrp << " " << chans.size() << endl;
        }
        grpcorrVals.clear();
        // Loop over the channels for this group
        for ( size_t icha : chans ) {
          unsigned int offlineChan = m_pChannelMap->Offline(pchans[icha]);
          if ( rawDigitMap.count(offlineChan) == 0 ) continue;
          AdcSignal sig = datamap[offlineChan].samples[isig];
          AdcFlag flag = datamap[offlineChan].samples[isig];
          if ( m_SkipStuckCodes && ( flag==AdcStuckOff || flag==AdcStuckOn ) ) continue;
          // The original code also skipped samples where the raw count or pedestal was less than 10.
          corrVals.push_back(sig);
        }
        // Evaluate the correction as the median of the selected signals.
        unsigned int corrValSize = corrVals.size();
        sort(corrVals.begin(),corrVals.end());
        double correction = 0;
        if ( corrValSize < 2 )
          correction = 0.0;
        else if ( (corrValSize % 2) == 0 )
          correction = (corrVals[corrValSize/2] + corrVals[(corrValSize/2)-1])/2.0;
        else
          correction = corrVals[(corrValSize-1)/2];
        // Loop over the channels for this group
        for ( size_t icha : chans ) {
          unsigned int offlineChan = m_pChannelMap->Offline(chans[icha]);
          if ( datamap.count(offlineChan) == 0 ) continue;
          if ( !m_CorrectStuckCodes && ( flag==AdcStuckOff || flag==AdcStuckOn ) ) continue;
          AdcSignal& sig = datamap[offlineChan].samples[isig];
          sig -= correction;
        }  // end loop over channels in the group
      }  // end loop over groups
    }  // end loop over orientations
  }  // end loop over time samples
  return 0;
}

//**********************************************************************

ostream& Dune35tNoiseRemovalService::
print(ostream& out, string prefix) const {
  out << prefix << "Dune35tNoiseRemovalService:"                   << endl;
  out << prefix << "               LogLevel: " << m_LogLevel              << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(Dune35tNoiseRemovalService, AdcMitigationService)

//**********************************************************************
