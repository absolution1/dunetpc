// Dune35tNoiseRemovalService_service.cc

#include "Dune35tNoiseRemovalService.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "lbne-raw-data/Services/ChannelMap/ChannelMapService.h"

using std::vector;
using std::string;
using std::ostream;
using std::cout;
using std::endl;
using std::ostringstream;
using std::setw;
using art::ServiceHandle;

//**********************************************************************

Dune35tNoiseRemovalService::
Dune35tNoiseRemovalService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_LogLevel(1) {
  const string myname = "Dune35tNoiseRemovalService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  m_GroupingFlag      = pset.get<int>("GroupingFlag");
  m_SkipStuckCodes    = pset.get<bool>("SkipStuckCodes");
  m_CorrectStuckCodes = pset.get<bool>("CorrectStuckCodes");
  m_ShowGroups        = pset.get<int>("ShowGroups");
  // Get services.
  m_pGeometry = &*ServiceHandle<geo::Geometry>();
  m_pChannelMap = &*ServiceHandle<lbne::ChannelMapService>();
  // Build the channel groupings.
  const unsigned int nchan = m_pGeometry->Nchannels();
  if ( m_GroupingFlag == 0 ) { //group by regulators
    m_GroupChannels.resize(3);
    for (size_t i = 0; i<3; ++i){
      m_GroupChannels[i].resize(32);
    }
    for (unsigned int i=0; i<nchan; ++i){//online channels
      unsigned int plane     = m_pChannelMap->PlaneFromOnlineChannel(i);
      unsigned int rce       = m_pChannelMap->RCEFromOnlineChannel(i);
      unsigned int regulator = m_pChannelMap->RegulatorFromOnlineChannel(i);
      m_GroupChannels[plane][rce*2+regulator].push_back(i);
    }
  } else if ( m_GroupingFlag == 1 ) {
    m_GroupChannels.resize(3);
    for (size_t i = 0; i<3; ++i){
      m_GroupChannels[i].resize(128);
    }
    for (unsigned int i=0; i<nchan; ++i){//online channels
      unsigned int plane     = m_pChannelMap->PlaneFromOnlineChannel(i);
      unsigned int rce       = m_pChannelMap->RCEFromOnlineChannel(i);
      unsigned int asic      = m_pChannelMap->ASICFromOnlineChannel(i);
      //cout<<i<<" "<<plane<<" "<<rce<<" "<<asic<<endl;
      m_GroupChannels[plane][rce*8+asic].push_back(i);
    }
  } else {
    cout << myname << "ERROR: Invalid GroupingFlag: " << m_GroupingFlag << endl;
    return;
  }
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
  if ( m_GroupChannels.size() == 0 ) {
    cout << myname << "ERROR: Channel groupings have not been defined." << endl;
    return 4;
  }
  if ( m_LogLevel >= 3 ) {
    cout << myname << "Entering..." << endl;
    cout << myname << "  # Channels: " << datamap.size() << endl;
    cout << myname << "     # Ticks: " << nsig << endl;
  }
  // loop through time slices
  for ( unsigned int isig=0; isig<nsig; ++isig) {
    // Loop over wire orientations.
    for ( size_t iori=0; iori<m_GroupChannels.size(); ++iori ){
      // Loop over channel groups.
      for ( size_t igrp=0; igrp<m_GroupChannels[iori].size(); ++igrp ){
        const AdcChannelVector& chans = m_GroupChannels[iori][igrp];
        if ( m_LogLevel >= 3 && isig==0 ) {
          cout << myname << setw(3) << iori << setw(4) << igrp << " " << chans.size() << endl;
        }
        // Loop over the channels for this time slice, orientation and group
        vector<double> corrVals;
        for ( AdcChannel chan : chans ) {
          unsigned int offlineChan = m_pChannelMap->Offline(chan);
          if ( m_LogLevel >= 4 && isig==0 ) {
            cout << myname << "Channel on/off: " << chan << "/" << offlineChan << endl;
          }
          AdcChannelDataMap::const_iterator iacd = datamap.find(offlineChan);
          if ( iacd == datamap.end() ) continue;
          const AdcChannelData& acd = iacd->second;
          AdcSignal sig = acd.samples[isig];
          AdcFlag flag = acd.flags.size() ? acd.flags[isig] : AdcGood;
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
        for ( AdcChannel chan : chans ) {
          unsigned int offlineChan = m_pChannelMap->Offline(chan);
          AdcChannelDataMap::iterator iacd = datamap.find(offlineChan);
          if ( iacd == datamap.end() ) continue;
          AdcChannelData& acd = iacd->second;
          AdcSignal& sig = acd.samples[isig];
          AdcFlag flag = acd.flags.size() ? acd.flags[isig] : AdcGood;
          if ( m_SkipStuckCodes && ( flag==AdcStuckOff || flag==AdcStuckOn ) ) continue;
          if ( datamap.count(offlineChan) == 0 ) continue;
          if ( !m_CorrectStuckCodes && ( flag==AdcStuckOff || flag==AdcStuckOn ) ) continue;
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
  out << prefix << "Dune35tNoiseRemovalService:"                          << endl;
  out << prefix << "               LogLevel: " << m_LogLevel              << endl;
  out << prefix << "           GroupingFlag: " << m_GroupingFlag          << endl;
  out << prefix << "         SkipStuckCodes: " << m_SkipStuckCodes        << endl;
  out << prefix << "      CorrectStuckCodes: " << m_CorrectStuckCodes     << endl;
  out << prefix << "             ShowGroups: " << m_ShowGroups            << endl;
  if ( m_ShowGroups ) {
    int wori = 10;
    int wgrp = 10;
    int wcha = 5;
    if ( m_ShowGroups == 1 ) {
      out << prefix << setw(wori) << "Orient." << setw(wgrp) << "Group" << ": " << "  Online channels" << endl;
    } else if ( m_ShowGroups == 2 ) {
      out << prefix << setw(wgrp+4) << "Group" << ": " << " Online channels" << endl;
    }
    vector<vector<string>> gochans;   // Store channel lists as strings ordering group before orient
    unsigned int nori = m_GroupChannels.size();
    for ( unsigned int iori=0; iori<nori; ++iori ) {
      const std::vector<AdcChannelVector>& gchans = m_GroupChannels[iori];
      for ( unsigned int igrp=0; igrp<gchans.size(); ++igrp ) {
        const AdcChannelVector& chans = gchans[igrp];
        if ( chans.size() == 0 ) continue;
        ostringstream sout;
        for ( AdcChannel chan : chans ) sout << setw(wcha) << chan;
        if ( m_ShowGroups == 1 ) {
          out << prefix << setw(wori) << iori << setw(wgrp) << igrp << ": " << sout.str() << endl;
        } else if ( m_ShowGroups == 2 ) {
          if ( gochans.size() < igrp+1 ) gochans.resize(igrp+1, std::vector<string>(nori));
          gochans[igrp][iori] = sout.str();
        }
      }
    }
    if ( m_ShowGroups == 2 ) {
      for ( unsigned int igrp=0; igrp<gochans.size(); ++igrp ) {
        for ( unsigned int iori=0; iori<nori; ++iori ) {
          string schans = gochans[igrp][iori];
          if ( schans.size() ) out << prefix << setw(wgrp+2) << igrp << "-" << iori << ": " << schans << endl;
        }
      }
    }
      
  }
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(Dune35tNoiseRemovalService, AdcNoiseRemovalService)

//**********************************************************************
