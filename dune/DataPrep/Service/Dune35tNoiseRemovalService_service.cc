// Dune35tNoiseRemovalService_service.cc

#include "Dune35tNoiseRemovalService.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
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
: m_LogLevel(1), m_nwarn(0) {
  const string myname = "Dune35tNoiseRemovalService::ctor: ";
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  m_GroupingFlag      = pset.get<int>("GroupingFlag");
  m_SkipStuckCodes    = pset.get<bool>("SkipStuckCodes");
  m_SkipSignals       = pset.get<bool>("SkipSignals");
  m_CorrectStuckCodes = pset.get<bool>("CorrectStuckCodes");
  m_ShowGroups        = pset.get<int>("ShowGroups");
  m_ShowGroupsChannel = pset.get<int>("ShowGroupsChannel");
  // Get services.
  m_pGeometry = &*ServiceHandle<geo::Geometry>();
  m_pChannelMap = &*ServiceHandle<lbne::ChannelMapService>();
  // Build the map of offline channel labels.
  const unsigned int nchan = m_pGeometry->Nchannels();
  for ( unsigned int chanon=0; chanon<nchan; ++chanon ) {
    unsigned int chanoff = m_pChannelMap->Offline(chanon);
    geo::View_t gview = m_pGeometry->View(chanoff);
    unsigned int iview = gview;
#if 0
    // Switch to this when larcore GeometryCore supports ROPs.
    // https://cdcvs.fnal.gov/redmine/issues/9264
    readout::ROPID ropid = m_pGeometry->ChannelToROP(chanoff);
    unsigned int rop = ropid.ROP;
    unsigned int chan0 = m_pGeometry->FirstChannelInROP(ropid);
    unsigned int ropchan = chanoff - chan0;
#else
    // Hard wire 35 ton convention.
    unsigned int ropchan = 999;
    {
      geo::WireID wid = m_pGeometry->ChannelToWire(chanoff).at(0);
      unsigned int tpcoff = wid.TPC;
      if ( gview == geo::kW && tpcoff%2 ) iview = 3;
      unsigned int apachan = chanoff%512;
      if ( iview == 0 ) ropchan = apachan;
      if ( iview == 1 ) ropchan = apachan - 144;
      if ( iview == 2 ) ropchan = apachan - 288;
      if ( iview == 3 ) ropchan = apachan - 400;
    }
#endif
    static string sview[4] = {"u", "v", "z1", "z2"};
    unsigned int apa = m_pChannelMap->APAFromOnlineChannel(chanon);
    ostringstream ssrop;
    ssrop << apa << sview[iview];
    string srop = ssrop.str();   // Name for the ROP
    ssrop << "-";
    if ( ropchan < 100 ) ssrop << "0";
    if ( ropchan <  10 ) ssrop << "0";
    ssrop << ropchan;
    string sropchan = ssrop.str();   // Name for the ROP channel
    m_sRopChannelMap[chanon] = sropchan;
  }
  // Build the channel groupings.
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
    if ( m_LogLevel >= 2 ) {
      cout << myname << "Online channel: plane:RCE:ASIC ROP-ropchan" << endl;
    }
    for ( unsigned int chanon=0; chanon<nchan; ++chanon ) { //online channels
      unsigned int plane     = m_pChannelMap->PlaneFromOnlineChannel(chanon);
      unsigned int rce       = m_pChannelMap->RCEFromOnlineChannel(chanon);
      unsigned int asic      = m_pChannelMap->ASICFromOnlineChannel(chanon);
      if ( m_LogLevel >= 2 ) {
        cout << myname << chanon << ": " << plane << ":" << rce << ":" << asic
             << "  " << m_sRopChannelMap.find(chanon)->second
             << endl;
      }
      m_GroupChannels[plane][rce*8+asic].push_back(chanon);
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
          // Skip stuck codes.
          // The original code also skipped samples where the raw count or pedestal was less than 10.
          if ( m_SkipStuckCodes && ( flag==AdcStuckOff || flag==AdcStuckOn ) ) continue;
          // Skip signals.
          bool isSignal = false;
          if ( m_SkipSignals ) {
            if ( acd.signal.size() > isig ) {
              isSignal = acd.signal[isig];
            } else if ( m_nwarn < 100 ) {
              ++m_nwarn;
              cout << myname << "WARNING: Signal skip requested without signal finding." << endl;
            }
          }
          if ( m_SkipSignals && isSignal ) continue;
          // Add current sample to correction.
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
  out << prefix << "      ShowGroupsChannel: " << m_ShowGroupsChannel     << endl;
  if ( m_ShowGroups ) {
    int wori = 10;
    int wgrp = 10;
    int wcha = 5;
    if ( m_ShowGroups == 1 ) {
      out << prefix << setw(wori) << "Orient." << setw(wgrp) << "Group" << ": " << "  Online channels" << endl;
    } else if ( m_ShowGroups == 2 ) {
      out << prefix << setw(wgrp+4) << "Group" << ": ";
      if      ( m_ShowGroupsChannel == 0 ) out << " Channels" << endl;
      else if ( m_ShowGroupsChannel == 1 ) out << " Online channels" << endl;
      else if ( m_ShowGroupsChannel <= 4 ) out << " Offline channels" << endl;
    }
    vector<vector<string>> gochans;   // Store channel lists as strings ordering group before orient
    unsigned int nori = m_GroupChannels.size();
    for ( unsigned int iori=0; iori<nori; ++iori ) {
      const std::vector<AdcChannelVector>& gchans = m_GroupChannels[iori];
      for ( unsigned int igrp=0; igrp<gchans.size(); ++igrp ) {
        const AdcChannelVector& chans = gchans[igrp];
        if ( chans.size() == 0 ) continue;
        std::vector<string> schans;
        bool usesetw = false;
        bool addspace = false;
        for ( AdcChannel chan : chans ) {
          if ( m_ShowGroupsChannel == 0 ) {
            schans.push_back(".");
          } else if ( m_ShowGroupsChannel == 1 ) {
            ostringstream sout;
            sout << chan;
            schans.push_back(sout.str());
            usesetw = true;
          } else if ( m_ShowGroupsChannel == 2 ) {
            ostringstream sout;
            sout << m_pChannelMap->Offline(chan);
            schans.push_back(sout.str());
            usesetw = true;
          } else if ( m_ShowGroupsChannel == 3 ) {
            schans.push_back(m_sRopChannelMap.find(chan)->second);
            addspace = true;
          } else if ( m_ShowGroupsChannel == 4 ) {
            string schan = m_sRopChannelMap.find(chan)->second;
            string::size_type pos = schan.find("z1-");
            if ( pos != string::npos ) schan.replace(pos, 2, "z");
            pos = schan.find("z2-");
            if ( pos != string::npos ) schan.replace(pos, 2, "Z");
            schans.push_back(schan);
            addspace = true;
          }
        }
        std::sort(schans.begin(), schans.end());
        ostringstream sschanlist;
        for ( string schan : schans ) {
          if ( usesetw ) sschanlist << setw(wcha);
          else if ( addspace && sschanlist.str().size() ) sschanlist << " ";
          sschanlist << schan;
        }
        string schanlist = sschanlist.str();
        if ( m_ShowGroups == 1 ) {
          out << prefix << setw(wori) << iori << setw(wgrp) << igrp << ": " << schanlist << endl;
        } else if ( m_ShowGroups == 2 ) {
          if ( gochans.size() < igrp+1 ) gochans.resize(igrp+1, std::vector<string>(nori));
          gochans[igrp][iori] = schanlist;
        }
      }
    }
    if ( m_ShowGroups == 2 ) {
      for ( unsigned int igrp=0; igrp<gochans.size(); ++igrp ) {
        for ( unsigned int iori=0; iori<nori; ++iori ) {
          string schanlist = gochans[igrp][iori];
          if ( schanlist.size() ) out << prefix << setw(wgrp+2) << igrp << "-" << iori << ": " << schanlist << endl;
        }
      }
    }
      
  }
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(Dune35tNoiseRemovalService, AdcNoiseRemovalService)

//**********************************************************************
