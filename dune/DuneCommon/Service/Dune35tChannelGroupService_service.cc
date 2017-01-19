// Dune35tChannelGroupService.cxx

#include "dune/DuneCommon/Service/Dune35tChannelGroupService.h"
#include <sstream>
#include <iomanip>

using std::cout;
using std::ostream;
using std::endl;
using std::string;
using std::ostringstream;
using std::setfill;
using std::setw;
using std::vector;

typedef ChannelGroupService::Index Index;
typedef ChannelGroupService::Name Name;
typedef ChannelGroupService::ChannelVector ChannelVector;
typedef vector<ChannelVector> ChannelVectorVector;
typedef vector<ChannelVectorVector> ChannelVectorVectorVector;

//**********************************************************************

Dune35tChannelGroupService::
Dune35tChannelGroupService(fhicl::ParameterSet const& pset)
: m_UseOffline(true), m_LogLevel(1) {
  const string myname = "Dune35tChannelGroupService::ctor: ";
  const Index nAsicPerRce = 8;
  const Index nRegPerRce = 2;
  const string plname[4] = {"u", "v", "z1", "z2"};
  // Fetch parameters.
  m_UseOffline = pset.get<bool>("UseOffline");
  m_SplitByPlane = pset.get<bool>("SplitByPlane");
  m_Grouping = pset.get<string>("Grouping");
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  const Index ncha = 2048;
  bool useAPA = m_Grouping == "APA";
  bool useASIC = m_Grouping == "ASIC";
  bool useRegulator = m_Grouping == "Regulator";
  bool split = m_SplitByPlane;
  bool off = m_UseOffline;
  Index ngrp = 0;
  ChannelVectorVectorVector vvvchans;
  Index npla = split ? 4 : 1;
  for ( Index icha=0; icha<ncha; ++icha ) {
    Index igrp = 9999;
    Index chon  = off ? m_pcms->Online(icha) : icha;
    Index irce = m_pcms->RCEFromOnlineChannel(chon);
    Index ipla = m_pcms->PlaneFromOnlineChannel(chon);
    if ( useAPA ) {
      igrp = m_pcms->APAFromOnlineChannel(chon);
    } else if ( useASIC ) {
      Index iasc = m_pcms->ASICFromOnlineChannel(chon);
      igrp = nAsicPerRce*irce + iasc;
    } else if ( useRegulator ) {
      Index ireg = m_pcms->RegulatorFromOnlineChannel(chon);
      igrp = nRegPerRce*irce + ireg;
    } else {
      cout << myname << "ERROR: Invalid group name: " << m_Grouping << endl;
      return;
    }
    if ( igrp >= ngrp ) {
      ngrp = igrp + 1;
      vvvchans.resize(ngrp, ChannelVectorVector(npla));
    }
    Index iview = split ? ipla : 0;
    vvvchans[igrp][iview].push_back(icha);
  }
  Name base = useAPA       ?  "apa" :
              useASIC      ? "asic" :
              useRegulator ?  "reg" : "Unknown";
  //if ( useAPA ) base = "apa";
  //if ( useASIC ) base = "asic";
  //if ( useRegulator ) base = "reg";
  for ( Index igrp=0; igrp<ngrp; ++igrp ) {
    for ( Index ipla=0; ipla<npla; ++ipla ) {
      ChannelVector& chans = vvvchans[igrp][ipla];
      if ( chans.size() ) {
        ostringstream ssnam;
        ssnam << base;
        if ( ngrp > 100 && igrp < 100 ) ssnam << "0";
        if ( ngrp >  10 && igrp <  10 ) ssnam << "0";
        ssnam << igrp;
        if ( split ) ssnam << plname[ipla];
        m_names.push_back(ssnam.str());
        m_chanvecs.push_back(chans);
      }
    }
  }
  if ( m_LogLevel > 0 ) print();
}

//**********************************************************************

Dune35tChannelGroupService::
Dune35tChannelGroupService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: Dune35tChannelGroupService(pset) { }

//**********************************************************************

Index Dune35tChannelGroupService::size() const {
  return m_chanvecs.size();
}

//**********************************************************************

Name Dune35tChannelGroupService::name(Index igrp) const {
  if ( igrp < m_names.size() ) return m_names[igrp];
  return "NoSuchGroup";
}

//**********************************************************************

const ChannelVector& Dune35tChannelGroupService::channels(Index igrp) const {
  if ( igrp < m_chanvecs.size() ) return m_chanvecs[igrp];
  static ChannelVector empty;
  return empty;
}

//**********************************************************************

ostream& Dune35tChannelGroupService::print(ostream& out, string prefix) const {
  out << prefix << "Dune35tChannelGroupService: " << endl;
  out << prefix << "    UseOffline: " << ( m_UseOffline ? "true" : "false") << endl;
  out << prefix << "      Grouping: " << m_Grouping << endl;
  out << prefix << "  SplitByPlane: " << m_SplitByPlane << endl;
  out << prefix << "      # groups: " << size() << endl;
  out << prefix << "      LogLevel: " << m_LogLevel << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(Dune35tChannelGroupService, ChannelGroupService)

//**********************************************************************
