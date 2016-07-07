// FileChannelMappingService.cxx

#include "dune/Daq/Service/FileChannelMappingService.h"
#include <limits>
#include <fstream>

using std::cout;
using std::ostream;
using std::endl;
using std::string;
using std::ifstream;

#undef UseSeedService

typedef ChannelMappingService::Channel Channel;

namespace {

// Bad channel value.
Channel bad() {
  Channel val = std::numeric_limits<Channel>::max();
  return val;
}

}

//**********************************************************************

FileChannelMappingService::
FileChannelMappingService(fhicl::ParameterSet const& pset)
: m_FilePathEnv("FW_SEARCH_PATH"), m_LogLevel(1) {
  const string myname = "FileChannelMappingService::ctor: ";
  m_FileName = pset.get<string>("FileName");
  pset.get_if_present<string>("FilePathEnv", m_FilePathEnv);
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  if ( m_LogLevel > 0 ) {
    print(cout, myname);
    cout << myname << "Reading map file." << endl;
  }
  // Read the map.
  cet::search_path sp(m_FilePathEnv);
  string fullname;
  sp.find_file(m_FileName, fullname);
  if ( fullname.empty() ) {
    cout << "ERROR: " << "Input file not found: " << m_FileName << endl;
    cout << "ERROR: " << "Search path: $" << m_FilePathEnv << endl;
    throw cet::exception("File not found");
  }
  std::ifstream infile(fullname);
  Channel chon;
  Channel choff;
  unsigned int count = 0;
  while ( infile.good() ) {
    infile >> chon >> choff;
    if ( m_onMap.size() < choff+1 ) m_onMap.resize(choff+1, bad());
    if ( m_offMap.size() < chon+1 ) m_offMap.resize(chon+1, bad());
    m_onMap[choff] = chon;
    m_offMap[chon] = choff;
    ++count;
  }
  if ( m_LogLevel > 0 ) {
    cout << myname << "Channel map file: " << fullname << endl;
    cout << myname << "    Number of map entries: " << count << endl;
    cout << myname << "   Maximum online channel: " << m_offMap.size()-1 << endl;
    cout << myname << "  Maximum offline channel: " << m_onMap.size()-1 << endl;
  }
}

//**********************************************************************

FileChannelMappingService::
FileChannelMappingService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: FileChannelMappingService(pset) { }

//**********************************************************************

Channel FileChannelMappingService::offline(Channel chin) const {
  const string myname  = "FileChannelMappingService::offline: ";
  Channel chout = bad();
  if ( chin < m_offMap.size() ) chout = m_offMap[chin];
  if ( chout == bad() ) {
    if ( m_LogLevel > 1 ) {
      cout << myname << "ERROR: " << "No offline channel mapped to online channel " << chin << endl;
    }
    throw cet::exception("Online channel not found.");
  }
  return chout;
}

//**********************************************************************

Channel FileChannelMappingService::online(Channel chin) const {
  const string myname  = "FileChannelMappingService::online: ";
  Channel chout = bad();
  if ( chin < m_onMap.size() ) chout = m_onMap[chin];
  if ( chout == bad() ) {
    if ( m_LogLevel > 1 ) {
      cout << myname << "ERROR: " << "No online channel mapped to offline channel " << chin << endl;
    }
    throw cet::exception("Offline channel not found.");
  }
  return chout;
}

//**********************************************************************

ostream& FileChannelMappingService::print(ostream& out, string prefix) const {
  string myname = prefix + "  ";
  cout << myname << "   FileName: " << m_FileName << endl;
  cout << myname << "FilePathEnv: " << m_FilePathEnv << endl;
  cout << myname << "   LogLevel: " << m_LogLevel << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(FileChannelMappingService, ChannelMappingService)

//**********************************************************************
