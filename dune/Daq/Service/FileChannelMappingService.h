// FileChannelMappingService

// David Adams
// February 2016
//
// Implementation of online-offline channel mapping reading from a file.
//
// Parameters:
//   FileName - Name of the channel map file. Search done in FW_SEARCH_PATH.
//   FilePathEnv ["FW_SEARCH_PATH"] - Name of the env variable nolding the file search path.
//   LogLevel [1] - Message logging level: 0=none, 1=init, 2+=every call to interface
//
// The parameters FilePathEnv and LogLevel are optional and have the indicated default values.

#ifndef FileChannelMappingService_H
#define FileChannelMappingService_H

#include "dune/DuneInterface/ChannelMappingService.h"
#include <vector>
#include <iostream>

class TH1;
namespace CLHEP {
class HepRandomEngine;
}

class FileChannelMappingService : public ChannelMappingService {

  typedef std::vector<Channel> ChannelMap;

public:

  // Ctor.
  FileChannelMappingService(fhicl::ParameterSet const& pset);

  // Ctor.
  FileChannelMappingService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Map online to offline.
  Channel offline(Channel onlineChannel) const;

  // Map offline to online.
  Channel online(Channel offlineChannel) const;

  // Print the parameters.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:
 
  // Parameters.
  std::string m_FileName;
  std::string m_FilePathEnv;
  int m_LogLevel;

  // Maps.
  ChannelMap m_onMap;
  ChannelMap m_offMap;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(FileChannelMappingService, ChannelMappingService, LEGACY)

#endif
