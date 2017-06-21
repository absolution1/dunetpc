// DuneToolManager.h
//
// David Adams
// April 2017
//
// This singleton provides access to DUNE tools. It is a quick and dirty
// implementation of a tool manager intended for use until we get something
// better from larsoft or art.

#ifndef DuneToolManager_H
#define DuneToolManager_H

#include <string>
#include <vector>
#include <memory>
#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/make_tool.h"

class DuneToolManager {

public:

  using Name = std::string;
  using NameVector = std::vector<Name>;

  // Return the one instance of this class.
  // The name is ignored once the tool manager is set.
  // if fclname is blank and the tool manager is not yet set, an attempt
  // is made to find the fcl name on the command line following "-c"..
  static DuneToolManager* instance(std::string fclname ="");

  // Return a private (not shared) copy of a tool.
  template<class T>
  std::unique_ptr<T> getPrivate(std::string name) {
    fhicl::ParameterSet psTool = m_pstools.get<fhicl::ParameterSet>(name);
    return art::make_tool<T>(psTool);
  }

  // Return the list of available tool names.
  const std::vector<std::string>& toolNames() const;

  // Display the list of available tools.
  void print() const;

private:

  // Ctor from FCL file name.
  explicit DuneToolManager(std::string fclname);

  std::string m_fclname;
  fhicl::ParameterSet m_pstools;
  NameVector m_toolNames;

};

#endif
