// DuneToolManager.h
//
// David Adams
// April 2017
//
// This class provides access to DUNE tools. It is a quick and dirty
// implementation of a tool manager intended for use until we get something
// better from larsoft or art.
//
// Tools are typically obtained from the primary instance of this class which
// is accessed via DuneToolManager::instance(). The first call to this method
// locates the top-level fcl file and loads all the tool definitions. Tools
// are constructed when getPrivate or getPublic is called.
//
// This first call typically occurs during job initialization. The fcl file
// name may be provided or, if not, it is obtained from the command line
// (option -c). If that is also absent, the default tools_dune.fcl is used.
//
// It is also possible to construct additional tool managers from fcl file names.
// Tool instances are not shared between different tool managers. These
// managers might be assigned to different threads or used to locate tools
// that are not defined in the primary fcl file.

#ifndef DuneToolManager_H
#define DuneToolManager_H

#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/make_tool.h"

class DuneToolManager {

public:  // Type aliases and subclasses.

  using Name = std::string;
  using NameVector = std::vector<Name>;

  class SharedToolEntry {
  public:
    virtual ~SharedToolEntry() =default;
  };

  template<class T>
  class TSharedToolEntry : public SharedToolEntry {
  public:
    using ToolPtr = std::unique_ptr<T>;
    TSharedToolEntry(ToolPtr& a_ptr) : m_ptr(std::move(a_ptr)) { }
    ~TSharedToolEntry() =default;
    T* get() { return m_ptr.get(); }
  private:
    ToolPtr m_ptr;
  };

  using SharedToolPtr = std::unique_ptr<SharedToolEntry>;
  using SharedToolMap = std::map<Name, SharedToolPtr>;

public:

  // Return the primary instance of this class.
  // The name is ignored once the tool manager is set.
  // if fclname is blank and the tool manager is not yet set, an attempt
  // is made to find the fcl name on the command line following "-c"..
  static DuneToolManager* instance(std::string fclname ="");

  // Ctor from FCL file name.
  explicit DuneToolManager(std::string fclname);

  // Return a private (not shared) copy of a tool.
  template<class T>
  std::unique_ptr<T> getPrivate(std::string name) {
    std::string myname = "DuneToolManager::getPrivate: ";
    fhicl::ParameterSet psTool;
    if ( name.size() == 0 ) {
      std::cout << myname << "Tool name is blank" << std::endl;
      return nullptr;
    } else if ( name[0] == '{' ) {    // Name is a tool cfg string
      makeParameterSet(name, psTool);
    } else if ( std::find(m_toolNames.begin(), m_toolNames.end(), name) != m_toolNames.end() ) {
      psTool = m_pstools.get<fhicl::ParameterSet>(name);
    } else {
      std::cout << myname << "No such tool name: " << name << std::endl;
      return nullptr;
    }
    return art::make_tool<T>(psTool);
  }

  // Return a shared tool.
  template<class T>
  T* getShared(std::string name) {
    SharedToolMap::iterator itoo = m_sharedTools.find(name);
    if ( itoo != m_sharedTools.end() ) {
      SharedToolPtr& pent = itoo->second;
      TSharedToolEntry<T>* ptent = dynamic_cast<TSharedToolEntry<T>*>(pent.get());
      return ptent->get();
    }
    fhicl::ParameterSet psTool = m_pstools.get<fhicl::ParameterSet>(name);
    std::unique_ptr<T> ptoo = art::make_tool<T>(psTool);
    TSharedToolEntry<T>* ptent(new TSharedToolEntry<T>(ptoo));
    SharedToolEntry* pent = ptent;
    //m_sharedTools[name] = pent;
    m_sharedTools.emplace(name, pent);
    return ptent->get();
  }

  // Return the list of available tool names.
  const std::vector<std::string>& toolNames() const;

  // Display the list of available tools.
  void print() const;

private:

  std::string m_fclname;
  fhicl::ParameterSet m_pstools;
  NameVector m_toolNames;
  SharedToolMap m_sharedTools;

  // Convert a string into a paramter set.
  int makeParameterSet(std::string scfg, fhicl::ParameterSet& ps);
  
};

#endif
