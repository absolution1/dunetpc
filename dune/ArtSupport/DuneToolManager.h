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

  // Return the one instance of this class.
  // The name is ignored once the tool manager is set.
  // if fclname is blank and the tool manager is not yet set, an attempt
  // is made to find the fcl name on the command line following "-c"..
  static DuneToolManager* instance(std::string fclname ="");

  // Return a private (not shared) copy of a tool.
  template<class T>
  std::unique_ptr<T> getPrivate(std::string name) {
    std::string myname = "DuneToolManager::getPrivate: ";
    if ( std::find(m_toolNames.begin(), m_toolNames.end(), name) == m_toolNames.end() ) {
      std::cout << myname << "No such tool name: " << name << std::endl;
      return nullptr;
    }
    fhicl::ParameterSet psTool = m_pstools.get<fhicl::ParameterSet>(name);
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

  // Ctor from FCL file name.
  explicit DuneToolManager(std::string fclname);

  std::string m_fclname;
  fhicl::ParameterSet m_pstools;
  NameVector m_toolNames;
  SharedToolMap m_sharedTools;

};

#endif
