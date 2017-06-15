// DuneToolManager.cxx

#include "DuneToolManager.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include <iostream>

using std::string;
using std::cout;
using std::endl;

//**********************************************************************

DuneToolManager& DuneToolManager::instance(string fclname) {
  static std::unique_ptr<DuneToolManager> pins(new DuneToolManager(fclname));
  return *pins.get();
}

//**********************************************************************

DuneToolManager::DuneToolManager(std::string fclname)
: m_fclname(fclname) {
  cet::filepath_lookup policy("FHICL_FILE_PATH");
  fhicl::intermediate_table tbl;
  fhicl::parse_document(fclname, policy, tbl);
  fhicl::ParameterSet psTop;
  fhicl::make_ParameterSet(tbl, psTop);
  m_pstools = psTop.get<fhicl::ParameterSet>("tools");
  m_toolNames = m_pstools.get_pset_names();
}

//**********************************************************************

const std::vector<std::string>& DuneToolManager::toolNames() const {
  return m_toolNames;
}

//**********************************************************************

void DuneToolManager::print() const {
  cout << "Tools for " << m_fclname << ":" << endl;
  for ( string name : toolNames() ) {
    cout << "  " << name << endl;
  }
}

//**********************************************************************
