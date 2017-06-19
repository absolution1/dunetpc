// DuneToolManager.cxx

#include "DuneToolManager.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>

using std::string;
using std::cout;
using std::endl;
using std::ifstream;
using std::ostringstream;
using std::isspace;
using Index = unsigned int;
using NameVector = std::vector<string>;

//**********************************************************************

DuneToolManager* DuneToolManager::instance(string a_fclname) {
  const string myname = "DuneToolManager::instance: ";
  const int dbg = 0;
  if ( dbg > 0 ) cout << myname << "Called with " << a_fclname << endl;
  static string fclname;
  static std::unique_ptr<DuneToolManager> pins;
  if ( !pins ) {
    fclname = a_fclname;
    if ( fclname == "" ) {
      // Use ps to discover the command line.
      Index pid = getpid();
      ostringstream ssftmp;
      ssftmp << "tmpproc" << pid << ".tmp";
      string sftmp = ssftmp.str();
      ostringstream sscom;
      sscom << "ps -fp " << pid << " >" << sftmp;
      string scom = sscom.str();
      system(scom.c_str());
      ifstream fin(sftmp.c_str());
      string hdrline;
      getline(fin, hdrline);   // Skip header
      string::size_type iposCom = hdrline.find("CMD");
      string longline;
      getline(fin, longline);
      string line = longline.substr(iposCom);
      scom = "rm " + sftmp;
      system(scom.c_str());
      cout << myname << "Taking fcl name from command line: " << line << endl;
      // Split command line into words.
      NameVector words;
      bool newWord = true;
      for ( string::size_type ipos=0; ipos<line.size(); ++ipos ) {
        char ch = line[ipos];
        //if ( ch == '\0' || ch == 0 || isspace(ch) ) {
        if ( isspace(ch) ) {
          newWord = true;
        } else {
          if ( newWord ) {
            newWord = false;
            words.emplace_back();
          }
          words.back() += ch;
        }
      }
      // Find the flag "-c".
      Index iflg = words.size();
      for ( Index iwrd=0; iwrd<words.size(); ++iwrd ) {
        if ( words[iwrd] == "-c" ) iflg = iwrd;
      }
      Index ifcl = iflg + 1;
      if ( ifcl < words.size() ) fclname = words[ifcl];
      else {
        cout << myname << "ERROR: unable to retrieve configuration file name from command line." << endl;
        fclname = "tools_dune.fcl";
      }
    }
    if ( fclname.size() ) {
      cout << myname << "Configuration file: " << fclname << endl;
      pins.reset(new DuneToolManager(fclname));
    }
  } else {
    if ( a_fclname.size() && a_fclname != fclname ) {
      cout << myname << "WARNING: Ignoring inconsistent configuration name: " << a_fclname << endl;
    }
  }
  return pins.get();
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
