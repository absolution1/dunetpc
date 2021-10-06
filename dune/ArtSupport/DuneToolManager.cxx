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

string DuneToolManager::fclFilename(string a_fclname, int dbg) {
  const string myname = "DuneToolManager::fclFilename: ";
  if ( dbg >= 2 ) cout << myname << "Called with " << a_fclname << endl;
  static string fclname;
  bool haveName = fclname.size();
  bool setName = a_fclname.size();
  if ( !haveName && !setName ) {
    // Use ps to discover the command line.
    Index pid = getpid();
    ostringstream ssftmp;
    ssftmp << "tmpproc" << pid << ".tmp";
    string sftmp = ssftmp.str();
    ostringstream sscom;
    sscom << "ps -fwwp " << pid << " >" << sftmp;
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
    if ( dbg >= 1 ) cout << myname << "Taking fcl name from command line: " << line << endl;
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
    // Find the flag "-c" or "--config" and use following text as the fcl name.
    for ( Index iwrd=0; iwrd<words.size(); ++iwrd ) {
      string word = words[iwrd];
      if ( word == "-c" || word == "--config") {
        if ( words.size() > iwrd+1 ) fclname = words[iwrd+1];
      } else if ( word.substr(0,2) == "-c" ) {
        fclname = word.substr(2);
      }
    }
    // If name was not found, switch to a default.
    if ( fclname.empty() ) {
      cout << myname << "ERROR: unable to retrieve configuration file name from command line." << endl;
      fclname = "tools_dune.fcl";
    }
  } else if ( !haveName && setName ) {
    fclname = a_fclname;
  } else if ( haveName && setName ) {
    if ( a_fclname != fclname ) {
      cout << myname << "WARNING: Ignoring inconsistent configuration name: "
           << a_fclname << " != " << fclname << endl;
    }
  }
  if ( fclname.empty() ) {
    cout << myname << "ERROR: Unexpected empty file name!" << endl;
    fclname = "nosuchfile.fcl";
  }
  return fclname;
}

//**********************************************************************

DuneToolManager* DuneToolManager::instance(string a_fclname, int dbg) {
  const string myname = "DuneToolManager::instance: ";
  if ( dbg >= 2 ) cout << myname << "Called with " << a_fclname << endl;
  static std::unique_ptr<DuneToolManager> pins;
  string fclname = fclFilename(a_fclname, dbg);
  if ( !pins ) pins.reset(new DuneToolManager(fclname));
  return pins.get();
}

//**********************************************************************

void DuneToolManager::help() {
  cout << "DuneToolManager provides access to named tools." << endl;
  cout << "It is typically accessed a singleton, e.g." << endl;
  cout << "  DuneToolManager* ptm = DuneToolManger::instance();" << endl;
  cout << "The name of the fcl file defining the tools may be provided as an argument." << endl;
  cout << endl;
  cout << "Available tools may be listed with:" << endl;
  cout << "  ptm->print()      // Displays the names of available tools." << endl;
  cout << "  ptm->toolNames()  // Returns a vector of the tool names." << endl;
  cout << endl;
  cout << "A private instance of a tool with type MyTool may be otained with:" << endl;
  cout << "  ptm->getPrivate<MyTool>(\"tnam\")" << endl;
  cout << "The tool configuration is taken from the fcl block tools.tnam." << endl;
  cout << "A unique pointer is returned so the tool instance is destroyed with that pointer." << endl;
  cout << endl;
  cout << "A shared instance of a tool may be obtained with:" << endl;
  cout << "  ptm->getShared<MyTool>(\"tnam\")" << endl;
  cout << "A bare pointer is returned and subsequent callers receive the same pointer." << endl;
  cout << "The tool manager deletes all shared tools when it is destroyed." << endl;
  cout << endl;
  cout << "In either case, tool instances are created only when they are obtained." << endl;
}

//**********************************************************************

DuneToolManager::DuneToolManager(std::string fclname)
: m_fclname(fclname) {
  cet::filepath_lookup_nonabsolute policy("FHICL_FILE_PATH");
  m_pstools =
    fhicl::ParameterSet::make(fhicl::parse_document(fclname, policy)).
    get<fhicl::ParameterSet>("tools");
  m_toolNames = m_pstools.get_pset_names();
}

//**********************************************************************

int DuneToolManager::deleteShared(std::string tnam) {
  SharedToolMap::iterator itoo = m_sharedTools.find(tnam);
  if ( itoo == m_sharedTools.end() ) return 1;
  m_sharedTools.erase(tnam);
  return 0;
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

int DuneToolManager::makeParameterSet(std::string scfgin, fhicl::ParameterSet& ps) {
  // Strip surrounding braces.
  if ( scfgin.size() < 1 ) return 1;
  string scfg;
  if ( scfgin[0] == '{' ) {
    if ( scfgin[scfgin.size()-1] != '}' ) return 2;
    scfg = scfgin.substr(1, scfgin.size()-2);
  } else {
    scfg = scfgin;
  }
  // Parse the configuration string.
  ps = fhicl::ParameterSet::make(scfg);
  return 0;
}

//**********************************************************************
