// liblist.cxx
//
// David Adams
// September 2016
//
// List all libraries accessible on the library path.

#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <memory>
#include "cetlib/LibraryManager.h"

using std::string;
using std::cout;
using std::endl;
using std::ostream;
using std::istringstream;
using cet::LibraryManager;

typedef std::vector<std::string> NameList;

int main(int argc, char** argv) {
  const string myname = "listlibs";
  bool help = false;
  string libtype = "service";
  string pattern;
  if ( argc > 1 ) {
    string arg = argv[1];
    if ( arg == "-h" ) help = true;
    else libtype = arg;
  }
  if ( argc > 2 ) {
    string arg = argv[2];
    if ( arg == "-h" ) help = true;
    else pattern = arg;
  }
  if ( help ) {
    cout << "Usage: " << argv[0] << " [type] [pattern]" << endl;
    cout << "  type [service] - library type (service, module, ...)" << endl;
    cout << "  pattern [] - pattern in library name" << endl;
    return 0;
  }
  // Find the libs.
  NameList libnames;
  LibraryManager libmgr(libtype);
  unsigned int nlib = libmgr.getLoadableLibraries(libnames);
  if ( pattern.size() ) {
    NameList newnames;
    for ( string libname : libnames ) {
      if ( libname.find(pattern) != string::npos ) newnames.push_back(libname);
    }
    libnames = newnames;
  }
  unsigned int nsel = libnames.size();
  if ( pattern.size() && nlib > 0 ) cout << "Selected " << nsel << " of ";
  else cout << "Found ";
  cout << nlib << " " << libtype << " libraries." << endl;
  for ( string libname : libnames ) {
    cout << "  " << libname << endl;
  }
  return 0;
}
