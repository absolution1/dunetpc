// fcldump.cxx
//
// David Adams
// September 2015

#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <memory>
#include "cetlib/filepath_maker.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

using std::string;
using std::cout;
using std::endl;
using std::ostream;
using std::istringstream;
using std::make_shared;
using cet::filepath_lookup;
using fhicl::ParameterSet;
using fhicl::make_ParameterSet;

typedef std::shared_ptr<ParameterSet> PSPtr;
typedef std::vector<PSPtr> PSPtrVector;
typedef std::vector<int> IntVector;
typedef std::vector<double> DoubleVector;
typedef std::vector<string> StringVector;
typedef std::vector<DoubleVector> DVVector;
typedef std::vector<DVVector> DVVVector;
typedef std::vector<ParameterSet> PSVector;
typedef std::vector<PSPtr> PSPVector;

//**********************************************************************

// Delimiter.

template<class T>
string delim() { return ""; }

template<>
string delim<string>() { return "\""; }

//**********************************************************************

// Print a sequence.

template<class T>
ostream& operator<<(ostream& out, const std::vector<T>& vec) {
  out << "[";
  bool first = true;
  for ( auto val : vec ) {
    if ( ! first ) cout << ", ";
    out << delim<T>() << val << delim<T>();
    first = false;
  }
  out << "]";
  return out;
}

//**********************************************************************

// Print a block.

void print_block(string prefix, PSPtr pcfgs[], unsigned int nlevrem) {
  const string myname = "fcldump:print_block: ";
  ParameterSet* pcfg = pcfgs[0].get();
  // First print the values.
  for ( string key : pcfg->get_names() ) {
    if ( pcfg->is_key_to_table(key) ) continue;
    cout << prefix << key << ": ";
    if ( pcfg->is_key_to_sequence(key) ) {
      try {
        cout << pcfg->get<IntVector>(key);
      } catch (...) {
        try {
          cout << pcfg->get<DoubleVector>(key);
        } catch (...) {
          try {
            cout << pcfg->get<DVVector>(key);
          } catch (...) {
            try {
              cout << pcfg->get<DVVVector>(key);
            } catch (...) {
              try {
                cout << pcfg->get<StringVector>(key);
              } catch (...) {
                try {
                  // Trick case: sequnce of paremter sets.
                  PSVector psets = pcfg->get<PSVector>(key);
                  cout << "[";
                  int submaxlev = nlevrem;
                  if ( submaxlev > 0 ) {
                    for ( const ParameterSet& pset : psets ) {
                      cout << endl;
                      PSPtrVector subcfgs(submaxlev, std::make_shared<ParameterSet>(pset));
                      subcfgs[0] = std::make_shared<ParameterSet>(pset);
                      cout << prefix << "{" << endl;
                      print_block(prefix + "  ", &subcfgs[0], submaxlev-1);
                      cout << prefix << "}";
                    }
                    cout << endl << prefix << "]";
                  } else {
                    for ( unsigned int ips=0; ips<psets.size(); ++ips ) cout << ".";
                    cout << "]";
                  }
                } catch (fhicl::exception& exc) {
                  cout << "ERROR: Unknown data type for sequence key " << key << endl;
                  throw exc;
                }
              }
            }
          }
        }
      }
    } else {
      try {
        int val = pcfg->get<int>(key);
        cout << val;
      } catch (...) {
        try {
          double val = pcfg->get<double>(key);
          cout << val;
        } catch (...) {
          try {
            cout << "\"" << pcfg->get<string>(key) << "\"";
          } catch (fhicl::exception& exc) {
            cout << "ERROR: Unknown data type for non-sequence key " << key << endl;
            throw exc;
          }
        }
      }
    }
    cout << endl;
  }
  // Next the blocks.
  for ( string key : pcfg->get_pset_names() ) {
    cout << prefix << key << ": {";
    PSPtr pcfgnext(new ParameterSet);
    if ( ! pcfg->get_if_present<ParameterSet>(key, *pcfgnext) ) {
      cout << " ERROR!!!!!!!!!!!!!!!!" << endl;
      break;
    }
    int nKeysNext = pcfgnext->get_names().size();
    if ( ! pcfgnext->get_names().size() ) {
      cout << " }";
    } else if ( nlevrem > 0 ) {
      std::swap(pcfgs[1], pcfgnext);
      cout << endl;
      print_block(prefix + "  ", &pcfgs[1], nlevrem-1);
      cout << prefix << "}";
    } else {
      for ( int ikey=0; ikey<nKeysNext; ++ikey ) cout << ".";
      cout << "}";
    }
    cout << endl;
  }
}
  
int main(int argc, char** argv) {
  const string myname = "fcldump: ";
  bool help = argc == 1;
  unsigned int maxlev = 0;
  string fname;
  if ( argc > 1 ) {
    string arg = argv[1];
    if ( arg == "-h" ) help = true;
    else fname = argv[1];
  }
  if ( argc > 2 ) {
    istringstream ssarg(argv[2]);
    ssarg >> maxlev;
  }
  if ( help ) {
    cout << "Usage: " << argv[0] << " FCLFILE [DEPTH]" << endl;
    cout << "  DEPTH > 0 dumps fcl contents to that depth." << endl;
    cout << "  Full path to the fcl file is displayed if DEPTH is omitted" << endl;
    return 0;
  }

  // Find the file.
  string path = getenv("FHICL_FILE_PATH");
  if ( path.size() == 0 ) path = ".";
  filepath_lookup fpm(path);
  string filepath;
  try { filepath = fpm(fname); }
  catch(...) {
    cout << myname << "ERROR: Unable to find file " << fname << endl;
    cout << myname << "Search path:" << endl;
    string::size_type ipos = 0;
    while ( ipos != string::npos ) {
      string::size_type jpos = path.find(":", ipos+1);
      cout << "  " << path.substr(ipos, jpos-ipos) << endl;
      if ( jpos == string::npos ) break;
      ipos = jpos + 1;
    }
    return 8;
  }

  cout << filepath << endl;
  if ( maxlev <= 0 ) return 0;
  cout << endl;

  // Fetch top-level configuration and print.
  PSPtrVector cfgs(maxlev, std::make_shared<ParameterSet>());
  make_ParameterSet(fname, fpm, *cfgs[0]);
  string prefix;
  print_block(prefix, &cfgs[0], maxlev-1);
  return 0;
}
