// ArtServiceHelper.cxx

#include "dune/ArtSupport/ArtServiceHelper.h"
#include <iostream>
#include <iomanip>
#include <memory>
#include "cetlib/filepath_maker.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceToken.h"
#include "art/Framework/EventProcessor/ServiceDirector.h"
#include "art/Framework/Services/Registry/ServiceRegistry.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "art/Framework/Services/System/CurrentModule.h"

using std::string;
using std::ostream;
using std::cout;
using std::endl;
using std::setw;
using std::find;
using std::unique_ptr;
using std::vector;
using cet::filepath_lookup;
using fhicl::ParameterSet;
using fhicl::make_ParameterSet;
using art::ServiceToken;
using art::ServiceDirector;
using art::ActivityRegistry;
using art::ServiceRegistry;
using art::TriggerNamesService;
using art::CurrentModule;

typedef ArtServiceHelper::NameList NameList;
typedef ArtServiceHelper::ConfigurationMap ConfigurationMap;

//**********************************************************************
// Local helpers.
//**********************************************************************

namespace {

// Return the art activity register.
ActivityRegistry& ar() {
  static unique_ptr<ActivityRegistry> par(new ActivityRegistry);
  return *par.get();
}

class ArtServiceHelperCloser {
public:
  ~ArtServiceHelperCloser() {
    cout << "ArtServiceHelperCloser::dtor: Closing ArtServiceHelper." << endl;
    ArtServiceHelper::close();
    cout << "ArtServiceHelperCloser::dtor: Closed." << endl;
  }
};

}  // end unnamed namespace

//**********************************************************************
// Public members.
//**********************************************************************

ArtServiceHelper& ArtServiceHelper::instance() {
  return *instancePtr().get();
}
  
//**********************************************************************

ArtServiceHelper& ArtServiceHelper::load(string fname) {
  const string myname = "ArtServiceHelper::instance: ";
  ArtServiceHelper& ins = *instancePtr().get();
  if ( ins.fileNames().size() == 0 ) {
    if ( fname.size() ) ins.loadServices(fname);
  } else {
    if ( ins.fileNames().size() > 1 ||
         ins.fileNames()[0] != fname ) {
      cout << myname << "ArtServiceHelper is already configured with different services." << endl;
      cout << myname << "The existing instance is returned." << endl;
    }
  }
  return ins;
}
  
//**********************************************************************

void ArtServiceHelper::close() {
  if ( instance().m_load == 3 ) return;
  if ( instance().m_load != 1 ) return;
  // Close existing services and registry.
  cout << "ArtServiceHelper::close: Closing art services." << endl;
  ArtServiceHelper& obj = instance();
  ServiceRegistry::Operate* pop = static_cast<ServiceRegistry::Operate*>(obj.m_poperate);
  delete pop;
  instance().m_poperate = nullptr;            // Reset the old operate.
  instancePtr().reset(new ArtServiceHelper);  // Delete the old service helper.
  instance().m_load = 3;                      // Put new instance in deleted state.
}
  
//**********************************************************************

ArtServiceHelper::ArtServiceHelper()
: m_LogLevel(0) { }

//**********************************************************************

ArtServiceHelper::~ArtServiceHelper() { }

//**********************************************************************

int ArtServiceHelper::addService(string name, string sval, bool isFile) {
  const string myname = "ArtServiceHelper::addService: ";
  if ( m_LogLevel > 1 ) {
    cout << myname << "Adding service " << name;
    if ( isFile ) cout << " from FCL file " << sval << "." << endl;
    else cout << " from FCL text." << endl;
  }
  if ( m_load == 1 || m_load == 2 ) {
    cout << myname << myname << "ERROR: Services may not be added after services are loaded." << endl;
    return 1;
  } else if ( m_load == 3 ) {
    cout << myname << myname << "ERROR: Services may not be added after service helper is closed." << endl;
    return 3;
  } else if ( m_load != 0 ) {
    cout << myname << myname << "ERROR: Unexpected load status: " << m_load << "." << endl;
    return 4;
  }
  if ( find(m_names.begin(), m_names.end(), name) != m_names.end() ) {
    cout << myname << myname << "ERROR: Service " << name << " is already registered" << endl;
    return 5;
  }
  if ( name == "TriggerNamesService" ) {
    cout << myname << myname << "ERROR: TriggerNamesService is configured automatically." << endl;
    return 6;
  }
  if ( name == "CurrentModule" ) {
    cout << myname << myname << "ERROR: CurrentModule service is configured automatically." << endl;
    return 7;
  }
  string scfg;
  if ( isFile ) {
    string fname = sval;
    string path = getenv("FHICL_FILE_PATH");
    if ( path.size() == 0 ) path = ".";
    filepath_lookup fpm(path); 
    string filepath;
    try { filepath = fpm(fname); }
    catch(...) {
      cout << myname << "ERROR: Unable to find file " << fname << endl;
      cout << myname << "Search path is " << path << endl;
      return 8;
    }
    if ( m_LogLevel > 1 ) cout << myname << "FCL file path: " << filepath << endl;
    // Fetch the cfg for full file.
    ParameterSet cfg_file;
    make_ParameterSet(fname, fpm, cfg_file);
    int fail = 0;
    bool fileHasName = false;
    bool servicesHasName = false;
    bool userHasName = false;
    const ParameterSet* pcfg = nullptr;
    ParameterSet cfgsrv;
    if ( cfg_file.is_empty() ) {
      fail = 1;
    } else {
      // Get configuration for services and determine
      // where the service configuration is present.
      ParameterSet cfg_services;
      ParameterSet cfg_user;
      fileHasName= cfg_services.has_key(name);
      if ( cfg_file.get_if_present<ParameterSet>("services", cfg_services) ) {
        servicesHasName = cfg_services.has_key(name);
        if ( cfg_services.get_if_present<ParameterSet>("user", cfg_user) ) {
          userHasName = cfg_user.has_key(name);
        }
      }
      if ( userHasName ) {
        pcfg = &cfg_user;
        if ( servicesHasName )
          if ( m_LogLevel ) cout << myname << "WARNING: Service " << name << " found in services blocks." << endl;
        if ( fileHasName )
          if ( m_LogLevel ) cout << myname << "WARNING: Service " << name << " found in services block and at file level." << endl;
      } else if ( servicesHasName ) {
        pcfg = &cfg_services;
        if ( fileHasName )
          if ( m_LogLevel ) cout << myname << "WARNING: Service " << name << " found in services block and at file level." << endl;
      } else if ( fileHasName ) {
        pcfg = &cfg_file;
      } else {
        fail = 2;
      }
      // Get the configuration for this service.
      if ( fail == 0 ) {
        if ( pcfg == nullptr ) {
          fail = 3;
        } else {
          if ( pcfg->get_if_present<ParameterSet>(name, cfgsrv) ) {
            scfg = cfgsrv.to_string();
          } else {
            fail = 4;
          }
        }
      }
    }
    if ( fail ) {
      cout << myname << "ERROR: Unable to find service " << name << " in file " << fname << endl;
      cout << myname << "ERROR: Resolved file path: " << fpm(fname) << endl;
      cout << myname << "ERROR: ";
      if ( fail == 1 ) {
        cout << "File is empty." << endl;
      } else if ( fail == 2 ) {
        cout << "Block not found in file, services: " << name << endl;
      } else if ( fail ) {
        cout << "Unexpected error " << fail << endl;
      }
      return 9;
    }
  } else {
    scfg = sval;
  }
  // scfgload embeds the the configuration in a named block and adds the field service_type
  // if it is not already present.
  ParameterSet cfg;
  make_ParameterSet(scfg, cfg);
  string scfgload = name + ": { " + scfg;
  if ( ! cfg.has_key("service_type") ) {
    scfgload += " service_type: " + name;
    scfgload += " }";
  } else {
    scfgload += " }";
    if ( m_LogLevel ) {
      cout << myname << "WARNING: Service " << name << " has predefined value for service_type." << endl;
      cout << myname << scfgload << endl;
    }
  }
  m_names.push_back(name);
  m_cfgmap[name] = scfg;
  if ( m_scfgs.size() ) m_scfgs += " ";
  m_scfgs += scfgload;
  if ( name == "TFileService" ) m_needTriggerNamesService = true;
  if ( name == "RandomNumberGenerator" ) m_needCurrentModuleService = true;
  return 0;
}

//**********************************************************************

int ArtServiceHelper::addServices(string sval, bool isFile) {
  const string myname = "ArtServiceHelper::addServices: ";
  if ( m_LogLevel > 1 ) {
    cout << myname << "Adding services";
    if ( isFile ) cout << " from FCL file " << sval << "." << endl;
    else cout << " from FCL text." << endl;
  }
  if ( m_load == 1 || m_load == 2 ) {
    cout << myname << myname << "ERROR: Services may not be added after services are loaded." << endl;
    return 1;
  } else if ( m_load == 3 ) {
    cout << myname << myname << "ERROR: Services may not be added after service helper is closed." << endl;
    return 3;
  } else if ( m_load != 0 ) {
    cout << myname << myname << "ERROR: Unexpected load status: " << m_load << "." << endl;
    return 4;
  }
  // Make vector of parameter sets.
  vector<ParameterSet> psets;
  psets.reserve(2);
  if ( isFile ) {
    string fname = sval;
    string path = getenv("FHICL_FILE_PATH");
    if ( path.size() == 0 ) path = ".";
    filepath_lookup fpm(path); 
    string filepath;
    try { filepath = fpm(fname); }
    catch(...) {
      cout << myname << "ERROR: Unable to find file " << fname << endl;
      cout << myname << "Search path is " << path << endl;
      return 8;
    }
    if ( m_LogLevel > 1 ) cout << myname << "FCL file path: " << filepath << endl;
    // Fetch the cfg for full file.
    ParameterSet cfg_file;
    make_ParameterSet(fname, fpm, cfg_file);
    if ( cfg_file.is_empty() ) {
      cout << myname << "ERROR: Empty configuration file: " << fname << endl;
      cout << myname << "ERROR: Resolved file path: " << fpm(fname) << endl;
      return 5;
    } else {
      // Get configuration for services and determine
      // where the service configuration is present.
      psets.push_back(ParameterSet());
      const ParameterSet& pset_services = psets.back();
      if ( cfg_file.get_if_present<ParameterSet>("services", psets.back())) {
        psets.push_back(ParameterSet());
        pset_services.get_if_present<ParameterSet>("user", psets.back());
      }
    }
    if ( psets.size() == 0 ) {
      cout << myname << "ERROR: Unable to find service block in file " << fname << endl;
      cout << myname << "ERROR: Resolved file path: " << fpm(fname) << endl;
      return 9;
    }
    m_fnames.push_back(fname);
  } else {
    psets.push_back(ParameterSet());
    make_ParameterSet(sval, psets.back());
    m_fnames.push_back("STRING");
  }
  // Loop over the service parameter sets and add each service description to m_scfgs.
  for ( const ParameterSet& pset : psets ) {
    std::vector<std::string> snames = pset.get_pset_names();
    for ( string sname : snames ) {
      if ( sname == "user" ) continue;
      if ( sname == "message" ) continue;
      if ( std::find(m_names.begin(), m_names.end(), sname) != m_names.end() ) {
        cout << myname << "WARNING: Ignoring duplicate service entry: " << sname << endl;
        continue;
      }
      cout << myname << "Adding service " << sname << endl;
      ParameterSet cfg = pset.get<ParameterSet>(sname);
      string scfg = cfg.to_string();
      string scfgload = sname + ": { " + scfg;
      if ( ! cfg.has_key("service_type") ) {
        scfgload += " service_type: " + sname;
      }
      scfgload += " }";
      ParameterSet cfgload;
      make_ParameterSet(scfgload, cfgload);
      m_names.push_back(sname);
      m_cfgmap[sname] = scfg;
      if ( m_scfgs.size() ) m_scfgs += " ";
      m_scfgs += scfgload;
      if ( sname == "TFileService" ) m_needTriggerNamesService = true;
      if ( sname == "RandomNumberGenerator" ) m_needCurrentModuleService = true;
    }
  }
  return 0;
}

//**********************************************************************

int ArtServiceHelper::loadServices() {
  string myname = "ArtServiceHelper::loadServices: ";
  if ( m_load == 1 ) {
    cout << myname << "ERROR: Services are already loaded." << endl;
    return 1;
  }
  if ( m_load == 2 ) {
    cout << myname << "ERROR: Service load failed earlier." << endl;
    return 2;
  }
  if ( m_load == 3 ) {
    cout << myname << "ERROR: Service manager is closed." << endl;
    return 2;
  }
  // Register standard services.
  ParameterSet cfgServices;
  make_ParameterSet(fullServiceConfiguration(), cfgServices);
  ServiceToken serviceToken = ServiceToken::createInvalid();
  ServiceDirector director(std::move(cfgServices), ar(), serviceToken);
  // Register special services.
  if ( m_needTriggerNamesService ) {
    ParameterSet cfgTriggerNamesService;
    string scfgTriggerNamesService = "process_name: \"myproc\"";
    make_ParameterSet(scfgTriggerNamesService, cfgTriggerNamesService);
    std::vector<string> tns;
    director.addSystemService<TriggerNamesService>(cfgTriggerNamesService, tns);
  }
  if ( m_needCurrentModuleService ) {
    director.addSystemService<CurrentModule>(ar());
  }
  // Make the services available
  // We need to keep this object around if we want to use the services.
  //// m_poperate.reset(new ServiceRegistry::Operate(serviceToken));
  m_poperate = new ServiceRegistry::Operate(serviceToken);
  return m_load = 1;
}

//**********************************************************************

int ArtServiceHelper::loadServices(string fclfile) {
  string myname = "ArtServiceHelper::loadServices: ";
  if ( m_fnames.size() ) {
    cout << myname << "ERROR: Services have already been added." << endl;
    return 100;
  }
  int astat = addServices(fclfile, true);
  if ( astat != 0 ) return 100 + astat;
  return loadServices();
}

//**********************************************************************

void ArtServiceHelper::setLogLevel(int lev) {
  const string myname = "ArtServiceHelper::setLogLevel: ";
  if ( m_LogLevel > 1 || lev > 1 ) cout << myname << "Setting log level to " << lev << endl;
  m_LogLevel = lev;
}

//**********************************************************************

NameList ArtServiceHelper::serviceNames() const {
  return m_names;
}

//**********************************************************************

NameList ArtServiceHelper::fileNames() const {
  return m_fnames;
}

//**********************************************************************

string ArtServiceHelper::serviceConfiguration(std::string name) const {
  ConfigurationMap::const_iterator icfg = m_cfgmap.find(name);
  if ( icfg == m_cfgmap.end() ) return "";
  return icfg->second;
}

//**********************************************************************

string ArtServiceHelper::fullServiceConfiguration() const {
  return m_scfgs;
}

//**********************************************************************

int ArtServiceHelper::serviceStatus() const {
  return m_load;
}

//**********************************************************************

bool ArtServiceHelper::isReady() const {
  return m_load == 1;
}

//**********************************************************************

void ArtServiceHelper::print(ostream& out) const {
  if ( m_load == 3 ) {
    out << "Service manager is closed" << endl;
    return;
  }
  out << "ArtServiceHelper has " << m_names.size() << " "
      << (m_names.size()==1 ? "entry" : "entries")
      << (m_names.size()==0 ? "" : ":") << endl;
  unsigned int wnam = 8;
  for ( string name : m_names ) if ( name.size() > wnam ) wnam = name.size();
  for ( string name : m_names ) {
    out << setw(wnam+2) << name << ": { " << m_cfgmap.find(name)->second << " }" << endl;
  }
  if ( m_load == 0 ) out << "Services have not been loaded." << endl;
  else if ( m_load == 1 ) out << "Services have been loaded and are available." << endl;
  else if ( m_load == 2 ) out << "Service load failed." << endl;
  else out << "Unrecognized load status: " << m_load << "." << endl;
}

//**********************************************************************
// Private members.
//**********************************************************************

std::unique_ptr<ArtServiceHelper>& ArtServiceHelper::instancePtr() {
  static unique_ptr<ArtServiceHelper> psh(new ArtServiceHelper);
  // DLA Jan 2017: I did a test with Root and this close comes too late.
  // See https://cdcvs.fnal.gov/redmine/issues/15162#note-9.
  // static ArtServiceHelperCloser closer;
  return psh;
}
  
//**********************************************************************
