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
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "art/Framework/Services/System/CurrentModule.h"

using std::string;
using std::ostream;
using std::cout;
using std::endl;
using std::setw;
using std::find;
using std::unique_ptr;
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

}  // end unnamed namespace

//**********************************************************************
// Public members.
//**********************************************************************

ArtServiceHelper& ArtServiceHelper::instance() {
  return *instancePtr().get();
}
  
//**********************************************************************

void ArtServiceHelper::close() {
  if ( instance().m_load == 3 ) return;
  delete instance().m_poperate;               // Close existing services and registry.
  instance().m_poperate = nullptr;            // Reset the old operatre.
  instancePtr().reset(new ArtServiceHelper);  // Delete the old service helper.
  instance().m_load = 3;                      // Put new instance in deleted state.
}
  
//**********************************************************************

int ArtServiceHelper::addService(string name, string sval, bool isFile) {
  string myname = "ArtServiceHelper::addService: ";
  if ( m_load == 1 || m_load == 2 ) {
    cout << myname << "ERROR: Services may not be added after services are loaded." << endl;
    return 1;
  } else if ( m_load == 3 ) {
    cout << myname << "ERROR: Services may not be added after service hlper is closed." << endl;
    return 3;
  } else if ( m_load != 0 ) {
    cout << myname << "ERROR: Unexpected load status: " << m_load << "." << endl;
    return 4;
  }
  if ( find(m_names.begin(), m_names.end(), name) != m_names.end() ) {
    cout << myname << "ERROR: Service " << name << " is already registered" << endl;
    return 5;
  }
  if ( name == "TriggerNamesService" ) {
    cout << myname << "ERROR: TriggerNamesService is configured automatically." << endl;
    return 6;
  }
  if ( name == "CurrentModule" ) {
    cout << myname << "ERROR: CurrentModule service is configured automatically." << endl;
    return 7;
  }
  string scfg;
  if ( isFile ) {
    string fname = sval;
    string path = getenv("FHICL_FILE_PATH");
    if ( path.size() == 0 ) path = ".";
    filepath_lookup fpm(path); 
    try { string filepath = fpm(fname); }
    catch(...) {
      cout << myname << "ERROR: Unable to find file " << fname << endl;
      cout << myname << "Search path is " << path << endl;
      return 8;
    }
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
      // Get configuration for services and services.user and determine
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
          cout << myname << "WARNING: Service " << name << " found in services.user and services blocks." << endl;
        if ( fileHasName )
          cout << myname << "WARNING: Service " << name << " found in services.user block and at file level." << endl;
      } else if ( servicesHasName ) {
        pcfg = &cfg_services;
        if ( fileHasName )
          cout << myname << "WARNING: Service " << name << " found in services block and at file level." << endl;
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
        cout << "Block not found in file, services or services.user: " << name << endl;
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
    cout << myname << "WARNING: Service " << name << " has predefined value for service_type." << endl;
    cout << myname << scfgload << endl;
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
  ServiceToken serviceToken;
  ServiceDirector director(cfgServices, ar(), serviceToken);
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

NameList ArtServiceHelper::serviceNames() const {
  return m_names;
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
  return psh;
}
  
//**********************************************************************
