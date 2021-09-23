// ArtServiceHelper.h

#ifndef ArtServiceHelper_H
#define ArtServiceHelper_H

// =================================================================
// The purpose of the 'ArtServiceHelper' class is to construct and
// manage a set of art services that normally appear in an art job.
// In some circumstances, the available means of testing a new
// algorithm may rely on functionality that has been expressed in
// terms of an art service.  In such a case, the ArtServiceHelper
// class can be used to initialize the required services and allow
// users to create art::ServiceHandle objects for the configured
// services.
//
// Unlike the art framework's services manager, the ArtServiceHelper
// is not aware of framework transitions.  Because of that, using some
// services outside of the context of the framework will not reflect
// the same behavior as using them within it.  For that reason, using
// such services outside of the framework is circumspect.  This must
// be taken into account when using this class.
//
// Initialization
// ==============
//
// The 'ArtServiceHelper' is initialized by specifying the same kind
// of configuration one might specify within art (see below under the
// Configuration layout section).  Allowed initialization patterns
// include:
//
// (1) string-based configuration:
//
//       std::string const config{"MyService: {...}"};
//       ArtServiceHelper::load_services(config);
//
//     The string is allowed to '#include' any files that will be
//     accessible via the FHICL_FILE_PATH environment variable.
//
// (2) filename-based configuration:
//
//       ArtServiceHelper::load_services("config.fcl",
//                                       ArtServiceHelper::FileOnPath);
//
//     where 'config.fcl' is a filename relative to one of the
//     directories on the FHICL_FILE_PATH environment variable.
//
// (3) stream-based configuration:
//
//     It is permissible to specify an input stream that will be used
//     to parse the configuration.  This can be helpful when a
//     configuration file is required that is not accessible via
//     FHICL_FILE_PATH:
//
//       std::ifstream in{"/some/path/to_some_file.fcl"};
//       ArtServiceHelper::load_services(in);
//
//     Note however, that FHICL_FILE_PATH lookup is *still* enabled in
//     this mode--i.e. configurations within the file can '#include'
//     files accessible via FHICL_FILE_PATH.  Stream-based
//     configuration can also be used to assemble more complicated
//     configurations via an std::stringstream object:
//
//       std::stringstream config;
//       config << "#include \"some_config.fcl\"\n"
//              << "PedestalCalibration.some_value: " << 15 << '\n;
//       ArtServiceHelper::load_services(config);
//
// (4) ParameterSet-based configuration:
//
//     If a fhicl::ParameterSet object is already available, it can be
//     used to initialize the ArtServiceHelper:
//
//       auto const pset = get_parameter_set(...);
//       ArtServiceHelper::load_services(pset);
//
// Configuration layout
// ====================
//
// The allowed configuration layout is of two forms.  The minimal form
// lists all desired services by themselves:
//
//   Service1: {...}
//   Service2: {...}
//   ServiceN: {...}
//
// The other form nests these services inside of a 'services' table,
// which supports cases when it is desired to use a configuration file
// that would be used for an art job:
//
//   services: {
//     Service1: {...}
//     Service2: {...}
//     ServiceN: {...}
//   }
//
// As in art, it is not necessary to specify the 'service_type'
// parameter for each service--the ArtServiceHelper class inserts
// those configuration parameters automatically.
//
// =================================================================

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServicesManager.h"
#include "art/Utilities/SharedResource.h"

#include <iosfwd>
#include <map>
#include <string>

class ArtServiceHelper {
public:
  struct FileOnPath_t {};
  constexpr static FileOnPath_t FileOnPath{};

  static void load_services(std::string const& config);
  static void load_services(std::string const& filename, FileOnPath_t);
  static void load_services(std::istream& config);
  static void load_services(fhicl::ParameterSet const& pset);

  // For backward compatibility.
  static void load(std::string const& filename) { load_services(filename, FileOnPath_t{}); }

private:
  explicit ArtServiceHelper(fhicl::ParameterSet&& pset);
  art::ActivityRegistry activityRegistry_;
  art::detail::SharedResources sharedResources_;
  art::ServicesManager servicesManager_;
};

#endif
