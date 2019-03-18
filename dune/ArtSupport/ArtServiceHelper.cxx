#include "art/Framework/Services/Registry/ServiceRegistry.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/intermediate_table.h"

#include "dune/ArtSupport/ArtServiceHelper.h"

// The following is an ugly hack, but it is necessary for being able
// to set the services manager.
namespace art {
  class EventProcessor {
  public:
    static void set_services_manager(ServicesManager* manager)
    {
      ServiceRegistry::instance().setManager(manager);
    }
  };
}

ArtServiceHelper::ArtServiceHelper(fhicl::ParameterSet&& pset) :
  activityRegistry_{},
  servicesManager_{std::move(pset), activityRegistry_}
{
  art::EventProcessor::set_services_manager(&servicesManager_);
  servicesManager_.forceCreation();
}

void ArtServiceHelper::load_services(std::string const& config)
{
  std::istringstream iss{config};
  load_services(iss);
}

void ArtServiceHelper::load_services(std::istream& config)
{
  cet::filepath_lookup lookup{"FHICL_FILE_PATH"};
  fhicl::intermediate_table table;
  fhicl::parse_document(config, lookup, table);
  fhicl::ParameterSet pset;
  fhicl::make_ParameterSet(table, pset);
  load_services(std::move(pset));
}

void ArtServiceHelper::load_services(std::string const& filename, FileOnPath_t)
{
  cet::filepath_lookup lookup{"FHICL_FILE_PATH"};
  fhicl::intermediate_table table;
  fhicl::parse_document(filename, lookup, table);
  fhicl::ParameterSet pset;
  fhicl::make_ParameterSet(table, pset);
  load_services(std::move(pset));
}

void ArtServiceHelper::load_services(fhicl::ParameterSet const& pset)
{
  // Make sure each service has a value for the "service_type" parameter
  fhicl::ParameterSet post_processed_pset;
  auto const& input = pset.has_key("services") ? pset.get<fhicl::ParameterSet>("services") : pset;
  auto const service_names = input.get_pset_names();
  for (auto const& service_name : service_names) {
    auto service_pset = input.get<fhicl::ParameterSet>(service_name);
    service_pset.put("service_type", service_name);
    post_processed_pset.put(service_name, service_pset);
  }

  static ArtServiceHelper helper{std::move(post_processed_pset)};
}
