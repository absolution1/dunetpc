#include "art/Framework/Services/Registry/ServiceRegistry.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/intermediate_table.h"

#include "dune/ArtSupport/ArtServiceHelper.h"

ArtServiceHelper::ArtServiceHelper(fhicl::ParameterSet&& pset) :
  activityRegistry_{},
  sharedResources_{},
  servicesManager_{std::move(pset), activityRegistry_, sharedResources_}
{
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
  load_services(fhicl::ParameterSet::make
                (fhicl::parse_document(config, lookup)));
}

void ArtServiceHelper::load_services(std::string const& filename, FileOnPath_t)
{
  cet::filepath_lookup lookup{"FHICL_FILE_PATH"};
  load_services(fhicl::ParameterSet::make(filename, lookup));
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
