////////////////////////////////////////////////////////////////////////
// Name:  FileCatalogMetadataDUNE_service.cc.  
//
// Purpose:  Art service adds microboone-specific per-job sam metadata.
//
//           This service does not have user-callable methods.  Simply
//           add to an art configuration in services.user block of job
//           file.
//
// Created:  3-Dec-2014,  T. Yang
// Copied FileCatalogMetadataMicroBooNE_service.cc by H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"

namespace util {

  // Class declaration.

  class FileCatalogMetadataDUNE
  {
  public:

    // Constructor, destructor.

    FileCatalogMetadataDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~FileCatalogMetadataDUNE() = default;

  private:

    // Callbacks.

    void postBeginJob();

    // Data members.

    std::string fMCGenerators;
    std::string fMCOscillationP;
    std::string fMCTriggerListVersion;
    std::string fMCBeamEnergy;
    std::string fMCBeamFluxID;
    std::string fMCName;
    std::string fMCDetectorType;
    std::string fMCNeutrinoFlavors;
    std::string fMCMassHierarchy;
    std::string fMCMiscellaneous;
    std::string fMCGeometryVersion;
    std::string fMCOverlay;
    std::string fDataRunMode;
    std::string fDataDetectorType;
    std::string fDataName;
    std::string fStageName;
  };

  //--------------------------------------------------------------------
  // Constructor.

  FileCatalogMetadataDUNE::
  FileCatalogMetadataDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  {
    // Get parameters.

    fMCGenerators = pset.get<std::string>("MCGenerators","");
    fMCOscillationP = pset.get<std::string>("MCOscillationP","");
    fMCTriggerListVersion = pset.get<std::string>("MCTriggerListVersion","");
    fMCBeamEnergy = pset.get<std::string>("MCBeamEnergy","");
    fMCBeamFluxID = pset.get<std::string>("MCBeamFluxID","");
    fMCName = pset.get<std::string>("MCName","");
    fMCDetectorType = pset.get<std::string>("MCDetectorType","");
    fMCNeutrinoFlavors = pset.get<std::string>("MCNeutrinoFlavors","");
    fMCMassHierarchy = pset.get<std::string>("MCMassHierarchy","");
    fMCMiscellaneous = pset.get<std::string>("MCMiscellaneous","");
    fMCGeometryVersion = pset.get<std::string>("MCGeometryVersion","");
    fMCOverlay = pset.get<std::string>("MCOverlay","");
    fDataRunMode = pset.get<std::string>("DataRunMode","");
    fDataDetectorType = pset.get<std::string>("DataDetectorType","");
    fDataName = pset.get<std::string>("DataName","");
    fStageName = pset.get<std::string>("StageName","");
    // Register for callbacks.

    reg.sPostBeginJob.watch(this, &FileCatalogMetadataDUNE::postBeginJob);
  }

  //--------------------------------------------------------------------
  // PostBeginJob callback.
  // Insert per-job metadata via FileCatalogMetadata service.
  void util::FileCatalogMetadataDUNE::postBeginJob()
  {
    // Get art metadata service.

    art::ServiceHandle<art::FileCatalogMetadata> mds;

    // Add metadata.

    if (fMCGenerators!="") mds->addMetadata("duneMCGenerators", fMCGenerators);
    if (fMCOscillationP!="") mds->addMetadata("duneMCOscillationP", fMCOscillationP);
    if (fMCTriggerListVersion!="") mds->addMetadata("duneMCTriggerListVersion", fMCTriggerListVersion);
    if (fMCBeamEnergy!="") mds->addMetadata("duneMCBeamEnergy", fMCBeamEnergy);
    if (fMCBeamFluxID!="") mds->addMetadata("duneMCBeamFluxID", fMCBeamFluxID);
    if (fMCName!="") mds->addMetadata("duneMCName", fMCName);
    if (fMCDetectorType!="") mds->addMetadata("duneMCDetectorType", fMCDetectorType);
    if (fMCNeutrinoFlavors!="") mds->addMetadata("duneMCNeutrinoFlavors", fMCNeutrinoFlavors);
    if (fMCMassHierarchy!="") mds->addMetadata("duneMCMassHierarchy", fMCMassHierarchy);
    if (fMCMiscellaneous!="") mds->addMetadata("duneMCMiscellaneous", fMCMiscellaneous);
    if (fMCGeometryVersion!="") mds->addMetadata("duneMCGeometryVersion", fMCGeometryVersion);
    if (fMCOverlay!="") mds->addMetadata("duneMCOverlay", fMCOverlay);
    if (fDataRunMode!="") mds->addMetadata("duneDataRunMode", fDataRunMode);
    if (fDataDetectorType!="") mds->addMetadata("duneDataDetectorType", fDataDetectorType);
    if (fDataName!="") mds->addMetadata("duneDataName", fDataName);
    if (fStageName!="") mds->addMetadata("StageName",fStageName);
  }

} // namespace util

DECLARE_ART_SERVICE(util::FileCatalogMetadataDUNE, LEGACY)
DEFINE_ART_SERVICE(util::FileCatalogMetadataDUNE)
