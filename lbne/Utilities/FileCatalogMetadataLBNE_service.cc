////////////////////////////////////////////////////////////////////////
// Name:  FileCatalogMetadataLBNE_service.cc.  
//
// Purpose:  Art service adds microboone-specific per-job sam metadata.
//
//           FCL parameters:
//
//           FCLName        - FCL file name.
//           FCLVersion     - FCL file version.
//           ProjectName    - Project name.
//           ProjectStage   - Project stage.
//           ProjectVersion - Project version.
//
//           Above values are recorded in internal sam metadata generated
//           by art program.
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

  class FileCatalogMetadataLBNE
  {
  public:

    // Constructor, destructor.

    FileCatalogMetadataLBNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~FileCatalogMetadataLBNE() = default;

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
  };

  //--------------------------------------------------------------------
  // Constructor.

  FileCatalogMetadataLBNE::
  FileCatalogMetadataLBNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
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

    // Register for callbacks.

    reg.sPostBeginJob.watch(this, &FileCatalogMetadataLBNE::postBeginJob);
  }

  //--------------------------------------------------------------------
  // PostBeginJob callback.
  // Insert per-job metadata via FileCatalogMetadata service.
  void util::FileCatalogMetadataLBNE::postBeginJob()
  {
    // Get art metadata service.

    art::ServiceHandle<art::FileCatalogMetadata> mds;

    // Add metadata.

    mds->addMetadata("lbneMCGenerators", fMCGenerators);
    mds->addMetadata("lbneMCOscillationP", fMCOscillationP);
    mds->addMetadata("lbneMCTriggerListVersion", fMCTriggerListVersion);
    mds->addMetadata("lbneMCBeamEnergy", fMCBeamEnergy);
    mds->addMetadata("lbneMCBeamFluxID", fMCBeamFluxID);
    mds->addMetadata("lbneMCName", fMCName);
    mds->addMetadata("lbneMCDetectorType", fMCDetectorType);
    mds->addMetadata("lbneMCNeutrinoFlavors", fMCNeutrinoFlavors);
    mds->addMetadata("lbneMCMassHierarchy", fMCMassHierarchy);
    mds->addMetadata("lbneMCMiscellaneous", fMCMiscellaneous);
    mds->addMetadata("lbneMCGeometryVersion", fMCGeometryVersion);
    mds->addMetadata("lbneMCOverlay", fMCOverlay);
    mds->addMetadata("lbneDataRunMode", fDataRunMode);
    mds->addMetadata("lbneDataDetectorType", fDataDetectorType);
    mds->addMetadata("lbneDataName", fDataName);
  }

} // namespace util

DECLARE_ART_SERVICE(util::FileCatalogMetadataLBNE, LEGACY)
DEFINE_ART_SERVICE(util::FileCatalogMetadataLBNE)
