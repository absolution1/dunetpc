////////////////////////////////////////////////////////////////////////
// Name:  FileCatalogMetadataDUNE_service.cc.  
//
// Purpose:  Art service adds dune-specific per-job sam metadata.
//
//           This service does not have user-callable methods.  Simply
//           add to an art configuration in services.user block of job
//           file.
//
// Created:  1-Nov-2017
//  split from the original by T. Yang so we can include this
// Copied FileCatalogMetadataMicroBooNE_service.cc by H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#ifndef FILECATALOGMETADATADUNE_H
#define FILECATALOGMETADATADUNE_H

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

    // accessors

    const std::string&  MCGenerators()		const   { return fMCGenerators; }		  
    const std::string&  MCOscillationP()        const   { return fMCOscillationP; }	  
    const std::string&  MCTriggerListVersion()	const   { return fMCTriggerListVersion; }	  
    const std::string&  MCBeamEnergy()		const   { return fMCBeamEnergy; }		  
    const std::string&  MCBeamFluxID()		const   { return fMCBeamFluxID; }		  
    const std::string&  MCName()	        const   { return fMCName; }		  
    const std::string&  MCDetectorType()        const   { return fMCDetectorType; }	  
    const std::string&  MCNeutrinoFlavors()	const   { return fMCNeutrinoFlavors; }	  
    const std::string&  MCMassHierarchy()       const   { return fMCMassHierarchy; }	  
    const std::string&  MCMiscellaneous()       const   { return fMCMiscellaneous; }	  
    const std::string&  MCGeometryVersion()	const   { return fMCGeometryVersion; }	  
    const std::string&  MCOverlay()		const   { return fMCOverlay; }		  
    const std::string&  DataRunMode()		const   { return fDataRunMode; }		  
    const std::string&  DataDetectorType()	const   { return fDataDetectorType; }	  
    const std::string&  DataName()		const   { return fDataName; }		  
    const std::string&  StageName()             const   { return fStageName; }               

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


} // namespace util

DECLARE_ART_SERVICE(util::FileCatalogMetadataDUNE, LEGACY)

#endif
