# dunetpc/dune/DataPrep/Service/CMakeLists.txt

# David Adams
# May 2016

# Add plugin for each service.

include_directories(AFTER $ENV{VALGRIND_INC})

cet_build_plugin(ConfigurableAdcChannelDataCopyService   art::service
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(InterpolatingAdcMitigationService   art::service
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(StandardRawDigitExtractService   art::service
                dune_ArtSupport
                lardataobj_RawData
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(ThresholdNoiseRemovalService   art::service
                larcorealg_Geometry
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(MultiChannelNoiseRemovalService   art::service
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(DuneDPhase3x1x1NoiseRemovalService   art::service
                dune_ArtSupport
                larcorealg_Geometry
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
                ROOT_FFTW
                lardataobj_RawData
             )

cet_build_plugin(MedianPedestalService   art::service
		ART_FRAMEWORK_SERVICES_REGISTRY
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcSuppressSignalFindingService art::service
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(DuneDeconvolutionService art::service
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
                SignalShapingServiceDUNE_service
             )

cet_build_plugin(DuneDPhaseDeconvolutionService art::service
		lardata_Utilities
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
                SignalShapingServiceDUNEDPhase_service
             )

cet_build_plugin(KeepAllRoiBuildingService art::service
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(DuneRoiBuildingService art::service
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
                SignalShapingServiceDUNE_service
             )

cet_build_plugin(DuneDPhaseRoiBuildingService art::service
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
                SignalShapingServiceDUNEDPhase_service
             )

cet_build_plugin(StandardAdcWireBuildingService art::service
                lardataobj_RecoBase
                lardata_ArtDataHelper
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(ToolBasedRawDigitPrepService   art::service
                dune_ArtSupport
                lardataobj_RawData
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(TpcToolBasedRawDigitPrepService   art::service
                dune_DuneInterface_Data
                dune_ArtSupport
                lardataobj_RawData
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(StandardRawDigitPrepService   art::service
                dune_ArtSupport
                lardataobj_RawData
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

add_subdirectory(test)

install_source()
