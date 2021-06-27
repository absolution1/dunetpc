# dunetpc/dune/DataPrep/Service/CMakeLists.txt

# David Adams
# May 2016

# Add plugin for each service.

include_directories(AFTER $ENV{VALGRIND_INC})

simple_plugin(ConfigurableAdcChannelDataCopyService   "service"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(InterpolatingAdcMitigationService   "service"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(StandardRawDigitExtractService   "service"
                dune_ArtSupport
                lardataobj_RawData
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ThresholdNoiseRemovalService   "service"
                larcorealg_Geometry
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(MultiChannelNoiseRemovalService   "service"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(DuneDPhase3x1x1NoiseRemovalService   "service"
                dune_ArtSupport
                larcorealg_Geometry
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
                ${ROOT_FFTW}
                lardataobj_RawData
             )

simple_plugin(MedianPedestalService   "service"
		${ART_FRAMEWORK_SERVICES_REGISTRY}
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcSuppressSignalFindingService "service"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(DuneDeconvolutionService "service"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
                SignalShapingServiceDUNE_service
             )

simple_plugin(DuneDPhaseDeconvolutionService "service"
		lardata_Utilities
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
                SignalShapingServiceDUNEDPhase_service
             )

simple_plugin(KeepAllRoiBuildingService "service"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(DuneRoiBuildingService "service"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
                SignalShapingServiceDUNE_service
             )

simple_plugin(DuneDPhaseRoiBuildingService "service"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
                SignalShapingServiceDUNEDPhase_service
             )

simple_plugin(StandardAdcWireBuildingService "service"
                lardataobj_RecoBase
                lardata_ArtDataHelper
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ToolBasedRawDigitPrepService   "service"
                dune_ArtSupport
                lardataobj_RawData
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(TpcToolBasedRawDigitPrepService   "service"
                dune_DuneInterface_Data
                dune_ArtSupport
                lardataobj_RawData
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(StandardRawDigitPrepService   "service"
                dune_ArtSupport
                lardataobj_RawData
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)

install_source()
