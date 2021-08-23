# dunetpc/dune/DataPrep/Tool/CMakeLists.txt

# David Adams
# June 2017

# Add plugin for each tool.

cet_build_plugin(AdcThresholdSignalFinder dune::TpcDataTool
                 LIBRARIES PRIVATE
                 ROOT::Core
                 canvas::canvas
                 art::Utilities
                 dune::ArtSupport
                 dune::DuneCommon_Utility)

cet_build_plugin(StandardAdcChannelStringTool   art::tool
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcResultDumper   art::tool
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcChannelDumper   art::tool
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcChannelTrimmer  art::tool
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcChannelSplitter  art::tool
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcChannelPlotter  art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
		ART_FRAMEWORK_SERVICES_REGISTRY
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcPedestalFitter  art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcCodeMitigator  art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcSampleScaler   art::tool
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcRangeSampleScaler   art::tool
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcTimingMitigator  art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcDataDumper   art::tool
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcTickModViewer   art::tool
                dune_DataPrep_Utility
                dune_ArtSupport
                dune_DuneCommon_Utility
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcDataPlotter   art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
		ART_FRAMEWORK_SERVICES_REGISTRY
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcDetectorPlotter   art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
                dune_Geometry
		ART_FRAMEWORK_SERVICES_REGISTRY
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcChannelMetric   art::tool
                dune_DuneCommon_Utility
                dune_ArtSupport
		ART_FRAMEWORK_SERVICES_REGISTRY
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcEventViewer   art::tool
                dune_DuneCommon_Utility
                dune_ArtSupport
		ART_FRAMEWORK_SERVICES_REGISTRY
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AcdDigitReader   art::tool
                dune_DuneCommon_Utility
                lardataobj_RawData
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AcdWireReader   art::tool
                dune_DuneCommon_Utility
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcSampleFiller   art::tool
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(SubtractBaseline   art::tool
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(UndershootCorr   art::tool
                        dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
                        art_Utilities canvas
                        messagefacility::MF_MessageLogger
                        ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(ExpTailRemover   art::tool
                dune_DuneCommon_Utility
                dune_ArtSupport
                dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                ART_FRAMEWORK_CORE
                ART_FRAMEWORK_PRINCIPAL
                ART_FRAMEWORK_SERVICES_REGISTRY
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                ART_ROOT_IO_TFILESERVICE_SERVICE
                art_Utilities canvas
                messagefacility::MF_MessageLogger
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(ExpTailPedRemover   art::tool
                dune_DuneCommon_Utility
                dune_ArtSupport
                dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                ART_FRAMEWORK_CORE
                ART_FRAMEWORK_PRINCIPAL
                ART_FRAMEWORK_SERVICES_REGISTRY
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                ART_ROOT_IO_TFILESERVICE_SERVICE
                art_Utilities canvas
                messagefacility::MF_MessageLogger
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(VintageDeconvoluter   art::tool
                SignalShapingServiceDUNE_service
                art_Framework_Services_Registry
                art_Utilities
                canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(PddpVintageDeconvoluter   art::tool
                SignalShapingServiceDUNEDPhase_service
		lardata_Utilities
                art_Framework_Services_Registry
                art_Utilities
                canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )


cet_build_plugin(AdcUnderflowFlagger   art::tool
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(FclStickyCodeFlagger   art::tool
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcKeepAllSignalFinder   art::tool
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcNoiseSignalFinder   art::tool
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcRoiShifter   art::tool
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(DuneAdcSignalFinder   art::tool
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcRegularSignalFinder   art::tool
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcRoiViewer   art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
		canvas
                cetlib 
		cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcDPhase3x1x1LocalRoiBuilder   art::tool
                dune_DuneCommon_Utility
                lardataobj_RawData
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(FembLinearCalibration   art::tool
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(FloatArrayGainCalibration   art::tool
                dune_DuneCommon_Utility
                dune_ArtSupport
		canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(FembMappedAdcModifier   art::tool
                dune_ArtSupport
		canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcChannelFFT   art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcDeconvoluteFFT   art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
		canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )
cet_build_plugin(Adc2dConvolute   art::tool
                dune_ArtSupport
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcChannelDftPlotter   art::tool
                dune_DataPrep_Utility
                dune_DuneCommon_Utility
                dune_ArtSupport
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcDPhase3x1x1RoiBuilder   art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
                lardataobj_RawData
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcMultiThreshSignalFinder   art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(PdspNoiseRemoval   art::tool
               dune_ArtSupport
               larcorealg_Geometry
               dune_DuneCommon_Utility
               lardataobj_RawData
               lardataobj_RecoBase
               ART_FRAMEWORK_SERVICES_REGISTRY
               art_Utilities canvas
               cetlib cetlib_except
               dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
               ROOT_BASIC_LIB_LIST
               ROOT_SPECTRUM
             )

cet_build_plugin(RemoveBadChannels   art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
		ART_FRAMEWORK_SERVICES_REGISTRY
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcRoiSlicer   art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
                art_Utilities canvas
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(BaselineDetrend   art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(AdcRoiToTree   art::tool
                dune_ArtSupport
                dune_DuneCommon_Utility
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

add_subdirectory(test)
install_source()
