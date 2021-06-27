# dunetpc/dune/DataPrep/Tool/CMakeLists.txt

# David Adams
# June 2017

# Add plugin for each tool.

simple_plugin(StandardAdcChannelStringTool   "tool"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcResultDumper   "tool"
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcChannelDumper   "tool"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcChannelTrimmer  "tool"
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcChannelSplitter  "tool"
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcChannelPlotter  "tool"
                dune_ArtSupport
                dune_DuneCommon
		${ART_FRAMEWORK_SERVICES_REGISTRY}
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcPedestalFitter  "tool"
                dune_ArtSupport
                dune_DuneCommon
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcCodeMitigator  "tool"
                dune_ArtSupport
                dune_DuneCommon
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcSampleScaler   "tool"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcRangeSampleScaler   "tool"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcTimingMitigator  "tool"
                dune_ArtSupport
                dune_DuneCommon
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcDataDumper   "tool"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcTickModViewer   "tool"
                dune_DataPrep_Utility
                dune_ArtSupport
                dune_DuneCommon
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcDataPlotter   "tool"
                dune_ArtSupport
                dune_DuneCommon
		${ART_FRAMEWORK_SERVICES_REGISTRY}
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcDetectorPlotter   "tool"
                dune_ArtSupport
                dune_DuneCommon
                dune_Geometry
		${ART_FRAMEWORK_SERVICES_REGISTRY}
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcChannelMetric   "tool"
                dune_DuneCommon
                dune_ArtSupport
		${ART_FRAMEWORK_SERVICES_REGISTRY}
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcEventViewer   "tool"
                dune_DuneCommon
                dune_ArtSupport
		${ART_FRAMEWORK_SERVICES_REGISTRY}
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AcdDigitReader   "tool"
                dune_DuneCommon
                lardataobj_RawData
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AcdWireReader   "tool"
                dune_DuneCommon
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcSampleFiller   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(SubtractBaseline   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(UndershootCorr   "tool"
              dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        art_Utilities canvas
                        ${CETLIB_LIBS}
                        ${MF_MESSAGELOGGER}
                        ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ExpTailRemover   "tool"
                dune_DuneCommon
                dune_ArtSupport
                dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                art_Utilities canvas
                ${CETLIB_LIBS}
                ${MF_MESSAGELOGGER}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ExpTailPedRemover   "tool"
                dune_DuneCommon
                dune_ArtSupport
                dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                art_Utilities canvas
                ${CETLIB_LIBS}
                ${MF_MESSAGELOGGER}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(VintageDeconvoluter   "tool"
                SignalShapingServiceDUNE_service
                art_Framework_Services_Registry
                art_Utilities
                canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(PddpVintageDeconvoluter   "tool"
                SignalShapingServiceDUNEDPhase_service
		lardata_Utilities
                art_Framework_Services_Registry
                art_Utilities
                canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )


simple_plugin(AdcUnderflowFlagger   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(FclStickyCodeFlagger   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcKeepAllSignalFinder   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcThresholdSignalFinder   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcNoiseSignalFinder   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcRoiShifter   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(DuneAdcSignalFinder   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcRegularSignalFinder   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcRoiViewer   "tool"
                dune_ArtSupport
                dune_DuneCommon
		canvas
                cetlib 
		cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcDPhase3x1x1LocalRoiBuilder   "tool"
                dune_DuneCommon
                lardataobj_RawData
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(FembLinearCalibration   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(FloatArrayGainCalibration   "tool"
                dune_ArtSupport
		canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(FembMappedAdcModifier   "tool"
                dune_ArtSupport
		canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcChannelFFT   "tool"
                dune_ArtSupport
                dune_DuneCommon
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcDeconvoluteFFT   "tool"
                dune_ArtSupport
                dune_DuneCommon
		canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )
simple_plugin(Adc2dConvolute   "tool"
                dune_ArtSupport
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcChannelDftPlotter   "tool"
                dune_DataPrep_Utility
                dune_DuneCommon
                dune_ArtSupport
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcDPhase3x1x1RoiBuilder   "tool"
                dune_ArtSupport
                dune_DuneCommon
                lardataobj_RawData
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcMultiThreshSignalFinder   "tool"
                dune_ArtSupport
                dune_DuneCommon
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(PdspNoiseRemoval   "tool"
               dune_ArtSupport
               larcorealg_Geometry
               dune_DuneCommon
               lardataobj_RawData
               lardataobj_RecoBase
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               art_Utilities canvas
               cetlib cetlib_except
               dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
               ${ROOT_BASIC_LIB_LIST}
               ${ROOT_SPECTRUM}
             )

simple_plugin(RemoveBadChannels   "tool"
                dune_ArtSupport
                dune_DuneCommon
		${ART_FRAMEWORK_SERVICES_REGISTRY}
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcRoiSlicer   "tool"
                dune_ArtSupport
                dune_DuneCommon
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(BaselineDetrend   "tool"
                dune_ArtSupport
                dune_DuneCommon
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcRoiToTree   "tool"
                dune_ArtSupport
                dune_DuneCommon
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)
install_source()

