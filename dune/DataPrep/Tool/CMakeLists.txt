# dunetpc/dune/DataPrep/Tool/CMakeLists.txt

# David Adams
# June 2017

# Add plugin for each tool.

simple_plugin(AdcChannelDumper   "tool"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcChannelPlotter  "tool"
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcPedestalFitter  "tool"
                dune_ArtSupport
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcDataDumper   "tool"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcDataPlotter   "tool"
                dune_DuneCommon
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

simple_plugin(AdcUnderflowFlagger   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcThresholdSignalFinder   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AdcRoiViewer   "tool"
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(AcdLocalRoiBuilder   "tool"
                dune_ArtSupport
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

simple_plugin(AcdDPhase3x1x1RoiBuilder   "tool"
                dune_ArtSupport
                dune_DuneCommon
                lardataobj_RawData
                lardataobj_RecoBase
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)

