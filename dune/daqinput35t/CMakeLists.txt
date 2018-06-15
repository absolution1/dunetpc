include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

ADD_SUBDIRECTORY(utilities)
ADD_SUBDIRECTORY(test)

art_make(BASENAME_ONLY
        LIBRARY_NAME lbne_daqinput35t
        LIB_LIBRARIES
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        art_Persistency_Provenance canvas
        art_Utilities canvas
        ${ROOT_CORE}
        ${MF_MESSAGELOGGER}
        
        ${Boost_SYSTEM_LIBRARY}
        RunHistoryDUNE
        larcore_Geometry_Geometry_service
        larcorealg_Geometry
        lbne_daqinput35t_utilities
        lbne-raw-data_Overlays
	lbne-raw-data_Services_ChannelMap_ChannelMapService_service
        lardataobj_RawData
        lardataobj_RecoBase
        dune_DuneCommon
        NO_PLUGINS
        )

simple_plugin(TpcMilliSliceDump "module"
lbne-raw-data_Overlays
larcore_Geometry_Geometry_service
larcorealg_Geometry
lardataobj_RawData
nusimdata_SimulationBase
# ${ART_FRAMEWORK_CORE}
# ${ART_FRAMEWORK_PRINCIPAL}
${ART_FRAMEWORK_SERVICES_REGISTRY}
${ART_FRAMEWORK_SERVICES_OPTIONAL}
${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
# ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
# art_Persistency_Common canvas
# art_Persistency_Provenance canvas
# art_Utilities canvas
# ${MF_MESSAGELOGGER}
# 
cetlib cetlib_except
${ROOT_HIST}
# ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(SSPDump "module"
lbne-raw-data_Overlays
larcore_Geometry_Geometry_service
larcorealg_Geometry
lardataobj_RawData
nusimdata_SimulationBase
art_Framework_Services_Registry
art_Framework_Services_Optional
art_Framework_Services_Optional_TFileService_service
cetlib cetlib_except
${ROOT_HIST}
)

simple_plugin(SSPToOffline "module"
lbne-raw-data_Overlays
larcore_Geometry_Geometry_service
larcorealg_Geometry
lbne_daqinput35t
lardataobj_RawData
lardataobj_RecoBase
nusimdata_SimulationBase
${MF_MESSAGELOGGER}
cetlib cetlib_except
)

simple_plugin(PTBToOffline "module"
lbne-raw-data_Overlays
larcore_Geometry_Geometry_service
larcorealg_Geometry
lbne_daqinput35t
lardataobj_RawData
nusimdata_SimulationBase
${MF_MESSAGELOGGER}
cetlib cetlib_except
)

simple_plugin(SSPTriggerAna "module"
lbne-raw-data_Overlays
larcore_Geometry_Geometry_service
larcorealg_Geometry
lbne_daqinput35t
lardataobj_RawData
nusimdata_SimulationBase
${MF_MESSAGELOGGER}
cetlib cetlib_except
art_Framework_Services_Registry
art_Framework_Services_Optional
art_Framework_Services_Optional_TFileService_service
${ROOT_BASIC_LIB_LIST}
${ROOT_HIST}
)

simple_plugin(SSPDiagnosticAna "module"
lbne-raw-data_Overlays
larcore_Geometry_Geometry_service
larcorealg_Geometry
lbne_daqinput35t
lardataobj_RawData
nusimdata_SimulationBase
${MF_MESSAGELOGGER}
cetlib cetlib_except
art_Framework_Services_Registry
art_Framework_Services_Optional
art_Framework_Services_Optional_TFileService_service
${ROOT_HIST}
${ROOT_SPECTRUM}
)

simple_plugin(TpcDAQToOffline "module"
lbne-raw-data_Overlays
lbne-raw-data_Services_ChannelMap_ChannelMapService_service
lbne_daqinput35t_utilities
lbne_daqinput35t
larcore_Geometry_Geometry_service
larcorealg_Geometry
lardataobj_RawData
art_Framework_Services_Registry
art_Framework_Services_Optional
${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
${FHICLCPP}
cetlib cetlib_except
${ROOT_BASIC_LIB_LIST}
)

simple_plugin(NoiseCorrelation "module"
			       lbne-raw-data_Overlays
			       lbne_daqinput35t_utilities
			       lbne_daqinput35t
			       larcore_Geometry_Geometry_service
			       larcorealg_Geometry
			       lardataobj_RawData
			       ${ART_FRAMEWORK_CORE}
			       ${ART_FRAMEWORK_PRINCIPAL}
			       ${ART_FRAMEWORK_SERVICES_REGISTRY}
			       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			       art_Persistency_Common canvas
			       art_Persistency_Provenance canvas
			       art_Utilities canvas
			       ${MF_MESSAGELOGGER}
			       
			       cetlib cetlib_except
			       ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(FilterWF "module"
 
 lbne_daqinput35t
 art_Framework_IO_Root
 art_Framework_IO_Sources
 lbne-raw-data_Overlays
 lbne_daqinput35t_utilities
 lbne_daqinput35t
 lardataobj_RawData
 ${ROOT_TREE}
 ${ROOT_RIO}
 ${MF_MESSAGELOGGER}
 cetlib cetlib_except
 art_Framework_Services_Registry
 art_Framework_Services_Optional
 ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
)

simple_plugin(CheckRawDigit "module"
 
 lbne_daqinput35t
 art_Framework_IO_Root
 art_Framework_IO_Sources
 lbne-raw-data_Overlays
 lbne_daqinput35t_utilities
 lbne_daqinput35t
 lardataobj_RawData
 ${ROOT_TREE}
 ${ROOT_RIO}
 ${MF_MESSAGELOGGER}
 cetlib cetlib_except
 art_Framework_Services_Registry
 art_Framework_Services_Optional
 ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
)

simple_plugin(Sigmoidfilter "module"

 lbne_daqinput35t
 art_Framework_IO_Root
 art_Framework_IO_Sources
 lbne-raw-data_Overlays
 lbne_daqinput35t_utilities
 lbne_daqinput35t
 lardataobj_RawData
 ${ROOT_TREE}
 ${ROOT_RIO}
 ${MF_MESSAGELOGGER}
 cetlib cetlib_except
 art_Framework_Services_Registry
 art_Framework_Services_Optional
 ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
)

simple_plugin(FilterAnalyzer "module"
 lbne_daqinput35t
 art_Framework_IO_Root
 art_Framework_IO_Sources
 lbne-raw-data_Overlays
 lbne_daqinput35t_utilities
 lbne_daqinput35t
 lardataobj_RawData
 ${ROOT_TREE}
 ${ROOT_RIO}
 ${MF_MESSAGELOGGER}
 cetlib cetlib_except
 art_Framework_Services_Registry
 art_Framework_Services_Optional
 ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
)

simple_plugin(FilterAnalyzerRunList "module"
 lbne_daqinput35t
 art_Framework_IO_Root
 art_Framework_IO_Sources
 lbne-raw-data_Overlays
 lbne_daqinput35t_utilities
 lbne_daqinput35t
 lardataobj_RawData
 ${ROOT_TREE}
 ${ROOT_RIO}
 ${MF_MESSAGELOGGER}
 art_Framework_Services_Registry
 cetlib cetlib_except
 art_Framework_Services_Optional
 ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
)

simple_plugin(GoodRun "module"
 lbne-raw-data_Overlays
 lbne_daqinput35t_utilities
 lbne_daqinput35t
 lardataobj_RawData
 ${ART_FRAMEWORK_CORE}
 ${ART_FRAMEWORK_PRINCIPAL}
 ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 cetlib cetlib_except
 ${ROOT_BASIC_LIB_LIST}
 ${ROOT_CORE}
)

simple_plugin(BadTimings "module"
 lbne-raw-data_Overlays
 lbne_daqinput35t_utilities
 lbne_daqinput35t
 lardataobj_RawData
 ${ART_FRAMEWORK_CORE}
 ${ART_FRAMEWORK_PRINCIPAL}
 ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 cetlib cetlib_except
 ${ROOT_BASIC_LIB_LIST}
 ${ROOT_CORE}
)

simple_plugin(SplitterInput "source"
 lbne_daqinput35t
 art_Framework_IO_Root
 art_Framework_IO_Sources
 larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
 ${ROOT_TREE}
 ${ROOT_RIO}
 ${MF_MESSAGELOGGER}
 cetlib cetlib_except
 nusimdata_SimulationBase
 ${SIMULATION}
)

install_fhicl()
install_source()
install_headers()
