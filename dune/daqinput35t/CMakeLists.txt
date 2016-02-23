include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

ADD_SUBDIRECTORY(utilities)

art_make(BASENAME_ONLY
        LIBRARY_NAME lbne_daqinput35t
        LIB_LIBRARIES
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_PERSISTENCY_PROVENANCE}
        ${ART_UTILITIES}
        ${ROOT_CORE}
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${Boost_SYSTEM_LIBRARY}
        RunHistoryDUNE
        larcore_Geometry_Geometry_service
        larcore_Geometry
        lbne_daqinput35t_utilities
        lbne-raw-data_Overlays
        lardata_RawData
        lardata_RecoBase
        NO_PLUGINS
        )

simple_plugin(TpcMilliSliceDump "module"
lbne-raw-data_Overlays
larcore_Geometry_Geometry_service
larcore_Geometry
lardata_RawData
${SIMULATIONBASE}
# ${ART_FRAMEWORK_CORE}
# ${ART_FRAMEWORK_PRINCIPAL}
${ART_FRAMEWORK_SERVICES_REGISTRY}
${ART_FRAMEWORK_SERVICES_OPTIONAL}
${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
# ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
# ${ART_PERSISTENCY_COMMON}
# ${ART_PERSISTENCY_PROVENANCE}
# ${ART_UTILITIES}
# ${MF_MESSAGELOGGER}
# ${MF_UTILITIES}
# ${CETLIB}
${ROOT_HIST}
# ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(SSPDump "module"
lbne-raw-data_Overlays
larcore_Geometry_Geometry_service
larcore_Geometry
lardata_RawData
${SIMULATIONBASE}
art_Framework_Services_Registry
art_Framework_Services_Optional
art_Framework_Services_Optional_TFileService_service
${ROOT_HIST}
)

simple_plugin(SSPToOffline "module"
lbne-raw-data_Overlays
larcore_Geometry_Geometry_service
larcore_Geometry
lbne_daqinput35t
lardata_RawData
lardata_RecoBase
${SIMULATIONBASE}
${MF_MESSAGELOGGER}
)

simple_plugin(PTBToOffline "module"
lbne-raw-data_Overlays
larcore_Geometry_Geometry_service
larcore_Geometry
lbne_daqinput35t
lardata_RawData
${SIMULATIONBASE}
${MF_MESSAGELOGGER}
)

simple_plugin(SSPTriggerAna "module"
lbne-raw-data_Overlays
larcore_Geometry_Geometry_service
larcore_Geometry
lbne_daqinput35t
lardata_RawData
${SIMULATIONBASE}
${MF_MESSAGELOGGER}
art_Framework_Services_Registry
art_Framework_Services_Optional
art_Framework_Services_Optional_TFileService_service
${ROOT_BASIC_LIB_LIST}
${ROOT_HIST}
)

simple_plugin(SSPDiagnosticAna "module"
lbne-raw-data_Overlays
larcore_Geometry_Geometry_service
larcore_Geometry
lbne_daqinput35t
lardata_RawData
${SIMULATIONBASE}
${MF_MESSAGELOGGER}
art_Framework_Services_Registry
art_Framework_Services_Optional
art_Framework_Services_Optional_TFileService_service
${ROOT_HIST}
${ROOT_SPECTRUM}
)

simple_plugin(TpcDAQToOffline "module"
lbne-raw-data_Overlays
lbne_daqinput35t_utilities
lbne_daqinput35t
larcore_Geometry_Geometry_service
larcore_Geometry
lardata_RawData
)

simple_plugin(NoiseCorrelation "module"
			       lbne-raw-data_Overlays
			       lbne_daqinput35t_utilities
			       lbne_daqinput35t
			       larcore_Geometry_Geometry_service
			       larcore_Geometry
			       lardata_RawData
			       ${ART_FRAMEWORK_CORE}
			       ${ART_FRAMEWORK_PRINCIPAL}
			       ${ART_FRAMEWORK_SERVICES_REGISTRY}
			       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			       ${ART_PERSISTENCY_COMMON}
			       ${ART_PERSISTENCY_PROVENANCE}
			       ${ART_UTILITIES}
			       ${MF_MESSAGELOGGER}
			       ${MF_UTILITIES}
			       ${CETLIB}
			       ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(FilterWF "module"
 
 lbne_daqinput35t
 art_Framework_IO_Root
 art_Framework_IO_Sources
 lbne-raw-data_Overlays
 lbne_daqinput35t_utilities
 lbne_daqinput35t
 lardata_RawData
 ${ROOT_TREE}
 ${ROOT_RIO}
 ${MF_MESSAGELOGGER}
 art_Framework_Services_Registry
 art_Framework_Services_Optional
)

simple_plugin(CheckTiming "module"
			  lbne-raw-data_Overlays
			  lbne_daqinput35t_utilities
			  lbne_daqinput35t
			  larcore_Geometry_Geometry_service
			  larcore_Geometry
			  lardata_RawData
			  ${ART_FRAMEWORK_CORE}
			  ${ART_FRAMEWORK_PRINCIPAL}
			  ${ART_FRAMEWORK_SERVICES_REGISTRY}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			  ${ART_PERSISTENCY_COMMON}
			  ${ART_PERSISTENCY_PROVENANCE}
			  ${ART_UTILITIES}
			  ${MF_MESSAGELOGGER}
			  ${MF_UTILITIES}
			  ${CETLIB}
			  ${ROOT_BASIC_LIB_LIST}
		)

simple_plugin(SplitterInput "source"
 lbne_daqinput35t
 art_Framework_IO_Root
 art_Framework_IO_Sources
 larsim_Simulation
 ${ROOT_TREE}
 ${ROOT_RIO}
 ${MF_MESSAGELOGGER}
 ${SIMULATIONBASE}
 ${SIMULATION}
)

install_fhicl()
install_source()
install_headers()
