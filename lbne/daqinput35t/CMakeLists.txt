include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

ADD_SUBDIRECTORY(utilities)

art_make(LIB_LIBRARIES
	${ART_FRAMEWORK_SERVICES_REGISTRY}
	${ART_PERSISTENCY_PROVENANCE}
	${ART_UTILITIES}
	${ROOT_CORE}
	Geometry_service
	Geometry
	lbne_daqinput35t_utilities
	lbne-raw-data_Overlays
        RawData
        NO_PLUGINS
        )

simple_plugin(TpcMilliSliceDump "module"
lbne-raw-data_Overlays
Geometry_service
Geometry
RawData
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
Geometry_service
Geometry
RawData
${SIMULATIONBASE}
art_Framework_Services_Registry
art_Framework_Services_Optional
art_Framework_Services_Optional_TFileService_service
${ROOT_HIST}
)

simple_plugin(SSPToOffline "module"
lbne-raw-data_Overlays
Geometry_service
Geometry
RawData
${SIMULATIONBASE}
# ${ART_FRAMEWORK_CORE}
# ${ART_FRAMEWORK_PRINCIPAL}
# ${ART_FRAMEWORK_SERVICES_REGISTRY}
# ${ART_FRAMEWORK_SERVICES_OPTIONAL}
# ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
# ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
# ${ART_PERSISTENCY_COMMON}
# ${ART_PERSISTENCY_PROVENANCE}
# ${ART_UTILITIES}
${MF_MESSAGELOGGER}
# ${MF_UTILITIES}
# ${CETLIB}
# ${ROOT_BASIC_LIB_LIST}
)

simple_plugin(TpcDAQToOffline "module"
lbne-raw-data_Overlays
lbne_daqinput35t_utilities
lbne_daqinput35t
Geometry_service
Geometry
RawData
)

simple_plugin(SplitterInput "source"
 lbne_daqinput35t
 art_Framework_IO_Root
 art_Framework_IO_Sources
 ${ROOT_TREE}
 ${ROOT_RIO}
 ${MF_MESSAGELOGGER}
)

install_fhicl()
install_source()
install_headers()
