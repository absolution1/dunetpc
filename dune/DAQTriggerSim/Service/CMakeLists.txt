# dunetpc/dune/DAQTriggerSim/Service/CMakeLists.txt
#
# Michael Baird
# March 2016
#
# CMake instructions for DAQTriggerSim/Services


#
# Addplugin for each service.
#

simple_plugin(TemplateTriggerService   "service"
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${SEEDSERVICE_SERVICE}
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas_Utilities
                ${MF_MESSAGELOGGER}
                ${MF_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(ActivityTriggerService   "service"
                larcore_Geometry
                larcore_Geometry_Geometry_service
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${SEEDSERVICE_SERVICE}
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas_Utilities
                ${MF_MESSAGELOGGER}
                ${MF_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
		${ROOT_GEOM}
             )

simple_plugin(HardwareMapperService   "service"
                larcore_Geometry
                dune_Geometry
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${SEEDSERVICE_SERVICE}
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas_Utilities
                ${MF_MESSAGELOGGER}
                ${MF_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
             )


install_fhicl()
install_headers()
install_source()
