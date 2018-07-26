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
                nutools_RandomUtils_NuRandomService_service
                ${ART_FRAMEWORK_CORE}
                art_Utilities canvas
                ${MF_MESSAGELOGGER}
                
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(HardwareMapperService   "service"
                larcorealg_Geometry
                dune_Geometry
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_PRINCIPAL}
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${MF_MESSAGELOGGER}
                
                ${ROOT_BASIC_LIB_LIST}
                art_Utilities canvas
             )


install_fhicl()
install_headers()
install_source()
