# dunetpc/dune/DAQTriggerSim/Service/test/CMakeLists.txt
#
# Michael Baird
# March 2016
#
# CMake instructions for DAQTriggerSim Services (the test module)

# Addplugin for each service.

simple_plugin(TriggerTestService   "service"
                ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                ${SEEDSERVICE_SERVICE}
                ${ART_FRAMEWORK_CORE}
                ${ART_UTILITIES}
                ${MF_MESSAGELOGGER}
                ${MF_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
             )
