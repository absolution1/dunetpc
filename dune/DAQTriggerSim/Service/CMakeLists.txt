# dunetpc/dune/DAQTriggerSim/Service/CMakeLists.txt
#
# Michael Baird
# March 2016
#
# CMake instructions for DAQTriggerSim/Services

cet_build_plugin(TemplateTriggerService   art::service
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
                nurandom_RandomUtils_NuRandomService_service
                ART_FRAMEWORK_CORE
                ART_FRAMEWORK_PRINCIPAL
                art_Utilities canvas
                messagefacility::MF_MessageLogger
                cetlib cetlib_except
                ROOT_BASIC_LIB_LIST
             )

cet_build_plugin(HardwareMapperService   art::service
                larcorealg_Geometry
                dune_Geometry
                ART_FRAMEWORK_CORE
                ART_FRAMEWORK_PRINCIPAL
                ART_FRAMEWORK_SERVICES_REGISTRY
                ART_ROOT_IO_TFILE_SUPPORT
                ROOT::Core
                ART_PERSISTENCY_PROVENANCE
                messagefacility::MF_MessageLogger
                ROOT_BASIC_LIB_LIST
                art_Utilities canvas
             )

install_fhicl()
install_headers()
install_source()
