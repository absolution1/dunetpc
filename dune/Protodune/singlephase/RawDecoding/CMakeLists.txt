

simple_plugin(TimingRawDecoder "module"
                        lardataobj_RawData
                        dune-raw-data_Overlays
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 #                       art_Persistency_Common canvas_Persistency_Common
 #                       art_Persistency_Provenance canvas_Persistency_Provenance
 #                       art_Utilities canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        BASENAME_ONLY
)

simple_plugin(SSPRawDecoder "module"
                        lardataobj_RawData
                        dune-raw-data_Overlays
			artdaq-core_Data
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 #                       art_Persistency_Common canvas_Persistency_Common
 #                       art_Persistency_Provenance canvas_Persistency_Provenance
 #                       art_Utilities canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        BASENAME_ONLY
)


include_directories( "$ENV{DUNEPDSPRCE_INC}" ) 
cet_find_library( RCEDAMLIB NAMES protodune-dam PATHS ENV DUNEPDSPRCE_LIB NO_DEFAULT_PATH )

# a cmake version of a test in JJ Russell's define_tools.mk.  Unfortunately this is needed
# whenever anyone includes protodune dam lib includes.
# haven't tested this with the clang compiler as we distribute gcc on mac

# it's gone in v0.0.4 of JJ's product.  Keep it here for v3's sake.

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCLASS_MEMBER_UNUSED=unused")
else(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCLASS_MEMBER_UNUSED=")
endif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")


simple_plugin(RceRawDecoder "module"
                        lardataobj_RawData
                        dune-raw-data_Overlays
			artdaq-core_Data
                        dune-raw-data_Services_ChannelMap_PdChannelMapService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 #                       art_Persistency_Common canvas_Persistency_Common
 #                       art_Persistency_Provenance canvas_Persistency_Provenance
 #                       art_Utilities canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${RCEDAMLIB}
                        BASENAME_ONLY
)

simple_plugin(FelixRawDecoder "module"
                        lardataobj_RawData
                        dune-raw-data_Overlays
                        artdaq-core_Data
                        dune-raw-data_Services_ChannelMap_PdChannelMapService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 #                       art_Persistency_Common canvas_Persistency_Common
 #                       art_Persistency_Provenance canvas_Persistency_Provenance
 #                       art_Utilities canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        ${RCEDAMLIB}
                        BASENAME_ONLY
)

art_make_library( LIBRARY_NAME dunetpc_protodune_RawDecodingUtils
                  SOURCE
                  RawDecodingUtils.cc             
                  LIBRARIES
                  lardataobj_RawData
                  dune-raw-data_Overlays
                  artdaq-core_Data
                  ${ART_FRAMEWORK_CORE}
                  ${ART_FRAMEWORK_PRINCIPAL}
                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 #                 art_Persistency_Common canvas_Persistency_Common
 #                 art_Persistency_Provenance canvas_Persistency_Provenance
 #                 art_Utilities canvas_Utilities
                  ${MF_MESSAGELOGGER}
                  ${MF_UTILITIES}
                  ${CETLIB}
                  ${ROOT_BASIC_LIB_LIST}
                  ${RCEDAMLIB}
)


install_headers()
install_fhicl()
install_source()
install_scripts()
