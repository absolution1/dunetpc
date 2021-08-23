include_directories( $ENV{LIBWDA_FQ_DIR}/include )
include_directories( $ENV{IFDH_ART_INC} )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )
include_directories( $ENV{IFBEAM_FQ_DIR}/include )

cet_find_library( WDALIB NAMES wda PATHS ENV LIBWDA_LIB NO_DEFAULT_PATH )
cet_find_library( IFBEAMSERVICE NAMES IFBeam_service PATHS ENV IFDH_ART_LIB NO_DEFAULT_PATH )
cet_find_library( IFBEAM NAMES ifbeam PATHS ENV IFBEAM_LIB   NO_DEFAULT_PATH )


art_make(  MODULE_LIBRARIES
                        dune_ArtSupport
                        IFBEAM
                        IFBEAMSERVICE
                        lardataobj_RawData
                        lardataobj_RecoBase
                        larreco_RecoAlg
                        larreco_Calorimetry
                        larreco_RecoAlg_PMAlg
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_ArtDataHelper
                        lardata_Utilities
                        larsim_Simulation 
                        lardataobj_Simulation
                        dune_DuneObj
                        dune_raw_data::dune-raw-data_Services_ChannelMap_PdspChannelMapService_service
                        larsim_MCCheater_PhotonBackTrackerService_service
                        larsim_MCCheater_BackTrackerService_service
                        larsim_MCCheater_ParticleInventoryService_service
                        nug4_ParticleNavigation 
                        nusimdata_SimulationBase
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        messagefacility::MF_MessageLogger
                        cetlib cetlib_except
                        ROOT_BASIC_LIB_LIST
	  #ProtoDUNEDataUtils
)

install_headers()
install_fhicl()
install_source()
install_scripts()
