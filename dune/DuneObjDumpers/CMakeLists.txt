art_make(
  LIB_LIBRARIES
    ROOT_EG
    ROOT::Core
  MODULE_LIBRARIES
  #dune_DuneObjDumpers
    dune_DuneObj
    nusimdata_SimulationBase
    ART_FRAMEWORK_CORE
    ART_FRAMEWORK_PRINCIPAL
    art_Persistency_Common
    art_Persistency_Provenance
    art_Utilities
    canvas
    messagefacility::MF_MessageLogger

  cetlib cetlib_except
    ROOT_EG
    ROOT_BASIC_LIB_LIST
  )

install_headers()
install_fhicl()
install_source()
