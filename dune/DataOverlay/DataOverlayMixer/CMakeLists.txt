include_directories($ENV{IFDHC_FQ_DIR}/inc)

simple_plugin(OverlayRawDataDUNE35t "module"
              dune_RawDigitMixer
              dune_OpDetWaveformMixer
              lardataobj_RawData
              nusimdata_SimulationBase
              lardataobj_MCBase
              larevt_CalibrationDBI_IOVData
              larevt_CalibrationDBI_Providers
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_IO_PRODUCTMIX}
              ${ART_FRAMEWORK_IO_ROOT}
              ${ART_FRAMEWORK_MODULES}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
              ${ART_PERSISTENCY_PROVENANCE}
              art_Framework_IO_Root_detail_sources
              art_Persistency_Common
              art_Utilities canvas
              ifdh
              ${MF_MESSAGELOGGER}
              ${BOOST_LIB}
              ${FHICLCPP}
              cetlib cetlib_except
              ${ROOT_TREEPLAYER}
              ${ROOT_BASIC_LIB_LIST})

simple_plugin(RawDigitOverlayDUNE35tAna "module"
              larcore_Geometry_Geometry_service
              larcorealg_Geometry
              dune_RawDigitAdderAna
              lardataobj_RawData
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_IO_ROOT}
              ${ART_FRAMEWORK_MODULES}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              ${ART_PERSISTENCY_PROVENANCE}
              ${BOOST_LIB}
              ${FHICLCPP}
              cetlib
              cetlib_except)

install_headers()
install_source()
install_fhicl()
