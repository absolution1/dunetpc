# basic source code CMakeLists.txt

include_directories( $ENV{IFDHC_FQ_DIR}/inc )

art_make(
	LIBRARY_NAME OverlayRawDataDUNE35t_module
	EXCLUDE RawDigitOverlayDUNE35tAna_module.cc
	MODULE_LIBRARIES
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
                        art_Persistency_Common canvas
                        art_Utilities canvas
			ifdh
                        ${MF_MESSAGELOGGER}
                        ${BOOST_LIB}
                        ${FHICLCPP}
                        cetlib cetlib_except
												${ROOT_TREEPLAYER}
												${ROOT_BASIC_LIB_LIST}
)							   

art_make(
	LIBRARY_NAME RawDigitOverlayDUNE35tAna_module
	EXCLUDE OverlayRawDataDUNE35t_module.cc
	MODULE_LIBRARIES
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
			${BOOST_LIB}
                        ${FHICLCPP}
                        cetlib cetlib_except
)							   

install_headers()
install_source()
install_fhicl()

# Libraries I didn't need
#			dune_RawDigitAdderAna
#
#			${ART_FRAMEWORK_IO_SOURCES}
#                        art_Persistency_Provenance canvas
#                        
#                        ${ROOT_BASIC_LIB_LIST}
