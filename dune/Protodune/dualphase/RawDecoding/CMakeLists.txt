simple_plugin(PDDPRawInput "source"
			PDDPRawInputDriver_service
                        lardataobj_RawData
                        lardata_Utilities
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_PERSISTENCY_COMMON}
                        art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${FHICLCPP}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        BASENAME_ONLY
)

simple_plugin(PDDPRawInputDriver "service"
			lardataobj_RawData
                        lardata_Utilities
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_PERSISTENCY_COMMON}
                        art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			pthread
			${FHICLCPP}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
                        BASENAME_ONLY
)



install_headers()
install_fhicl()
install_source()
install_scripts()

