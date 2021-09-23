cet_build_plugin(PDDPChannelMap art::service
                        ART_FRAMEWORK_CORE
			ART_FRAMEWORK_IO_SOURCES
                        ART_FRAMEWORK_PRINCIPAL
                        art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			FHICLCPP
                        messagefacility::MF_MessageLogger
                        CETLIB
                        ROOT_BASIC_LIB_LIST
                        BASENAME_ONLY
)

cet_build_plugin(PDDPRawInput art::source
			PDDPRawInputDriver_service
                        lardataobj_RawData
                        lardata_Utilities
                        ART_FRAMEWORK_CORE
			ART_FRAMEWORK_IO_SOURCES
                        ART_FRAMEWORK_PRINCIPAL
                        art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			FHICLCPP
                        messagefacility::MF_MessageLogger
                        CETLIB
                        ROOT_BASIC_LIB_LIST
                        BASENAME_ONLY
)

cet_build_plugin(PDDPRawInputDriver art::service
			PDDPChannelMap_service
			pthread
			lardataobj_RawData
                        lardata_Utilities
                        ART_FRAMEWORK_CORE
			ART_FRAMEWORK_IO_SOURCES
                        ART_FRAMEWORK_PRINCIPAL
                        art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			FHICLCPP
                        messagefacility::MF_MessageLogger
                        CETLIB
                        ROOT_BASIC_LIB_LIST
                        BASENAME_ONLY
)

cet_build_plugin(PDDPTPCRawDecoderTest art::module
			lardataobj_RawData
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
              		BASENAME_ONLY
)


install_headers()
install_fhicl()
install_source()
install_scripts()

