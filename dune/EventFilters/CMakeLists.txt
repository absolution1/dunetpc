
cet_build_plugin(DuneEventFilter art::module
			dune_ArtSupport
			lardataobj_RawData
                        lardataobj_RecoBase
	                larcorealg_Geometry
	                larcore_Geometry_Geometry_service
                        lardata_Utilities    
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
	                FHICLCPP
			cetlib cetlib_except
                        ROOT_GEOM
                        ROOT_XMLIO
                        ROOT_GDML
			ROOT_BASIC_LIB_LIST
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
