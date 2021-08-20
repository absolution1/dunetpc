cet_build_plugin(RawData311InputDriver art::service
			EventDecoder_service
			dlardaq_service

			lardataobj_RawData
			lardata_Utilities
			larcorealg_Geometry
			ART_FRAMEWORK_CORE
			ART_FRAMEWORK_IO_SOURCES
			ART_FRAMEWORK_PRINCIPAL
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			messagefacility::MF_MessageLogger
			
			CETLIB
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)

cet_build_plugin(dlardaq art::service
			lardataobj_RawData
			lardataobj_RecoBase
			larreco_RecoAlg
			larreco_Calorimetry
			larreco_RecoAlg_PMAlg
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
			
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
			
			CETLIB
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)


cet_build_plugin(EventDecoder art::service
			dlardaq_service
			HuffDataCompressor_service
			
			lardataobj_RawData
			lardataobj_RecoBase
			larreco_RecoAlg
			larreco_Calorimetry
			larreco_RecoAlg_PMAlg
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
			
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
			
			CETLIB
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)


cet_build_plugin(HuffDataCompressor art::service
			dlardaq_service

			lardataobj_RawData
			lardataobj_RecoBase
			larreco_RecoAlg
			larreco_Calorimetry
			larreco_RecoAlg_PMAlg
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
			
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
			
			CETLIB
			ROOT_BASIC_LIB_LIST
			BASENAME_ONLY
)




install_headers()
install_fhicl()
install_source()
install_scripts()
