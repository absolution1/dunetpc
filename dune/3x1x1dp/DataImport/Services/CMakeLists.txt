simple_plugin(RawData311InputDriver "service"
			EventDecoder_service
			dlardaq_service

			lardataobj_RawData
			lardata_Utilities
			larcorealg_Geometry
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_SOURCES}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_PERSISTENCY_COMMON}
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
			BASENAME_ONLY
)

simple_plugin(dlardaq "service"
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
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
			BASENAME_ONLY
)


simple_plugin(EventDecoder "service"
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
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
			BASENAME_ONLY
)


simple_plugin(HuffDataCompressor "service"
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
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
			BASENAME_ONLY
)




install_headers()
install_fhicl()
install_source()
install_scripts()
