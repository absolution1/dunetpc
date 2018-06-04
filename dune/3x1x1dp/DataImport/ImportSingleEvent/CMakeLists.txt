simple_plugin(ImportSingle311Event "module"
			EventDecoder_service

			lardataobj_RawData
			lardataobj_RecoBase
			larreco_RecoAlg
			larreco_Calorimetry
			larreco_RecoAlg_PMAlg
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			
			nusimdata_SimulationBase
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
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
