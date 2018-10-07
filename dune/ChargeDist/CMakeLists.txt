
simple_plugin(ChargeDist "module"
                        larcore_Geometry_Geometry_service
                        larcorealg_Geometry
                        lardataobj_RawData
                        lardataobj_RecoBase
                        larsim_Simulation 
                        nutools_ParticleNavigation 
                        lardataobj_Simulation                       
                        nusimdata_SimulationBase
                        lardata_Utilities
                        
                        larsim_MCCheater_BackTrackerService_service
						larsim_MCCheater_ParticleInventoryService_service
						
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
install_scripts()
