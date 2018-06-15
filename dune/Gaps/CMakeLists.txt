simple_plugin(GapFilter "module"
                        lardataobj_RawData
                        lardataobj_RecoBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        ${FHICLCPP}
                        cetlib cetlib_except
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
                        nusimdata_SimulationBase
              BASENAME_ONLY
)

simple_plugin(GapWidth "module"
                       lardataobj_RawData
                       lardataobj_RecoBase
                       lardata_ArtDataHelper
                       larcorealg_Geometry
                       larcore_Geometry_Geometry_service
                       lardata_Utilities
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
                       
                       cetlib cetlib_except
                       ${ROOT_BASIC_LIB_LIST}
                       ${ROOT_MINUIT}
                       ${ROOT_MINUIT2}
             BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
