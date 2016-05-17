
simple_plugin(ShowerAnalysis "module"
			     lardata_RawData
		      	     lardata_RecoBase
		      	     lardata_RecoBaseArt
		             larcore_Geometry
			     larcore_Geometry_Geometry_service
			     lardata_Utilities
			     larsim_MCCheater_BackTracker_service
			     ${SIMULATIONBASE}
                             ${ART_FRAMEWORK_CORE}
			     ${ART_FRAMEWORK_PRINCIPAL}
			     ${ART_FRAMEWORK_SERVICES_REGISTRY}
                             ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                             ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			     ${ART_PERSISTENCY_COMMON}
			     ${ART_PERSISTENCY_PROVENANCE}
			     ${ART_UTILITIES}
			     ${MF_MESSAGELOGGER}
			     ${MF_UTILITIES}
			     ${CETLIB}
			     ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()
install_scripts()


# art_make(LIBRARY_NAME  ShowerAna
# 	 LIB_LIBRARIES lardata_RawData
# 	 	       lardata_RecoBase	  
# 	 	       larcore_Geometry  
# 	 	       larcore_Geometry_Geometry_service
# 	 	       lardata_Utilities
# 		       larsim_MCCheater_BackTracker_service
# 	 	       ${SIMULATIONBASE}
#          	       ${ART_FRAMEWORK_CORE}
# 	 	       ${ART_FRAMEWORK_PRINCIPAL}
# 	 	       ${ART_FRAMEWORK_SERVICES_REGISTRY}
#          	       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
#          	       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
#  	 	       ${ART_PERSISTENCY_COMMON}
# 	 	       ${ART_PERSISTENCY_PROVENANCE}
# 	 	       ${ART_UTILITIES}
# 	 	       ${MF_MESSAGELOGGER}
# 	 	       ${MF_UTILITIES}
# 	 	       ${CETLIB}
# 	 	       ${ROOT_BASIC_LIB_LIST}
# 		       ShowerAna
# )

# install_headers()
# install_fhicl()
# install_source()
# install_scripts()
