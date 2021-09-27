
cet_build_plugin(ShowerAnalysis art::module
			     lardataobj_RawData
		      	     lardataobj_RecoBase
		      	     lardata_ArtDataHelper
		             larcorealg_Geometry
			     larcore_Geometry_Geometry_service
			     lardata_Utilities
			     larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
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
			     
			     cetlib cetlib_except
			     ROOT_BASIC_LIB_LIST
)


install_headers()
install_fhicl()
install_source()
install_scripts()


# art_make(LIBRARY_NAME  ShowerAna
# 	 LIB_LIBRARIES lardataobj_RawData
# 	 	       lardataobj_RecoBase	  
# 	 	       larcorealg_Geometry  
# 	 	       larcore_Geometry_Geometry_service
# 	 	       lardata_Utilities
# 		       larsim_MCCheater_BackTrackerService_service
#           larsim_MCCheater_ParticleInventoryService_service
# 	 	       nusimdata_SimulationBase
#          	       ${ART_FRAMEWORK_CORE}
# 	 	       ${ART_FRAMEWORK_PRINCIPAL}
# 	 	       ${ART_FRAMEWORK_SERVICES_REGISTRY}
#          	       ${ART_ROOT_IO_TFILE_SUPPORT} ROOT::Core
#          	       ${ART_ROOT_IO_TFILESERVICE_SERVICE}
#  	 	       art_Persistency_Common canvas
# 	 	       art_Persistency_Provenance canvas
# 	 	       art_Utilities canvas
# 	 	       ${messagefacility::MF_MessageLogger}
# 	 	       
# 	 	       cetlib cetlib_except
# 	 	       ROOT_BASIC_LIB_LIST
# 		       ShowerAna
# )

# install_headers()
# install_fhicl()
# install_source()
# install_scripts()
