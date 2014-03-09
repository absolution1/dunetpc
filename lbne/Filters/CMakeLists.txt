include_directories ( $ENV{GEANT4_FQ_DIR}/include )

art_make( BASENAME_ONLY
           MODULE_LIBRARIES Filters
                        Simulation
                        ${SIMULATIONBASE}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
                        ${G4_LIB_LIST}
		)

#install_headers()
install_fhicl()
install_source()

