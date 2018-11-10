# get libxml2
cet_find_library( XML2 NAMES xml2 PATHS $ENV{LIBXML2_FQ_DIR}/lib NO_DEFAULT_PATH )
# the GENIE libraries are already defined, but we need the lib list
set(GENIE_LIB_LIST  ${GALGORITHM}
                    ${GBARYONRESONANCE}
                    ${GBASE}
                    ${GBODEKYANG}
                    ${GCHARM}
                    ${GCOH}
                    ${GDFRC}
                    ${GDIS}
                    ${GCROSSSECTIONS}
                    ${GDECAY}
                    ${GELAS}
                    ${GELFF}
                    ${GHEP}
                    ${GEVGCORE}
                    ${GEVGMODULES}
                    ${GEVGDRIVERS}
                    ${GGIBUU}
                    ${GHADRONTRANSP}
                    ${GFRAGMENTATION}
                    ${GINTERACTION}
                    ${GLLEWELLYNSMITH}
                    ${GMEC}
                    ${GMESSENGER}
                    ${GNUGAMMA}
                    ${GNUE}
                    ${GNTUPLE}
                    ${GNUCLEAR}
                    ${GNUMERICAL}
                    ${GQPM}
                    ${GPDG}
                    ${GPDF}
                    ${GQEL}
                    ${GRES}
                    ${GREGISTRY}
                    ${GREINSEHGAL}
                    ${GUTILS}
                    ${GGEO}
                    ${GFLUXDRIVERS}
                    ${GMUELOSS}
                    ${GREWEIGHT}
                    ${GNUCLEONDECAY} )


include_directories ( $ENV{GENIE_INC}/GENIE )

SET(GENIE_EV_LIBS
  ${GENIE_LIB_LIST}
  ${LOG4CPP}
  ${XML2}
  ${LHAPDF}
  ${ROOT_BASIC_LIB_LIST}
  ${ROOT_EG}
  ${ROOT_TREEPLAYER}
  ${ROOT_FFTW}
  ${ROOT_REFLEX}
  ${ROOTSYS}/lib/libEGPythia6.so
  ${ROOT_GUI}
  ${PYTHIA6}
)

# prior to GENIE R-2_9_0 Sehgal's name was spelled incorrectly
# so if minor version # is 8 or less use GReinSeghal, otherwise GReinSehgal
# parse_ups_version returns VMAJ, VMIN, and VPRJ
parse_ups_version( ${GENIE_VERSION} )
set ( GRSSPELLING GReinSehgal )
if ( ${VMAJ} EQUAL 2 AND ${VMIN} LESS 9 )
  set ( GRSSPELLING GReinSeghal )
endif ( ${VMAJ} EQUAL 2 AND ${VMIN} LESS 9 )
message(STATUS "GENIE version ${GENIE_VERSION} major ${VMAJ} minor ${VMIN} ${GRSSPELLING}")
cet_find_library( GREINSEHGAL NAMES ${GRSSPELLING} PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
# nusystematics
cet_find_library( SYSTTOOLS NAMES systematicstools_interface systematicstools_interpreters systematicstools_utility PATHS ENV SYSTEMATICSTOOLS_LIB NO_DEFAULT_PATH )

SET(SYSTTOOLS_LIBS
  systematicstools_interface
  systematicstools_interpreters
  systematicstools_utility)

simple_plugin(CAFMaker   "module"
                        nutools_NuReweight
                        nutools_NuReweight_art
			nusimdata_SimulationBase
                        ${SYSTTOOLS_LIBS}
			${GENIE_EV_LIBS}
			${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities 
                        canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        cetlib cetlib_except
                        ${ROOT_BASIC_LIB_LIST}
                        CVNFunc
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
