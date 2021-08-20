include_directories ( $ENV{GENIE_INC}/GENIE )

# nusystematics
##cet_find_library( SYSTTOOLS NAMES systematicstools_interface systematicstools_interpreters systematicstools_utility PATHS ENV SYSTEMATICSTOOLS_LIB NO_DEFAULT_PATH )

cet_build_plugin(CAFMaker   art::module
                        nugen_NuReweight
                        nugen_NuReweight_art
			nusimdata_SimulationBase
                        GFWGHEP
                        GFWEG
			GFWPARDAT
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
			GFWALG
			GFWPARDAT
			GPHDISXS
			GPHDISXS
			GPHCHMXS
			GPHCOHEG
      GPHCOHXS
			GPHDFRCEG
			GPHDISEG
			GPHXSIG
			GPHDCY
			GPHQELXS
			GPHQELXS
			GFWGHEP
			GFWEG
			GPHCMN
			GFWEG
			GPHRESXS
			GPHHADTRANSP
			GPHHADNZ
			GFWINT
			GPHQELXS
			GPHMNUCEG
      GPHMNUCXS
			GPHAMNGEG
      GPHAMNGXS
			GPHNUELEG
      GPHNUELXS
			GFWNTP
			GPHNUCLST
			GFWNUM
			GFWPARDAT
			GPHPDF
			GPHQELEG
			GPHRESEG
      GPHRESXS
			GFWREG
			GPHRESXS
			GFWUTL
			GTLGEO
			GTLFLX
			GPHMEL
			GPHNDCY
			GFWMSG
                        systematicstools_interface
			systematicstools_interpreters
			systematicstools_utility
                        dune_CVN_func
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
