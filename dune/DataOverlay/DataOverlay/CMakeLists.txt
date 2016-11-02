# basic source code CMakeLists.txt

art_make_library( LIBRARY_NAME  dune_RawDigitAdderAna
	SOURCE RawDigitAdderAna.cxx
	LIBRARIES
                 lardataobj_RawData
                 ${ROOT_BASIC_LIB_LIST}
)

art_make_library( LIBRARY_NAME  dune_RawDigitAdder
	SOURCE RawDigitAdder.cxx
	             RawDigitAdder_HardSaturate.cxx
		     RawDigitAdder_35t.cxx
	LIBRARIES
                 ${ROOT_BASIC_LIB_LIST}
		 larsim_DetSim
)

art_make_library( LIBRARY_NAME  dune_RawDigitMixer
	SOURCE RawDigitMixer.cxx
	LIBRARIES
		dune_RawDigitAdder
		lardataobj_RawData
		${ROOT_BASIC_LIB_LIST}
)

art_make_library( LIBRARY_NAME  dune_OpDetWaveformMixer
	SOURCE OpDetWaveformMixer.cxx
	LIBRARIES
		dune_RawDigitAdder
		lardataobj_RawData
		${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_source()
