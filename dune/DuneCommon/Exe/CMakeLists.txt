# dunetpc/dune/DuneCommon/Exe
#
# David Adams
# May 2018
#
# Instructions to build and install duneRunData

cet_make_exec(duneRunData
  SOURCE duneRunData.cxx
  LIBRARIES
    dune_ArtSupport
    fhiclcpp
    cetlib cetlib_except
)

install_scripts(
  LIST duneHelp
)
