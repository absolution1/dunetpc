# dunetpc/dune/ArtSupport/exe
#
# David Adams
# August 2016
#
# Instructions to build and install fcldump.

cet_make_exec(fcldump
  SOURCE fcldump.cxx
  LIBRARIES
   fhiclcpp
    cetlib
)

install_source()
