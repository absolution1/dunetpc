# dunetpc/dune/DuneInterface/Data/CMakeLists.txt
#
# David Adams
# September 2017
#
# CMake instructions for DUNE service interface classes.

add_subdirectory(test)

art_make(BASENAME_ONLY
         NO_PLUGINS
        )

# Including DataMap.h confers a transitive dependency on ROOT::Hist
cet_make_library(LIBRARY_NAME dune_DuneInterface_Data_DataMap INTERFACE
  SOURCE DataMap.h
  LIBRARIES INTERFACE ROOT::Hist)

# Including RunData.h confers a transitive dependency on ROOT::Hist
cet_make_library(LIBRARY_NAME dune_DuneInterface_Data_RunData INTERFACE
  SOURCE RunData.h
  LIBRARIES INTERFACE ROOT::Hist)

install_headers()
install_source()
