#!/bin/bash

# with wires for larcore
perl generate_dunedphase10kt_v2.pl -o dphase_geo.xml
make_gdml.pl -i dphase_geo.xml -o dunedphase10kt_v2.gdml

# clean up
#rm dunedphase10kt_v2_*.gdml

# no wires for geant
perl generate_dunedphase10kt_v2.pl -o dphase_geo.xml -w 0
make_gdml.pl -i dphase_geo.xml -o dunedphase10kt_v2_nowires.gdml

# clean up
#rm dunedphase10kt_v2_nowires_*.gdml
#rm dphase_geo.xml


#root -l dunedphase10kt_geo.C
