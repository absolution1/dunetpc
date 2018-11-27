#!/bin/bash

#####################################
### 3x1x1 dual phase (drift in Y) ###
#####################################

# with wires for larcore
perl generate_3x1x1dphase_driftY.pl -w 1 -o 3x1x1dphase_driftY_geo.xml
perl make_gdml.pl -i 3x1x1dphase_driftY_geo.xml -o 3x1x1dphase_driftY.gdml

# clean up
rm 3x1x1dphase_driftY_geo.xml
rm 3x1x1dphase_driftY_*.gdml

# no wires for geant
perl generate_3x1x1dphase_driftY.pl -w 0 -o 3x1x1dphase_driftY_geo_nowires.xml
perl make_gdml.pl -i 3x1x1dphase_driftY_geo_nowires.xml -o 3x1x1dphase_driftY_nowires.gdml

# clean up
rm 3x1x1dphase_driftY_geo_nowires.xml
rm 3x1x1dphase_driftY_nowires_*.gdml

