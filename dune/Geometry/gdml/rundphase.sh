#!/bin/bash

####################
### Far detector ###
####################

echo "--- Generating Far Detector Geometry ---"

# with wires for larcore
perl generate_dunedphase10kt_v2.pl -w 1 -o dphase_geo.xml
perl make_gdml.pl -i dphase_geo.xml -o dunedphase10kt_v2.gdml

# clean up
rm dunedphase10kt_v2_*.gdml
rm dphase_geo.xml

# no wires for geant
perl generate_dunedphase10kt_v2.pl -w 0 -o dphase_geo_nowires.xml
perl make_gdml.pl -i dphase_geo_nowires.xml -o dunedphase10kt_v2_nowires.gdml

# clean up
rm dunedphase10kt_v2_nowires_*.gdml
rm dphase_geo_nowires.xml


#################################
### Far detector (drift in Y) ###
#################################

echo "--- Generating Far Detector Geometry Drift in Y---"

# with wires for larcore
perl generate_dunedphase10kt_v2_driftY.pl -w 1 -o dphase_geo_driftY.xml --pmtdensity 1 --TPBFoils 1
perl make_gdml.pl -i dphase_geo_driftY.xml -o dunedphase10kt_v2_driftY.gdml

# clean up
rm dunedphase10kt_v2_driftY_*.gdml
rm dphase_geo_driftY.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY.pl -w 0 -o dphase_geo_driftY_nowires.xml --pmtdensity 1 --TPBFoils 1
perl make_gdml.pl -i dphase_geo_driftY_nowires.xml -o dunedphase10kt_v2_driftY_nowires.gdml

# clean up
rm dunedphase10kt_v2_driftY_nowires_*.gdml
rm dphase_geo_driftY_nowires.xml

# with wires for larcore
perl generate_dunedphase10kt_v2_driftY.pl -w 1 -o dphase_geo_driftY_HalfFoil.xml --pmtdensity 1 --TPBFoils 2
perl make_gdml.pl -i dphase_geo_driftY_HalfFoil.xml -o dunedphase10kt_v2_driftY_HalfFoil.gdml

# clean up
rm dunedphase10kt_v2_driftY_HalfFoil_*.gdml
rm dphase_geo_driftY_HalfFoil.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY.pl -w 0 -o dphase_geo_driftY_HalfFoil_nowires.xml --pmtdensity 1 --TPBFoils 2
perl make_gdml.pl -i dphase_geo_driftY_HalfFoil_nowires.xml -o dunedphase10kt_v2_driftY_HalfFoil_nowires.gdml

# clean up
rm dunedphase10kt_v2_driftY_HalfFoil_nowires_*.gdml
rm dphase_geo_driftY_HalfFoil_nowires.xml

# with wires for larcore
perl generate_dunedphase10kt_v2_driftY.pl -w 1 -o dphase_geo_driftY_2PMTperM2.xml --pmtdensity 2 --TPBFoils 0
perl make_gdml.pl -i dphase_geo_driftY_2PMTperM2.xml -o dunedphase10kt_v2_driftY_2PMTperM2.gdml

# clean up
rm dunedphase10kt_v2_driftY_2PMTperM2_*.gdml
rm dphase_geo_driftY_2PMTperM2.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY.pl -w 0 -o dphase_geo_driftY_2PMTperM2_nowires.xml --pmtdensity 2 --TPBFoils 0
perl make_gdml.pl -i dphase_geo_driftY_2PMTperM2_nowires.xml -o dunedphase10kt_v2_driftY_2PMTperM2_nowires.gdml

# clean up
rm dunedphase10kt_v2_driftY_2PMTperM2_nowires_*.gdml
rm dphase_geo_driftY_2PMTperM2_nowires.xml


########################################
### Far detector workspace: 2x1 TPCs ###
########################################

echo "--- Generating Far Detector Geometry 2x1 workspace ---"

# with wires for larcore
perl generate_dunedphase10kt_v2.pl -w 1 -k 1 -o dphase_geo_workspace.xml
perl make_gdml.pl -i dphase_geo_workspace.xml -o dunedphase10kt_v2_workspace.gdml

# clean up
rm dunedphase10kt_v2_workspace_*.gdml
rm dphase_geo_workspace.xml

# no wires for geant
perl generate_dunedphase10kt_v2.pl -w 0 -k 1 -o dphase_geo_workspace_nowires.xml
perl make_gdml.pl -i dphase_geo_workspace_nowires.xml -o dunedphase10kt_v2_workspace_nowires.gdml

# clean up
rm dunedphase10kt_v2_workspace_nowires_*.gdml
rm dphase_geo_workspace_nowires.xml


#####################################################
### Far detector workspace: 2x1 TPCs (drift in Y) ###
#####################################################

echo "--- Generating Far Detector Geometry 2x1 workspace Drift in Y ---"

# with wires for larcore
perl generate_dunedphase10kt_v2_driftY.pl -w 1 -k 1 -o dphase_geo_driftY_workspace.xml --pmtdensity 1 --TPBFoils 0
perl make_gdml.pl -i dphase_geo_driftY_workspace.xml -o dunedphase10kt_v2_driftY_workspace.gdml

# clean up
rm dunedphase10kt_v2_driftY_workspace_*.gdml
rm dphase_geo_driftY_workspace.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY.pl -w 0 -k 1 -o dphase_geo_driftY_workspace_nowires.xml --pmtdensity 1 --TPBFoils 0
perl make_gdml.pl -i dphase_geo_driftY_workspace_nowires.xml -o dunedphase10kt_v2_driftY_workspace_nowires.gdml

# clean up
rm dunedphase10kt_v2_driftY_workspace_nowires_*.gdml
rm dphase_geo_driftY_workspace_nowires.xml

########################################
### Far detector workspace: 4x2 TPCs ###
########################################

echo "--- Generating Far Detector Geometry 4x2 workspace ---"

# with wires for larcore
perl generate_dunedphase10kt_v2.pl -w 1 -k 2 -o dphase_geo_workspace4x2.xml
perl make_gdml.pl -i dphase_geo_workspace4x2.xml -o dunedphase10kt_v2_workspace4x2.gdml

# clean up
rm dunedphase10kt_v2_workspace4x2_*.gdml
rm dphase_geo_workspace4x2.xml

# no wires for geant
perl generate_dunedphase10kt_v2.pl -w 0 -k 2 -o dphase_geo_workspace4x2_nowires.xml -w 0
perl make_gdml.pl -i dphase_geo_workspace4x2_nowires.xml -o dunedphase10kt_v2_workspace4x2_nowires.gdml

# clean up
rm dunedphase10kt_v2_workspace4x2_nowires_*.gdml
rm dphase_geo_workspace4x2_nowires.xml


#####################################################
### Far detector workspace: 4x2 TPCs (drift in Y) ###
#####################################################

echo "--- Generating Far Detector Geometry 4x2 workspace Drift in Y ---"

# with wires for larcore
perl generate_dunedphase10kt_v2_driftY.pl -w 1 -k 2 -o dphase_geo_driftY_workspace4x2.xml --pmtdensity 1 --TPBFoils 0
perl make_gdml.pl -i dphase_geo_driftY_workspace4x2.xml -o dunedphase10kt_v2_driftY_workspace4x2.gdml

# clean up
rm dunedphase10kt_v2_driftY_workspace4x2_*.gdml
rm dphase_geo_driftY_workspace4x2.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY.pl -w 0 -k 2 -o dphase_geo_driftY_workspace4x2_nowires.xml -w 0 --pmtdensity 1 --TPBFoils 0
perl make_gdml.pl -i dphase_geo_driftY_workspace4x2_nowires.xml -o dunedphase10kt_v2_driftY_workspace4x2_nowires.gdml

# clean up
rm dunedphase10kt_v2_driftY_workspace4x2_nowires_*.gdml
rm dphase_geo_driftY_workspace4x2_nowires.xml


#####################################################
### Far detector workspace: 3x3 TPCs (drift in Y) ###
#####################################################

echo "--- Generating Far Detector Geometry 3x3 workspace Drift in Y ---"

# with wires for larcore
perl generate_dunedphase10kt_v2_driftY.pl -w 1 -k 3 -o dphase_geo_driftY_workspace3x3.xml --pmtdensity 1 --TPBFoils 0
perl make_gdml.pl -i dphase_geo_driftY_workspace3x3.xml -o dunedphase10kt_v2_driftY_workspace3x3.gdml

# clean up
rm dunedphase10kt_v2_driftY_workspace3x3_*.gdml
rm dphase_geo_driftY_workspace3x3.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY.pl -w 0 -k 3 -o dphase_geo_driftY_workspace3x3_nowires.xml -w 0 --pmtdensity 1 --TPBFoils 0
perl make_gdml.pl -i dphase_geo_driftY_workspace3x3_nowires.xml -o dunedphase10kt_v2_driftY_workspace3x3_nowires.gdml

# clean up
rm dunedphase10kt_v2_driftY_workspace3x3_nowires_*.gdml
rm dphase_geo_driftY_workspace3x3_nowires.xml


####################
### protoDUNE-DP ###
####################

echo "--- Generating protoDUNEDP Geometry ---"

# with wires for larcore
perl generate_protodunedphase.pl -w 1 -opt 0 -o protodunedphase_geo.xml
perl make_gdml.pl -i protodunedphase_geo.xml -o protodunedphase.gdml

# clean up
rm protodunedphase_geo.xml
rm protodunedphase_*.gdml

# no wires for geant
perl generate_protodunedphase.pl -w 0 -opt 0 -o protodunedphase_geo_nowires.xml
perl make_gdml.pl -i protodunedphase_geo_nowires.xml -o protodunedphase_nowires.gdml

# clean up
rm protodunedphase_geo_nowires.xml
rm protodunedphase_nowires_*.gdml


#########################################
### protoDUNE-DP rotated (drift in Y) ###
#########################################

echo "--- Generating protoDUNEDP Geometry Drift in Y ---"

# with wires for larcore
perl generate_protodunedphase_driftY.pl -w 1 -opt 0 -o protodunedphase_driftY_geo.xml
perl make_gdml.pl -i protodunedphase_driftY_geo.xml -o protodunedphase_driftY.gdml

# clean up
rm protodunedphase_driftY_geo.xml
rm protodunedphase_driftY_*.gdml

# no wires for geant
perl generate_protodunedphase_driftY.pl -w 0 -opt 0 -o protodunedphase_driftY_nowires_geo.xml
perl make_gdml.pl -i protodunedphase_driftY_nowires_geo.xml -o protodunedphase_driftY_nowires.gdml

# clean up
rm protodunedphase_driftY_nowires_geo.xml
rm protodunedphase_driftY_nowires_*.gdml


##################################################
### protoDUNE-DP for photon library generation ###
##################################################

echo "--- Generating protoDUNEDP OPTICAL Geometry ---"

# all detail to generate photon libraries
perl generate_protodunedphase.pl -w 1 -opt 1 -o protodunedphase_optical_geo.xml
perl make_gdml.pl -i protodunedphase_optical_geo.xml -o protodunedphase_optical.gdml

# clean up
rm protodunedphase_optical_geo.xml
rm protodunedphase_optical_*.gdml

# no wires for geant
perl generate_protodunedphase.pl -w 0 -opt 1 -o protodunedphase_optical_nowires_geo.xml
perl make_gdml.pl -i protodunedphase_optical_nowires_geo.xml -o protodunedphase_optical_nowires.gdml

# clean up
rm protodunedphase_optical_nowires_geo.xml
rm protodunedphase_optical_nowires_*.gdml

echo "--- Generating protoDUNEDP OPTICAL Geometry Drift in Y ---"

# all detail to generate photon libraries
perl generate_protodunedphase_driftY.pl -w 1 -opt 1 -o protodunedphase_driftY_optical_geo.xml
perl make_gdml.pl -i protodunedphase_driftY_optical_geo.xml -o protodunedphase_driftY_optical.gdml

# clean up
rm protodunedphase_driftY_optical_geo.xml
rm protodunedphase_driftY_optical_*.gdml

# no wires for geant
perl generate_protodunedphase_driftY.pl -w 0 -opt 1 -o protodunedphase_driftY_optical_nowires_geo.xml
perl make_gdml.pl -i protodunedphase_driftY_optical_nowires_geo.xml -o protodunedphase_driftY_optical_nowires.gdml

# clean up
rm protodunedphase_driftY_optical_nowires_geo.xml
rm protodunedphase_driftY_optical_nowires_*.gdml

########################
### 3x1x1 dual phase ###
########################

echo "--- Generating 3x1x1 Geometry ---"

# with wires for larcore
perl generate_3x1x1dphase.pl -w 1 -o 3x1x1dphase_geo.xml
perl make_gdml.pl -i 3x1x1dphase_geo.xml -o 3x1x1dphase.gdml

# clean up
rm 3x1x1dphase_geo.xml
rm 3x1x1dphase_*.gdml

# no wires for geant
perl generate_3x1x1dphase.pl -w 0 -o 3x1x1dphase_geo_nowires.xml
perl make_gdml.pl -i 3x1x1dphase_geo_nowires.xml -o 3x1x1dphase_nowires.gdml

# clean up
rm 3x1x1dphase_geo_nowires.xml
rm 3x1x1dphase_nowires_*.gdml



#####################################
### 3x1x1 dual phase (drift in Y) ###
#####################################

echo "--- Generating 3x1x1 Geometry Drift in Y ---"

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

