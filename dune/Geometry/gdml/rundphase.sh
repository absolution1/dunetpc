#!/bin/bash

####################
### Far detector ###
####################

echo "--- Generating Far Detector Geometry ---"

# with wires for larcore
perl generate_dunedphase10kt_v2_legacy.pl -w 1 -o dphase_geo.xml
perl make_legacy.gdml.pl -i dphase_geo.xml -o dunedphase10kt_v2_legacy.gdml

# clean up
rm dunedphase10kt_v2_*_legacy.gdml
rm dphase_geo.xml

# no wires for geant
perl generate_dunedphase10kt_v2_legacy.pl -w 0 -o dphase_geo_nowires.xml
perl make_legacy.gdml.pl -i dphase_geo_nowires.xml -o dunedphase10kt_v2_nowires_legacy.gdml

# clean up
rm dunedphase10kt_v2_nowires_*_legacy.gdml
rm dphase_geo_nowires.xml


#################################
### Far detector (drift in Y) ###
#################################

echo "--- Generating Far Detector Geometry Drift in Y--- No Foils"

# with wires for larcore
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 1 -o dphase_geo_driftY.xml --pmtdensity 1 --TPBFoils 0
perl make_legacy.gdml.pl -i dphase_geo_driftY.xml -o dunedphase10kt_v2_driftY_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_*_legacy.gdml
rm dphase_geo_driftY.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 0 -o dphase_geo_driftY_nowires.xml --pmtdensity 1 --TPBFoils 0
perl make_legacy.gdml.pl -i dphase_geo_driftY_nowires.xml -o dunedphase10kt_v2_driftY_nowires_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_nowires_*_legacy.gdml
rm dphase_geo_driftY_nowires.xml


echo "--- Generating Far Detector Geometry Drift in Y--- Full Foil"

# with wires for larcore
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 1 -o dphase_geo_driftY_Foil.xml --pmtdensity 1 --TPBFoils 1
perl make_legacy.gdml.pl -i dphase_geo_driftY_Foil.xml -o dunedphase10kt_v2_driftY_Foil_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_Foil_*_legacy.gdml
rm dphase_geo_driftY_Foil.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 0 -o dphase_geo_driftY_Foil_nowires.xml --pmtdensity 1 --TPBFoils 1
perl make_legacy.gdml.pl -i dphase_geo_driftY_Foil_nowires.xml -o dunedphase10kt_v2_driftY_Foil_nowires_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_Foil_nowires_*_legacy.gdml
rm dphase_geo_driftY_Foil_nowires.xml

echo "--- Generating Far Detector Geometry Drift in Y--- Half Foil"
# with wires for larcore
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 1 -o dphase_geo_driftY_HalfFoil.xml --pmtdensity 1 --TPBFoils 2
perl make_legacy.gdml.pl -i dphase_geo_driftY_HalfFoil.xml -o dunedphase10kt_v2_driftY_HalfFoil_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_HalfFoil_*_legacy.gdml
rm dphase_geo_driftY_HalfFoil.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 0 -o dphase_geo_driftY_HalfFoil_nowires.xml --pmtdensity 1 --TPBFoils 2
perl make_legacy.gdml.pl -i dphase_geo_driftY_HalfFoil_nowires.xml -o dunedphase10kt_v2_driftY_HalfFoil_nowires_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_HalfFoil_nowires_*_legacy.gdml
rm dphase_geo_driftY_HalfFoil_nowires.xml

echo "--- Generating Far Detector Geometry Drift in Y--- No Foil - 2PMT per M2"
# with wires for larcore
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 1 -o dphase_geo_driftY_2PMTperM2.xml --pmtdensity 2 --TPBFoils 0
perl make_legacy.gdml.pl -i dphase_geo_driftY_2PMTperM2.xml -o dunedphase10kt_v2_driftY_2PMTperM2_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_2PMTperM2_*_legacy.gdml
rm dphase_geo_driftY_2PMTperM2.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 0 -o dphase_geo_driftY_2PMTperM2_nowires.xml --pmtdensity 2 --TPBFoils 0
perl make_legacy.gdml.pl -i dphase_geo_driftY_2PMTperM2_nowires.xml -o dunedphase10kt_v2_driftY_2PMTperM2_nowires_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_2PMTperM2_nowires_*_legacy.gdml
rm dphase_geo_driftY_2PMTperM2_nowires.xml



echo "--- Generating Far Detector Geometry Drift in Y--- Half Foil - 2PMT per M2"
# with wires for larcore
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 1 -o dphase_geo_driftY_2PMTperM2_HalfFoil.xml --pmtdensity 2 --TPBFoils 2
perl make_legacy.gdml.pl -i dphase_geo_driftY_2PMTperM2_HalfFoil.xml -o dunedphase10kt_v2_driftY_2PMTperM2_HalfFoil_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_2PMTperM2_HalfFoil_*_legacy.gdml
rm dphase_geo_driftY_2PMTperM2_HalfFoil.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 0 -o dphase_geo_driftY_2PMTperM2_HalfFoil_nowires.xml --pmtdensity 2 --TPBFoils 2
perl make_legacy.gdml.pl -i dphase_geo_driftY_2PMTperM2_HalfFoil_nowires.xml -o dunedphase10kt_v2_driftY_2PMTperM2_HalfFoil_nowires_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_2PMTperM2_HalfFoil_nowires_*_legacy.gdml
rm dphase_geo_driftY_2PMTperM2_HalfFoil_nowires.xml

########################################
### Far detector workspace: 2x1 TPCs ###
########################################

echo "--- Generating Far Detector Geometry 2x1 workspace ---"

# with wires for larcore
perl generate_dunedphase10kt_v2_legacy.pl -w 1 -k 1 -o dphase_geo_workspace.xml
perl make_legacy.gdml.pl -i dphase_geo_workspace.xml -o dunedphase10kt_v2_workspace_legacy.gdml

# clean up
rm dunedphase10kt_v2_workspace_*_legacy.gdml
rm dphase_geo_workspace.xml

# no wires for geant
perl generate_dunedphase10kt_v2_legacy.pl -w 0 -k 1 -o dphase_geo_workspace_nowires.xml
perl make_legacy.gdml.pl -i dphase_geo_workspace_nowires.xml -o dunedphase10kt_v2_workspace_nowires_legacy.gdml

# clean up
rm dunedphase10kt_v2_workspace_nowires_*_legacy.gdml
rm dphase_geo_workspace_nowires.xml


#####################################################
### Far detector workspace: 2x1 TPCs (drift in Y) ###
#####################################################

echo "--- Generating Far Detector Geometry 2x1 workspace Drift in Y ---"

# with wires for larcore
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 1 -k 1 -o dphase_geo_driftY_workspace.xml --pmtdensity 1 --TPBFoils 0
perl make_legacy.gdml.pl -i dphase_geo_driftY_workspace.xml -o dunedphase10kt_v2_driftY_workspace_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_workspace_*_legacy.gdml
rm dphase_geo_driftY_workspace.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 0 -k 1 -o dphase_geo_driftY_workspace_nowires.xml --pmtdensity 1 --TPBFoils 0
perl make_legacy.gdml.pl -i dphase_geo_driftY_workspace_nowires.xml -o dunedphase10kt_v2_driftY_workspace_nowires_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_workspace_nowires_*_legacy.gdml
rm dphase_geo_driftY_workspace_nowires.xml

########################################
### Far detector workspace: 4x2 TPCs ###
########################################

echo "--- Generating Far Detector Geometry 4x2 workspace ---"

# with wires for larcore
perl generate_dunedphase10kt_v2_legacy.pl -w 1 -k 2 -o dphase_geo_workspace4x2.xml
perl make_legacy.gdml.pl -i dphase_geo_workspace4x2.xml -o dunedphase10kt_v2_workspace4x2_legacy.gdml

# clean up
rm dunedphase10kt_v2_workspace4x2_*_legacy.gdml
rm dphase_geo_workspace4x2.xml

# no wires for geant
perl generate_dunedphase10kt_v2_legacy.pl -w 0 -k 2 -o dphase_geo_workspace4x2_nowires.xml -w 0
perl make_legacy.gdml.pl -i dphase_geo_workspace4x2_nowires.xml -o dunedphase10kt_v2_workspace4x2_nowires_legacy.gdml

# clean up
rm dunedphase10kt_v2_workspace4x2_nowires_*_legacy.gdml
rm dphase_geo_workspace4x2_nowires.xml


#####################################################
### Far detector workspace: 4x2 TPCs (drift in Y) ###
#####################################################

echo "--- Generating Far Detector Geometry 4x2 workspace Drift in Y ---"

# with wires for larcore
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 1 -k 2 -o dphase_geo_driftY_workspace4x2.xml --pmtdensity 1 --TPBFoils 0
perl make_legacy.gdml.pl -i dphase_geo_driftY_workspace4x2.xml -o dunedphase10kt_v2_driftY_workspace4x2_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_workspace4x2_*_legacy.gdml
rm dphase_geo_driftY_workspace4x2.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 0 -k 2 -o dphase_geo_driftY_workspace4x2_nowires.xml -w 0 --pmtdensity 1 --TPBFoils 0
perl make_legacy.gdml.pl -i dphase_geo_driftY_workspace4x2_nowires.xml -o dunedphase10kt_v2_driftY_workspace4x2_nowires_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_workspace4x2_nowires_*_legacy.gdml
rm dphase_geo_driftY_workspace4x2_nowires.xml


#####################################################
### Far detector workspace: 3x3 TPCs (drift in Y) ###
#####################################################

echo "--- Generating Far Detector Geometry 3x3 workspace Drift in Y ---"

# with wires for larcore
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 1 -k 3 -o dphase_geo_driftY_workspace3x3.xml --pmtdensity 1 --TPBFoils 0
perl make_legacy.gdml.pl -i dphase_geo_driftY_workspace3x3.xml -o dunedphase10kt_v2_driftY_workspace3x3_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_workspace3x3_*_legacy.gdml
rm dphase_geo_driftY_workspace3x3.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 0 -k 3 -o dphase_geo_driftY_workspace3x3_nowires.xml -w 0 --pmtdensity 1 --TPBFoils 0
perl make_legacy.gdml.pl -i dphase_geo_driftY_workspace3x3_nowires.xml -o dunedphase10kt_v2_driftY_workspace3x3_nowires_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_workspace3x3_nowires_*_legacy.gdml
rm dphase_geo_driftY_workspace3x3_nowires.xml


####################
### protoDUNE-DP ###
####################

echo "--- Generating protoDUNEDP Geometry ---"

# with wires for larcore
perl generate_protodunedphase_legacy.pl -w 1 -opt 0 -o protodunedphase_geo.xml
perl make_legacy.gdml.pl -i protodunedphase_geo.xml -o protodunedphase_legacy.gdml

# clean up
rm protodunedphase_geo.xml
rm protodunedphase_*_legacy.gdml

# no wires for geant
perl generate_protodunedphase_legacy.pl -w 0 -opt 0 -o protodunedphase_geo_nowires.xml
perl make_legacy.gdml.pl -i protodunedphase_geo_nowires.xml -o protodunedphase_nowires_legacy.gdml

# clean up
rm protodunedphase_geo_nowires.xml
rm protodunedphase_nowires_*_legacy.gdml


#########################################
### protoDUNE-DP rotated (drift in Y) ###
#########################################

echo "--- Generating protoDUNEDP Geometry Drift in Y ---"

# with wires for larcore
perl generate_protodunedphase_driftY_legacy.pl -w 1 -opt 0 -o protodunedphase_driftY_geo.xml -pmt 0
perl make_legacy.gdml.pl -i protodunedphase_driftY_geo.xml -o protodunedphase_driftY_legacy.gdml

# clean up
rm protodunedphase_driftY_geo.xml
rm protodunedphase_driftY_*_legacy.gdml

# no wires for geant
perl generate_protodunedphase_driftY_legacy.pl -w 0 -opt 0 -o protodunedphase_driftY_nowires_geo.xml -pmt 0
perl make_legacy.gdml.pl -i protodunedphase_driftY_nowires_geo.xml -o protodunedphase_driftY_nowires_legacy.gdml

# clean up
rm protodunedphase_driftY_nowires_geo.xml
rm protodunedphase_driftY_nowires_*_legacy.gdml


echo "--- Generating protoDUNEDP Geometry Drift in Y - Uniform PMT Layout - to be deprecated ---"

# with wires for larcore
perl generate_protodunedphase_driftY_legacy.pl -w 1 -opt 0 -o protodunedphase_driftY_uniform_geo.xml -pmt 1
perl make_legacy.gdml.pl -i protodunedphase_driftY_uniform_geo.xml -o protodunedphase_driftY_uniform_legacy.gdml

# clean up
rm protodunedphase_driftY_uniform_geo.xml
rm protodunedphase_driftY_uniform_*_legacy.gdml

# no wires for geant
perl generate_protodunedphase_driftY_legacy.pl -w 0 -opt 0 -o protodunedphase_driftY_uniform_nowires_geo.xml -pmt 1
perl make_legacy.gdml.pl -i protodunedphase_driftY_uniform_nowires_geo.xml -o protodunedphase_driftY_uniform_nowires_legacy.gdml

# clean up
rm protodunedphase_driftY_uniform_nowires_geo.xml
rm protodunedphase_driftY_uniform_nowires_*_legacy.gdml


##################################################
### protoDUNE-DP for photon library generation ###
##################################################

echo "--- Generating protoDUNEDP OPTICAL Geometry ---"

# all detail to generate photon libraries
perl generate_protodunedphase_legacy.pl -w 1 -opt 1 -o protodunedphase_optical_geo.xml
perl make_legacy.gdml.pl -i protodunedphase_optical_geo.xml -o protodunedphase_optical_legacy.gdml

# clean up
rm protodunedphase_optical_geo.xml
rm protodunedphase_optical_*_legacy.gdml

# no wires for geant
perl generate_protodunedphase_legacy.pl -w 0 -opt 1 -o protodunedphase_optical_nowires_geo.xml
perl make_legacy.gdml.pl -i protodunedphase_optical_nowires_geo.xml -o protodunedphase_optical_nowires_legacy.gdml

# clean up
rm protodunedphase_optical_nowires_geo.xml
rm protodunedphase_optical_nowires_*_legacy.gdml

echo "--- Generating protoDUNEDP OPTICAL Geometry Drift in Y ---"

# all detail to generate photon libraries
perl generate_protodunedphase_driftY_legacy.pl -w 1 -opt 1 -o protodunedphase_driftY_optical_geo.xml -pmt 0
perl make_legacy.gdml.pl -i protodunedphase_driftY_optical_geo.xml -o protodunedphase_driftY_optical_legacy.gdml

# clean up
rm protodunedphase_driftY_optical_geo.xml
rm protodunedphase_driftY_optical_*_legacy.gdml

# no wires for geant
perl generate_protodunedphase_driftY_legacy.pl -w 0 -opt 1 -o protodunedphase_driftY_optical_nowires_geo.xml -pmt 0
perl make_legacy.gdml.pl -i protodunedphase_driftY_optical_nowires_geo.xml -o protodunedphase_driftY_optical_nowires_legacy.gdml

# clean up
rm protodunedphase_driftY_optical_nowires_geo.xml
rm protodunedphase_driftY_optical_nowires_*_legacy.gdml

########################
### 3x1x1 dual phase ###
########################

echo "--- Generating 3x1x1 Geometry ---"

# with wires for larcore
perl generate_3x1x1dphase_legacy.pl -w 1 -o 3x1x1dphase_geo.xml
perl make_legacy.gdml.pl -i 3x1x1dphase_geo.xml -o 3x1x1dphase_legacy.gdml

# clean up
rm 3x1x1dphase_geo.xml
rm 3x1x1dphase_*_legacy.gdml

# no wires for geant
perl generate_3x1x1dphase_legacy.pl -w 0 -o 3x1x1dphase_geo_nowires.xml
perl make_legacy.gdml.pl -i 3x1x1dphase_geo_nowires.xml -o 3x1x1dphase_nowires_legacy.gdml

# clean up
rm 3x1x1dphase_geo_nowires.xml
rm 3x1x1dphase_nowires_*_legacy.gdml



#####################################
### 3x1x1 dual phase (drift in Y) ###
#####################################

echo "--- Generating 3x1x1 Geometry Drift in Y ---"

# with wires for larcore
perl generate_3x1x1dphase_driftY_legacy.pl -w 1 -o 3x1x1dphase_driftY_geo.xml
perl make_legacy.gdml.pl -i 3x1x1dphase_driftY_geo.xml -o 3x1x1dphase_driftY_legacy.gdml

# clean up
rm 3x1x1dphase_driftY_geo.xml
rm 3x1x1dphase_driftY_*_legacy.gdml

# no wires for geant
perl generate_3x1x1dphase_driftY_legacy.pl -w 0 -o 3x1x1dphase_driftY_geo_nowires.xml
perl make_legacy.gdml.pl -i 3x1x1dphase_driftY_geo_nowires.xml -o 3x1x1dphase_driftY_nowires_legacy.gdml

# clean up
rm 3x1x1dphase_driftY_geo_nowires.xml
rm 3x1x1dphase_driftY_nowires_*_legacy.gdml

