#!/bin/bash

####################
### Far detector ###
####################

# with wires for larcore
#perl generate_dunedphase10kt_v2.pl -w 1 -o dphase_geo.xml
#perl make_gdml.pl -i dphase_geo.xml -o dunedphase10kt_v2.gdml

# clean up
#rm dunedphase10kt_v2_*.gdml
#rm dphase_geo.xml

# no wires for geant
#perl generate_dunedphase10kt_v2.pl -w 0 -o dphase_geo_nowires.xml
#perl make_gdml.pl -i dphase_geo_nowires.xml -o dunedphase10kt_v2_nowires.gdml

# clean up
#rm dunedphase10kt_v2_nowires_*.gdml
#rm dphase_geo_nowires.xml


#################################
### Far detector (drift in Y) ###
#################################

# with wires for larcore
#perl generate_dunedphase10kt_v2_driftY.pl -w 1 -o dphase_geo_driftY.xml
#perl make_gdml.pl -i dphase_geo_driftY.xml -o dunedphase10kt_v2_driftY.gdml

# clean up
#rm dunedphase10kt_v2_driftY_*.gdml
#rm dphase_geo_driftY.xml

# no wires for geant
#perl generate_dunedphase10kt_v2_driftY.pl -w 0 -o dphase_geo_driftY_nowires.xml
#perl make_gdml.pl -i dphase_geo_driftY_nowires.xml -o dunedphase10kt_v2_driftY_nowires.gdml

# clean up
#rm dunedphase10kt_v2_driftY_nowires_*.gdml
#rm dphase_geo_driftY_nowires.xml


########################################
### Far detector workspace: 2x1 TPCs ###
########################################

# with wires for larcore
#perl generate_dunedphase10kt_v2.pl -w 1 -k 1 -o dphase_geo_workspace.xml
#perl make_gdml.pl -i dphase_geo_workspace.xml -o dunedphase10kt_v2_workspace.gdml

# clean up
#rm dunedphase10kt_v2_workspace_*.gdml
#rm dphase_geo_workspace.xml

# no wires for geant
#perl generate_dunedphase10kt_v2.pl -w 0 -k 1 -o dphase_geo_workspace_nowires.xml
#perl make_gdml.pl -i dphase_geo_workspace_nowires.xml -o dunedphase10kt_v2_workspace_nowires.gdml

# clean up
#rm dunedphase10kt_v2_workspace_nowires_*.gdml
#rm dphase_geo_workspace_nowires.xml


#####################################################
### Far detector workspace: 2x1 TPCs (drift in Y) ###
#####################################################

# with wires for larcore
#perl generate_dunedphase10kt_v2_driftY.pl -w 1 -k 1 -o dphase_geo_driftY_workspace.xml
#perl make_gdml.pl -i dphase_geo_driftY_workspace.xml -o dunedphase10kt_v2_driftY_workspace.gdml

# clean up
#rm dunedphase10kt_v2_driftY_workspace_*.gdml
#rm dphase_geo_driftY_workspace.xml

# no wires for geant
#perl generate_dunedphase10kt_v2_driftY.pl -w 0 -k 1 -o dphase_geo_driftY_workspace_nowires.xml
#perl make_gdml.pl -i dphase_geo_driftY_workspace_nowires.xml -o dunedphase10kt_v2_driftY_workspace_nowires.gdml

# clean up
#rm dunedphase10kt_v2_driftY_workspace_nowires_*.gdml
#rm dphase_geo_driftY_workspace_nowires.xml


########################################
### Far detector workspace: 4x2 TPCs ###
########################################

# with wires for larcore
#perl generate_dunedphase10kt_v2.pl -w 1 -k 2 -o dphase_geo_workspace4x2.xml
#perl make_gdml.pl -i dphase_geo_workspace4x2.xml -o dunedphase10kt_v2_workspace4x2.gdml

# clean up
#rm dunedphase10kt_v2_workspace4x2_*.gdml
#rm dphase_geo_workspace4x2.xml

# no wires for geant
#perl generate_dunedphase10kt_v2.pl -w 0 -k 2 -o dphase_geo_workspace4x2_nowires.xml -w 0
#perl make_gdml.pl -i dphase_geo_workspace4x2_nowires.xml -o dunedphase10kt_v2_workspace4x2_nowires.gdml

# clean up
#rm dunedphase10kt_v2_workspace4x2_nowires_*.gdml
#rm dphase_geo_workspace4x2_nowires.xml


#####################################################
### Far detector workspace: 4x2 TPCs (drift in Y) ###
#####################################################

# with wires for larcore
#perl generate_dunedphase10kt_v2_driftY.pl -w 1 -k 2 -o dphase_geo_driftY_workspace4x2.xml
#perl make_gdml.pl -i dphase_geo_driftY_workspace4x2.xml -o dunedphase10kt_v2_driftY_workspace4x2.gdml

# clean up
#rm dunedphase10kt_v2_driftY_workspace4x2_*.gdml
#rm dphase_geo_driftY_workspace4x2.xml

# no wires for geant
#perl generate_dunedphase10kt_v2_driftY.pl -w 0 -k 2 -o dphase_geo_driftY_workspace4x2_nowires.xml -w 0
#perl make_gdml.pl -i dphase_geo_driftY_workspace4x2_nowires.xml -o dunedphase10kt_v2_driftY_workspace4x2_nowires.gdml

# clean up
#rm dunedphase10kt_v2_driftY_workspace4x2_nowires_*.gdml
#rm dphase_geo_driftY_workspace4x2_nowires.xml

####################
### protoDUNE-DP ###
####################

# with wires for larcore
#perl generate_protodunedphase.pl -w 1 -opt 0 -o protodunedphase_geo.xml
#perl make_gdml.pl -i protodunedphase_geo.xml -o protodunedphase.gdml

# clean up
#rm protodunedphase_geo.xml
#rm protodunedphase_*.gdml

# no wires for geant
#perl generate_protodunedphase.pl -w 0 -opt 0 -o protodunedphase_geo_nowires.xml
#perl make_gdml.pl -i protodunedphase_geo_nowires.xml -o protodunedphase_nowires.gdml

# clean up
#rm protodunedphase_geo_nowires.xml
#rm protodunedphase_nowires_*.gdml


#########################################
### protoDUNE-DP rotated (drift in Y) ###
#########################################

# with wires for larcore
#perl generate_protodunedphase_rot.pl -w 1 -o protodunedphase_rot_geo.xml
#perl make_gdml.pl -i protodunedphase_rot_geo.xml -o protodunedphase_rot.gdml

# clean up
#rm protodunedphase_rot_geo.xml
#rm protodunedphase_rot_*.gdml

# no wires for geant
#perl generate_protodunedphase_rot.pl -w 0 -o protodunedphase_rot_nowires_geo.xml
#perl make_gdml.pl -i protodunedphase_rot_nowires_geo.xml -o protodunedphase_rot_nowires.gdml

# clean up
#rm protodunedphase_rot_nowires_geo.xml
#rm protodunedphase_rot_nowires_*.gdml


##################################################
### protoDUNE-DP for photon library generation ###
##################################################

# all detail to generate photon libraries
#perl generate_protodunedphase.pl -w 1 -opt 1 -o protodunedphase_geo_optical.xml
#perl make_gdml.pl -i protodunedphase_geo_optical.xml -o protodunedphase_optical.gdml

# clean up
#rm protodunedphase_geo_optical.xml
#rm protodunedphase_optical_*.gdml

# no wires for geant
#perl generate_protodunedphase.pl -w 0 -opt 1 -o protodunedphase_geo_optical_nowires.xml
#perl make_gdml.pl -i protodunedphase_geo_optical_nowires.xml -o protodunedphase_optical_nowires.gdml

# clean up
#rm protodunedphase_geo_optical_nowires.xml
#rm protodunedphase_optical_nowires_*.gdml



########################
### 3x1x1 dual phase ###
########################

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

