#!/bin/bash

####################
### Far detector ###
####################

# with wires for larcore
#perl generate_dunedphase10kt_v2.pl -o dphase_geo.xml
#make_gdml.pl -i dphase_geo.xml -o dunedphase10kt_v2.gdml

# clean up
#rm dunedphase10kt_v2_*.gdml

# no wires for geant
#perl generate_dunedphase10kt_v2.pl -o dphase_geo.xml -w 0
#make_gdml.pl -i dphase_geo.xml -o dunedphase10kt_v2_nowires.gdml

# clean up
#rm dunedphase10kt_v2_nowires_*.gdml
#rm dphase_geo.xml


#root -l dunedphase10kt_geo.C


####################
### protoDUNE-DP ###
####################

# with wires for larcore
perl generate_protodunedphase.pl -w 1 -o protodunedphase_geo.xml
perl make_gdml.pl -i protodunedphase_geo.xml -o protodunedphase.gdml

# clean up
rm protodunedphase_geo.xml
rm protodunedphase_*.gdml

# no wires for geant
perl generate_protodunedphase.pl -w 0 -o protodunedphase_geo_nowires.xml
perl make_gdml.pl -i protodunedphase_geo_nowires.xml -o protodunedphase_nowires.gdml

# clean up
rm protodunedphase_geo_nowires.xml
rm protodunedphase_nowires_*.gdml


########################
### 3x1x1 dual phase ###
########################

# with wires for larcore
#perl generate_3x1x1dphase.pl -w 1 -o 3x1x1dphase_geo.xml
#perl make_gdml.pl -i 3x1x1dphase_geo.xml -o 3x1x1dphase.gdml

# clean up
#rm 3x1x1dphase_geo.xml
#rm 3x1x1dphase_*.gdml

# no wires for geant
#perl generate_3x1x1dphase.pl -w 0 -o 3x1x1dphase_geo_nowires.xml
#perl make_gdml.pl -i 3x1x1dphase_geo_nowires.xml -o 3x1x1dphase_nowires.gdml

# clean up
#rm 3x1x1dphase_geo_nowires.xml
#rm 3x1x1dphase_nowires_*.gdml

#########################################
### protoDUNE-DP rotated (drift in Y) ###
#########################################

# with wires for larcore
perl generate_protodunedphase_rot.pl -w 1 -o protodunedphase_rot_geo.xml
perl make_gdml.pl -i protodunedphase_rot_geo.xml -o protodunedphase_rot.gdml

# clean up
rm protodunedphase_rot_geo.xml
rm protodunedphase_rot_*.gdml

# no wires for geant
perl generate_protodunedphase_rot.pl -w 0 -o protodunedphase_rot_nowires_geo.xml
perl make_gdml.pl -i protodunedphase_rot_nowires_geo.xml -o protodunedphase_rot_nowires.gdml

# clean up
rm protodunedphase_rot_nowires_geo.xml
rm protodunedphase_rot_nowires_*.gdml
