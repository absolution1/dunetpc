#!/bin/bash

#################################
### Far detector (drift in Y) ###
#################################

# with wires for larcore
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 1 -o dphase_geo_driftY.xml --pmtdensity 1 --TPBFoils 1
perl make_legacy.gdml.pl -i dphase_geo_driftY.xml -o dunedphase10kt_v2_driftY_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_*_legacy.gdml
rm dphase_geo_driftY.xml

# no wires for geant
perl generate_dunedphase10kt_v2_driftY_legacy.pl -w 0 -o dphase_geo_driftY_nowires.xml --pmtdensity 1 --TPBFoils 1
perl make_legacy.gdml.pl -i dphase_geo_driftY_nowires.xml -o dunedphase10kt_v2_driftY_nowires_legacy.gdml

# clean up
rm dunedphase10kt_v2_driftY_nowires_*_legacy.gdml
rm dphase_geo_driftY_nowires.xml

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

