#!/bin/bash

perl generate_dunedphase10kt_v2.pl -o dphase_geo.xml -w 0
make_gdml.pl -i dphase_geo.xml -o dunedphase10kt_v2_nowires.gdml

#rm -rf dphase_geo.xml
#rm -rf dunedphase10kt_v2_*.gdml

#root -l dunedphase10kt_geo.C
