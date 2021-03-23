# 1x2x6 geometry v4, no wires
perl generate_dune10kt_v4_legacy.pl -o dune10kt_v4_1x2x6_nowires_legacy.gdml -w=0 -k=2
perl make_legacy.gdml.pl -i dune10kt_v4_1x2x6_nowires_legacy.gdml -o dune10kt_v4_1x2x6_nowires_legacy.gdml
# 1x2x6 geometry v4, with wires
perl generate_dune10kt_v4_legacy.pl -o dune10kt_v4_1x2x6_legacy.gdml -k=2
perl make_legacy.gdml.pl -i dune10kt_v4_1x2x6_legacy.gdml -o dune10kt_v4_1x2x6_legacy.gdml

#This should also run deacrylify to make sure new geometries are consistent.
. deacrylify10kt_legacy.sh
