# 1x2x6 geometry v3, no wires
perl generate_dune10kt_v3_legacy.pl -o dune10kt_v3_1x2x6_nowires_legacy.gdml -w=0 -k=2
perl make_legacy.gdml.pl -i dune10kt_v3_1x2x6_nowires_legacy.gdml -o dune10kt_v3_1x2x6_nowires_legacy.gdml
# 1x2x6 geometry v3, with wires
perl generate_dune10kt_v3_legacy.pl -o dune10kt_v3_1x2x6_legacy.gdml -k=2
perl make_legacy.gdml.pl -i dune10kt_v3_1x2x6_legacy.gdml -o dune10kt_v3_1x2x6_legacy.gdml

. deacrylify10kt_legacy.sh
