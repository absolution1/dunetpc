# 1x2x6 geometry v3, no wires
perl generate_dune10kt_v3.pl -o dune10kt_v3_1x2x6_nowires.gdml -w=0 -k=2
perl make_gdml.pl -i dune10kt_v3_1x2x6_nowires.gdml -o dune10kt_v3_1x2x6_nowires.gdml
# 1x2x6 geometry v3, with wires
perl generate_dune10kt_v3.pl -o dune10kt_v3_1x2x6.gdml -k=2
perl make_gdml.pl -i dune10kt_v3_1x2x6.gdml -o dune10kt_v3_1x2x6.gdml

. deacrylify10kt.sh
