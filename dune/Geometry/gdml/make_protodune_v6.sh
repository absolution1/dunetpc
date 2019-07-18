# protodune geometry v6, no wires
perl generate_dune10kt_v6.pl -o protodune_v6_nowires.gdml -w=0 -p=1
perl make_gdml.pl -i protodune_v6_nowires.gdml -o protodune_v6_nowires.gdml
# protodune geometry v6, with wires
perl generate_dune10kt_v6.pl -o protodune_v6.gdml -p=1
perl make_gdml.pl -i protodune_v6.gdml -o protodune_v6.gdml

. deacrylify10kt.sh
