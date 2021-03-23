# protodune geometry v7, no wires
perl generate_protodune-sp_v7.pl -o protodune_v7_nowires.gdml -w=0 -p=1
perl make_legacy_gdml.pl -i protodune_v7_nowires.gdml -o protodune_v7_nowires.gdml
# protodune geometry v7, with wires
perl generate_protodune-sp_v7.pl -o protodune_v7.gdml -p=1
perl make_legacy_gdml.pl -i protodune_v7.gdml -o protodune_v7.gdml

. deacrylify10kt.sh
