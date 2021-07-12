# protodune geometry v8, no wires
perl generate_protodune-sp_v8.pl -o protodune_v8_nowires.gdml -w=0 -p=1
perl make_gdml.pl -i protodune_v8_nowires.gdml -o protodune_v8_nowires.gdml
# protodune geometry v8, with wires
perl generate_protodune-sp_v8.pl -o protodune_v8.gdml -p=1
perl make_gdml.pl -i protodune_v8.gdml -o protodune_v8.gdml

. deacrylify10kt.sh
