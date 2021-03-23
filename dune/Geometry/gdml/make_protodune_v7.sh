# protodune geometry v7, no wires
perl generate_protodune-sp_v7_legacy.pl -o protodune_v7_nowires_legacy.gdml -w=0 -p=1
perl make_legacy.gdml.pl -i protodune_v7_nowires_legacy.gdml -o protodune_v7_nowires_legacy.gdml
# protodune geometry v7, with wires
perl generate_protodune-sp_v7_legacy.pl -o protodune_v7_legacy.gdml -p=1
perl make_legacy.gdml.pl -i protodune_v7_legacy.gdml -o protodune_v7_legacy.gdml

. deacrylify10kt_legacy.sh
