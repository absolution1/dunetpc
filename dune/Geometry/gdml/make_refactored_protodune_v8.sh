# - D.R. <drivera@fnal.gov>

# refactored protodune geometry v8, no wires
perl generate_refactored_protodune-sp_v8.pl -o protodune_v8_refactored_nowires.gdml -w=0 -p=1
perl make_refactored_gdml.pl -i protodune_v8_refactored_nowires.gdml -o protodune_v8_refactored_nowires.gdml
# protodune geometry v8, with wires
perl generate_refactored_protodune-sp_v8.pl -o protodune_v8_refactored.gdml -p=1
perl make_refactored_gdml.pl -i protodune_v8_refactored.gdml -o protodune_v8_refactored.gdml

# why? 
#. deacrylify10kt.sh
