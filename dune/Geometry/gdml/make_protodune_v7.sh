# - D.R. <drivera@fnal.gov>

# refactored protodune geometry v7, no wires
perl generate_refactored_protodune-sp_v7.pl -o protodune_v7_refactored_nowires.gdml -w=0 -p=1
perl make_refactored_gdml.pl -i protodune_v7_refactored_nowires.gdml -o protodune_v7_refactored_nowires.gdml
# protodune geometry v7, with wires
perl generate_refactored_protodune-sp_v7.pl -o protodune_v7_refactored.gdml -p=1
perl make_refactored_gdml.pl -i protodune_v7_refactored.gdml -o protodune_v7_refactored.gdml

# why? 
#. deacrylify10kt.sh
