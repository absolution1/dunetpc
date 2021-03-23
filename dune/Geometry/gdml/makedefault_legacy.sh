# make sure that:
# $Pitch3mmVersion = 0;
# $UVAngle45Option = 0;


# full geometry, no wires
perl generate_dune10kt_v1_legacy.pl -o dune10kt_v1_nowires_legacy.gdml -w=0
perl make_legacy.gdml.pl -i dune10kt_v1_nowires_legacy.gdml -o dune10kt_v1_nowires_legacy.gdml
# full geometry, with wires
perl generate_dune10kt_v1_legacy.pl -o dune10kt_v1_legacy.gdml
perl make_legacy.gdml.pl -i dune10kt_v1_legacy.gdml -o dune10kt_v1_legacy.gdml


# workspace geometry, no wires
perl generate_dune10kt_v1_legacy.pl -o dune10kt_v1_workspace_nowires_legacy.gdml -w=0 -k=1
perl make_legacy.gdml.pl -i dune10kt_v1_workspace_nowires_legacy.gdml -o dune10kt_v1_workspace_nowires_legacy.gdml
# workspace geometry, with wires
perl generate_dune10kt_v1_legacy.pl -o dune10kt_v1_workspace_legacy.gdml -k=1
perl make_legacy.gdml.pl -i dune10kt_v1_workspace_legacy.gdml -o dune10kt_v1_workspace_legacy.gdml


# 1x2x6 geometry, no wires
perl generate_dune10kt_v1_legacy.pl -o dune10kt_v1_1x2x6_nowires_legacy.gdml -w=0 -k=2
perl make_legacy.gdml.pl -i dune10kt_v1_1x2x6_nowires_legacy.gdml -o dune10kt_v1_1x2x6_nowires_legacy.gdml
# 1x2x6 geometry, with wires
perl generate_dune10kt_v1_legacy.pl -o dune10kt_v1_1x2x6_legacy.gdml -k=2
perl make_legacy.gdml.pl -i dune10kt_v1_1x2x6_legacy.gdml -o dune10kt_v1_1x2x6_legacy.gdml

# 1x2x6 geometry v2, no wires
perl generate_dune10kt_v2_legacy.pl -o dune10kt_v2_1x2x6_nowires_legacy.gdml -w=0 -k=2
perl make_legacy.gdml.pl -i dune10kt_v2_1x2x6_nowires_legacy.gdml -o dune10kt_v2_1x2x6_nowires_legacy.gdml
# 1x2x6 geometry v2, with wires
perl generate_dune10kt_v2_legacy.pl -o dune10kt_v2_1x2x6_legacy.gdml -k=2
perl make_legacy.gdml.pl -i dune10kt_v2_1x2x6_legacy.gdml -o dune10kt_v2_1x2x6_legacy.gdml

# ICEBERG, no wires
perl generate_iceberg_v1_legacy.pl -o iceberg_v1_nowires_legacy.gdml -w=0
perl make_legacy.gdml.pl -i iceberg_v1_nowires_legacy.gdml -o iceberg_v1_nowires_legacy.gdml
# full geometry, with wires
perl generate_iceberg_v1_legacy.pl -o iceberg_v1_legacy.gdml
perl make_legacy.gdml.pl -i iceberg_v1_legacy.gdml -o iceberg_v1_legacy.gdml

#This should also run deacrylify to make sure new geometries are consistent.
. deacrylify10kt_legacy.sh
