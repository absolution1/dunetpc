#! /bin/bash

### What are my baseline values?
BaseLifetime="3000"
BaseField="[ 0.5, 0.782, 1.734 ]"
BaseDiffusion="6.2e-9"

### Declare the electron lifetime arrays
declare -a Lifetimes=("1000" "2000" "3000" "5000" "8000")
declare -a LifeNames=("1ms" "2ms" "3ms" "5ms" "8ms")

### Declare the electric field arrays
declare -a FieldValue=("[ 0.25, 0.782, 1.734 ]" "[ 0.375, 0.782, 1.734 ]" "[ 0.5, 0.782, 1.734 ]")
declare -a FieldNames=("250V" "375V" "500V")

### Delcare the diffusion const arrays
declare -a Diffusion=("0" "3.1e-9" "6.2e-9" "12.4e-9")
declare -a DiffNames=("0Diff" "50Diff" "100Diff" "200Diff")

### Check that the length of the two electron lifetime arrays are equal
LifeLen=${#Lifetimes[@]}
LNamLen=${#LifeNames[@]}
if [ $LifeLen -ne $LNamLen ]; then
    echo "The size of Lifetimes is ${LifeLen} and the size of LifeNames is ${LNamLen}.....Fix it"
    exit
fi
### Check that the length of the two electron field arrays are equal
FieldLen=${#FieldValue[@]}
FiNamLen=${#FieldNames[@]}
if [ $FieldLen -ne $FiNamLen ]; then
    echo "The size of FieldValue is ${FieldLen} and the size of FieldNames is ${FiNamLen}.....Fix it"
    exit
fi
### Check that the length of the two electron field arrays are equal
DiffLen=${#Diffusion[@]}
DNamLen=${#DiffNames[@]}
if [ $DiffLen -ne $DNamLen ]; then
    echo "The size of Diffusion is ${DiffLen} and the size of DiffNames is ${DNamLen}.....Fix it"
    exit
fi

### Define where I am going to put the fcl files
MCTestFclDir=$MRB_SOURCE/dunetpc/fcl/dune35t/MCTests/
echo "Putting fcl files in this directory -> ${MCTestFclDir}"

### Define a function to make all of the fcl files ###
function MakeTheFcls(){
    ThisLife=$BaseLifetime
    ThisField=$BaseField
    ThisDiff=$BaseDiffusion
    
    LifeName="3ms"
    FieldName="500V"
    DiffName="100Diff"

    if [ $1 -eq 1 ]; then
	ThisLife=${Lifetimes[$2-1]}
	LifeName=${LifeNames[$2-1]}
    elif [ $1 -eq 2 ]; then
	ThisField=${FieldValue[$2-1]}
	FieldName=${FieldNames[$2-1]}
    elif [ $1 -eq 3 ]; then
	ThisDiff=${Diffusion[$2-1]}
	DiffName=${DiffNames[$2-1]}
    else
	echo "Something is wrong, you're asking for parameter "$2
    fi

    echo "Making files for ${ThisLife} ${LifeName}, ${ThisField} ${FieldName}, ${ThisDiff} ${DiffName}"

### Make the G4 file
### Note that this is using the counter muon GEANT4 fcl file!!!
    cat <<EOF > ${MCTestFclDir}MCTest_g4_${LifeName}_${FieldName}_${DiffName}.fcl
#include "standard_g4_dune35t_countermu.fcl"

services.DetectorPropertiesService.Electronlifetime: ${ThisLife}
services.DetectorPropertiesService.Efield: ${ThisField}
services.LArG4Parameters.LongitudinalDiffusion: ${ThisDiff}

EOF

### Make the detsim file
### Note that this is using the millblock detsim module!!!
    cat <<EOF > ${MCTestFclDir}MCTest_detsim_${LifeName}_${FieldName}_${DiffName}.fcl
#include "standard_detsim_dune35t_milliblock.fcl"

services.DetectorPropertiesService.Electronlifetime: ${ThisLife}
services.DetectorPropertiesService.Efield: ${ThisField}
services.LArG4Parameters.LongitudinalDiffusion: ${ThisDiff}

EOF
    if [ $1 -ne 4 ]; then
	echo "services.ChannelNoiseService: @local::chnoiseold" >> ${MCTestFclDir}MCTest_detsim_${LifeName}_${FieldName}_${DiffName}.fcl
    else
	if [ $2 -eq 1 ]; then
	    echo "" >> ${MCTestFclDir}MCTest_detsim_${LifeName}_${FieldName}_${DiffName}.fcl
	elif [ $2 -eq 2 ]; then
	    echo "" >> ${MCTestFclDir}MCTest_detsim_${LifeName}_${FieldName}_${DiffName}.fcl
	elif [ $2 -eq 3 ]; then
	    echo "" >> ${MCTestFclDir}MCTest_detsim_${LifeName}_${FieldName}_${DiffName}.fcl
	elif [ $2 -eq 4 ]; then
	    echo "" >> ${MCTestFclDir}MCTest_detsim_${LifeName}_${FieldName}_${DiffName}.fcl
	fi
    fi
    
### Make the G4 file
### Note that this is using the counter muon GEANT4 fcl file!!!
    cat <<EOF > ${MCTestFclDir}MCTest_reco_${LifeName}_${FieldName}_${DiffName}.fcl
#include "standard_reco_dune35tsim_milliblock.fcl"

services.DetectorPropertiesService.Electronlifetime: ${ThisLife}
services.DetectorPropertiesService.Efield: ${ThisField}

EOF

}

### Make the fcl files for the different electron lifetimes.
for (( life=1; life<${LifeLen}+1; life++ ));
do
    MakeTheFcls 1 ${life}
done

### Now make the fcl files for the different electric fields.
for (( field=1; field<${FieldLen}+1; field++ ));
do
    MakeTheFcls 2 ${field}
done

### Now make the fcl files for the different diffusion constants.
for (( diff=1; diff<${DiffLen}+1; diff++ ));
do
    MakeTheFcls 3 ${diff}
done