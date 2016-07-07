#!/bin/bash

args=("$@")
argssize=${#args[*]}
if [ $argssize -ne 1 ];then
    echo ""
    echo ""
    echo ""
    echo "Usage:    ./run_RED35_evd.sh {run number - in 6 digit format please!}"
    echo ""
    echo "To get the most recent run, use run number 0."
    echo ""
    echo ""
    exit
fi

run=${args[0]}

ORIGINAL_PWD=$PWD

# Setup necessary pathways:
export version=v04_35_00
export comp=e9
export RelDir=/dune/app/home/duneana/35t_EventDisplay/larsoft_${version}
export LPDir=${RelDir}/localProducts_larsoft_${version}_${comp}_prof
export RED35Dir=/dune/app/home/duneana/35t_EventDisplay/RED35
export InputPath=/dune/data2/users/lbnedaq/nearline_evd/${version}
export FileSearch='/dune/data2/users/lbnedaq/nearline_evd/*/*/*/sample_done.root'
export filepos=9


# Set up the LArSoft environment
echo ""
echo ""
echo "Using ROOT provided by LArSoft and DUNETPC (we don't actually use these two):"
echo "=========================="
echo ""
source /grid/fermiapp/products/dune/setup_dune.sh
cd $LPDir
source ./setup
cd ..
mrbslp
echo ""
echo ""




# Find the most recent file
for file in $( find ${FileSearch} -mmin +1 | sort -t "/" -k$filepos -r )
do {
	# This is a sloppy way of getting only the most recent file but I don't
	# know how else to do it...      ...don't judge me...
	break
    } done


if [ $run -ne 0 ];then
    bigRun=${run:0:3}
    file=${InputPath}/${bigRun}/${run}_01/sample_done.root    
fi

echo "Opening file:"
echo $file



# Run the evd!
cd $RED35Dir/scripts
root -l run2d.C\(\"$file\"\)



echo ""
echo ""
echo "Finished running the RED35 EVD..."
echo ""
echo ""