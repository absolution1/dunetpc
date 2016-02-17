#!/bin/bash

args=("$@")
argssize=${#args[*]}
if [ $argssize -ne 1 ];then
    echo ""
    echo ""
    echo ""
    echo "Usage:    ./run_35t_evd.sh {run number - in 6 digit format please!}"
    echo ""
    echo "To get the most recent run, use run number 0."
    echo ""
    echo ""
    exit
fi

run=${args[0]}


# Setup necessary pathways:
export version=v04_36_01
export comp=e9
export RelDir=/dune/app/home/duneana/35t_EventDisplay/larsoft_${version}
export LPDir=${RelDir}/localProducts_larsoft_${version}_${comp}_prof
export EVDfclPath=${RelDir}/srcs/dunetpc/dune/NearlineMonitor/evd
export InputPath=/dune/data2/users/lbnedaq/nearline_evd/${version}
export FileSearch='/dune/data2/users/lbnedaq/nearline_evd/*/*/*/sliced_pedestal.root'
export filepos=9



# Set up the LArSoft environment
echo ""
echo ""
echo "Setting up LArSoft/DUNETPC:"
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
    run=`printf %06i $run`
    bigRun=${run:0:3}
    file=${InputPath}/${bigRun}/${run}_01/sliced_pedestal.root    
fi

echo "Opening file:"
echo $file



# Run the evd!
lar -c ${EVDfclPath}/evd_dune35t_slicer.fcl $file




echo ""
echo ""
echo "Finished running the 35t EVD..."
echo ""
echo ""