#!/bin/bash

args=("$@")
argssize=${#args[*]}
if [ $argssize -ne 3 ];then
    echo ""
    echo "Three input arguments must be provided! Exiting..."
    echo ""
    echo ""
    exit
fi

RunDir=${args[0]}
infile=${args[1]}
LPDir=${args[2]}

INFILE=`basename $infile`

fileend=.root
outhistfile=${INFILE%$fileend}_evd.root

export LOCKFILE=${INFILE}EVD.LOCK
export DONEFILE=${INFILE}EVD.DONE



# Touch a lock file...
echo "Creating file ${infile}EVD.LOCK"
touch $RunDir/$LOCKFILE

# Create a hard link...
ln ${infile} /data/lbnedaq/data/nearline-monitoring-links/${INFILE}

# Setup the LArSoft environment...
echo ""
echo "Setting up LArSoft/DUNETPC:"
source /grid/fermiapp/products/dune/setup_dune.sh
cd $LPDir
source ./setup
cd ..
mrbslp
echo ""

# Move into the appropriate output directory...
cd $RunDir

export infilesize=`ls -l $infile | awk '{ print $5 }'`


# Skip files that are too small and probably DAQ junk...
if [ $infilesize -gt 500 ];
then
    echo "Processing /data/lbnedaq/data/nearline-monitoring-links/${INFILE}"
    lar -c ctreeraw35t_trigTPC.fcl -n 10 /data/lbnedaq/data/nearline-monitoring-links/${INFILE}

    # Rename the output file for the RED35 EVD so that the script that looks for the most recent file
    # knows that the processing is finished.
    mv sample.root sample_done.root

    # Don't know why this text file gets created. Probably could be turned off at the fcl level.
    # Just going to remove it manually for now...
    rm -fv WireGeometry.txt
fi



# Remove lock file...
rm -v $RunDir/$LOCKFILE    

# Remove hard link...
rm -f /data/lbnedaq/data/nearline-monitoring-links/${INFILE}

# Touch done file...
touch $RunDir/$DONEFILE

echo ""
echo "Done with file $infile..."
echo ""
echo ""
