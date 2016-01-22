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
    echo "Processing $infile"
    lar -c ctreeraw35t_trigTPC.fcl -n 10 ${infile} >> ${outhistfile}_output.txt 2>&1
    # lar -c ctreeraw35t_trigTPC.fcl -n 10 ${infile} -o $outhistfile >> ${outhistfile}_output.txt 2>&1
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
