#!/bin/bash

args=("$@")
RunDir=${args[0]}
infile=${args[1]}
LPDir=${args[2]}

INFILE=`basename $infile`

fileend=.root
outhistfile=${INFILE%$fileend}_nearline_hist.root

export LOCKFILE=$INFILE.LOCK
export DONEFILE=$INFILE.DONE

# source /home/mbaird42/.bashrc # may be necessary later to have access to setup functions



# Touch a lock file...
echo "Creating file $infile.LOCK"
touch $RunDir/$LOCKFILE

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
    lar -c test_daqtooffline_nearlineana.fcl ${infile} -T $outhistfile
fi



# Remove lock file...
rm -v $RunDir/$LOCKFILE    

# Touch done file...
touch $RunDir/$DONEFILE

echo ""
echo "Done with file $infile..."
echo ""
echo ""