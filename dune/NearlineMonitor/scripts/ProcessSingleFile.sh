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
    lar -c test_daqtooffline_nearlineana.fcl ${infile} -o /home/lbnedaq/nearline/temp/${INFILE}.temp -T $outhistfile
fi

# No need to keep the art output file. Remove it...
rm -f /home/lbnedaq/nearline/temp/${INFILE}.temp



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