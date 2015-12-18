#!/bin/bash

args=("$@")
RunDir=${args[0]}
infile=${args[1]}
RelDir=${args[2]}
version=${args[3]}
qual=${args[4]}
comp=${args[5]}

INFILE=`basename $infile`

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
cd $RelDir
source $RelDir/setup
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
    lar -c test_daqtooffline_nearlineana.fcl ${infile}
fi



# Remove lock file...
rm -v $RunDir/$LOCKFILE    

# Touch done file...
touch $RunDir/$DONEFILE

echo ""
echo "Done with file $infile..."
echo ""
echo ""