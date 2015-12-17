#!/bin/bash

args=("$@")
RunDir=${args[0]}
infile=${args[1]}

INFILE=`basename $infile`

export LOCKFILE=$INFILE.LOCK
export DONEFILE=$INFILE.DONE

# source /home/mbaird42/.bashrc # may be necessary later to have access to setup functions



# Touch a lock file...
echo "Creating file $infile.LOCK"
touch $RunDir/$LOCKFILE

# Setup the LArSoft environment...
echo "Setting up LArSoft:"
echo "insert all setup stuff here..."

# Move into the appropriate output directory...
cd $RunDir

export infilesize=`ls -l $infile | awk '{ print $5 }'`


# Skip files that are too small and probably DAQ junk...
if [ $infilesize -gt 500 ];
then
    echo $infile
    echo "Processing $infile"
    echo "lar -c whateverjob.fcl ${infile}"
fi



# Remove lock file...
rm -v $RunDir/$LOCKFILE    

# Touch done file...
touch $RunDir/$DONEFILE

echo ""
echo "Done with file $infile..."
echo ""
echo ""