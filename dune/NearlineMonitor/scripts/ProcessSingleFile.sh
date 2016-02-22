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
outhistfile=${INFILE%$fileend}_nearline_hist.root
outhistfilemuon=${INFILE%$fileend}_nearline_muon_counters.root

export LOCKFILE=$INFILE.LOCK
export DONEFILE=$INFILE.DONE

START_DATE=`date`

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
    echo "Processing /data/lbnedaq/data/nearline-monitoring-links/${INFILE}"

    PEDESTAL=`cat /data/lbnedaq/pedestals/current_run.txt`
    NEARLINE_PEDESTAL=/data/lbnedaq/pedestals/database_pedestals/offline_databaseRun_${PEDESTAL}.csv
    if [ -e $NEARLINE_PEDESTAL ];then
	export NEARLINE_PEDESTAL=$NEARLINE_PEDESTAL
    else
	export NEARLINE_PEDESTAL="/home/lbnedaq/nearline/pedestal_files/offline_databaseRun_9754.csv"
    fi

    echo "Setting pedestal to: $NEARLINE_PEDESTAL"

    lar -c test_stitcher_nearlineana.fcl -n 10 /data/lbnedaq/data/nearline-monitoring-links/${INFILE} -T $outhistfile
    
    END_NEARLINE_ANA=`date`

    lar -c nearline_muoncounter35t.fcl /data/lbnedaq/data/nearline-monitoring-links/${INFILE} -T $outhistfilemuon

    END_NEARLINE_MUON=`date`

fi



# Remove lock file...
rm -v $RunDir/$LOCKFILE    

# Remove hard link...
rm -f /data/lbnedaq/data/nearline-monitoring-links/${INFILE}

END_DATE=`date`

# Touch done file...
touch $RunDir/$DONEFILE

echo "START_DATE $START_DATE" >> $RunDir/$DONEFILE
echo "END_NEARLINE_ANA $END_NEARLINE_ANA" >> $RunDir/$DONEFILE
echo "END_NEARLINE_MUON $END_NEARLINE_MUON" >> $RunDir/$DONEFILE
echo "END_DATE $END_DATE" >> $RunDir/$DONEFILE
echo "NEARLINE_PEDESTAL $NEARLINE_PEDESTAL" >> $RunDir/$DONEFILE

echo ""
echo "Done with file $infile..."
echo ""
echo ""
