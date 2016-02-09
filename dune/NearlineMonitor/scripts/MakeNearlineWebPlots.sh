#!/bin/bash

#
# Make plots for the Nearline webpage...
#

args=("$@")
argssize=${#args[*]}
if [ $argssize -ne 3 ];then
    echo ""
    echo "Usage:   ./MakeNearlineWebPlots.sh {dunetpc release: vXX_YY_ZZ} {compiler} {time period in days}"
    echo ""
    echo "Example: ./MakeNearlineWebPlots.sh v04_32_01 e9 2"
    echo ""
    exit
fi

REL=${args[0]}
COMP=${args[1]}
DAYS=${args[2]}



# Test to see if this job is already running.
if [ -e /tmp/Batch-MakeNearlineWebPlots-${REL}-${DAYS}days.lock ];then
    echo "lock file /tmp/Batch-MakeNearlineWebPlots-${REL}-${DAYS}days.lock exists.  Exiting..."
    exit
fi
touch /tmp/Batch-MakeNearlineWebPlots-${REL}-${DAYS}days.lock



# Get a ticket.
KEYTAB=/var/adm/krb5/lbnedaq.keytab
KEYUSE=`/usr/krb5/bin/klist -k ${KEYTAB} | grep FNAL.GOV | head -1 | cut -c 5- | cut -f 1 -d /`
/usr/krb5/bin/kinit -5 -A  -kt ${KEYTAB} ${KEYUSE}/dune/`hostname`@FNAL.GOV



# Specify path to nearline output files...
SearchPath="/lbne/data2/users/lbnedaq/nearline/v*/*/*/"



# Setup the LArSoft/dunetpc environment to pick up root.
export RelDir=/home/lbnedaq/nearline/nearline_test_release_${REL}
export ScriptPath=${RelDir}/srcs/dunetpc/dune/NearlineMonitor/scripts
export LPDir=${RelDir}/localProducts_larsoft_${REL}_${COMP}_prof

echo ""
echo "Setting up LArSoft/DUNETPC:"
source /grid/fermiapp/products/dune/setup_dune.sh
cd $LPDir
source ./setup
cd ..
mrbslp
cd ~/nearline/temp
echo ""



# Make temporary list of files to be processed
rm -vf 35t_${DAYS}Day_Nearline_File_List.txt
find ${SearchPath} -mtime -${DAYS} -name "lbne_*_nearline*.root" | sort -r > 35t_${DAYS}Day_Nearline_File_List.txt



# Execute root script to make plots...
root -l -b -q "${ScriptPath}/NearlinePlotMaker.C+(${DAYS})"



# Done, remove the lock file...
rm -v /tmp/Batch-MakeNearlineWebPlots-${REL}-${DAYS}days.lock
