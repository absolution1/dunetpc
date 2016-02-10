#!/bin/bash

args=("$@")
argssize=${#args[*]}
if [ $argssize -ne 3 ];then
    echo ""
    echo "Usage:   ./ProcessNewFilesEVD.sh {maximum number of jobs to run} {version} {compiler}"
    echo ""
    echo "Example:  ./ProcessNewFilesEVD.sh 4 v04_30_03 e9"
    echo ""
    echo ""
    exit
fi

maxjobs=${args[0]}
version=${args[1]}
comp=${args[2]}

echo ""
echo "Processing 35t nearline Event Display with maxjobs: $maxjobs"
echo ""
echo ""



# Check Disk space before making a mess
export DiskUsage=`df /home -h |  awk 'NR==2 && match($5,"%"){print substr($5,RSTART-3,2)}'`
if [ $DiskUsage -gt 95 ];then
    echo ""
    echo "Disk too full..."
    echo "nearline home disk for 35t machine is $DiskUsage % full."
    echo "All processing will be stopped until you clear out some space."
    echo ""
    echo ""
    exit
fi

export DiskUsage2=`df /lbne/data2 -h | awk 'NR==3 && match($4,"%"){print substr($4,RSTART-3,2)}'`
if [ $DiskUsage2 -gt 95 ];then
    echo ""
    echo "Disk too full..."
    echo "nearline-data disk for 35t machine is $DiskUsage % full."
    echo "All processing will be stopped until you clear out some space."
    echo ""
    echo ""
    exit
fi



# Test to see if this job is already running.
if [ -e /tmp/Batch-35t-NearlineEVD.LOCK ];then
    if [ $(( $(date +%s) - $(stat -c %Y /tmp/Batch-35t-NearlineEVD.LOCK) )) -gt 14400 ];then
      echo ""
      echo ".LOCK file is older than 4 hours..."
      echo "There is a stale lock file that needs attention."
      echo "Machine: $(hostname)"
      echo "File: /tmp/Batch-35t-NearlineEVD.LOCK"
      echo "ProcessNewFilesEVD.sh skipped this file at $(date)."
    fi
    echo "lock file /tmp/Batch-35t-NearlineEVD.LOCK exists. Exiting..."
    echo ""
    echo ""
    exit
fi

touch /tmp/Batch-35t-NearlineEVD.LOCK




# Setup necessary pathways:
export RelDir=/home/lbnedaq/nearline/nearline_test_release_${version}
export ScriptPath=${RelDir}/srcs/dunetpc/dune/NearlineMonitor/scripts
export LPDir=${RelDir}/localProducts_larsoft_${version}_${comp}_prof
export OutputPath=/lbne/data2/users/lbnedaq/nearline_evd/${version}
export OutputSearchPath=/lbne/data2/users/lbnedaq/nearline_evd/v*
export BaseFileName=lbne_r
export FileSearch='/data/lbnedaq/data/transferred_files/lbne_r*.root'
export filepos=6



# Only look for new files between 10 min and 1 day old and sort them so that files
# with the largest run/subrun numbers are first in the list.
for file in $( find ${FileSearch} -mtime -0.5 -mmin +10 | sort -t "/" -k$filepos -r )
do {

	if [ `ps aux | grep ProcessSingleFileEVD | wc -l` -gt $maxjobs ];then
	    echo ""
	    echo "Reached the max number of running jobs.  Exiting..."
	    echo ""
	    echo ""
	    rm -fv /tmp/Batch-35t-NearlineEVD.LOCK
	    exit
	fi

	FILE=`basename $file`

	export bigrun=${FILE:6:3}
	export run=${FILE:6:6}
	export subrun=${FILE:15:2}
        export RunDir=$OutputPath/$bigrun/${run}_${subrun} #use subrun as well in case we have subruns
        export RunSearchDir=$OutputSearchPath/$bigrun/${run}_${subrun} #use subrun as well in case we have subruns

	mkdir -p $RunDir

	num_files=`find $RunSearchDir/${FILE}EVD.DONE 2>/dev/null | wc -l`

	# Check if this file has already been processed or is curently being processed.
	if [ $num_files -gt 0 ];then
	    # echo "SKIP: $FILE has already been processed."
	    continue
	fi

	num_files=`find $RunSearchDir/${FILE}EVD.LOCK 2>/dev/null | wc -l`

	if [ $num_files -gt 0 ];then
	    # echo "SKIP: $FILE is currently being processed."
	    continue
	fi

      	# execute script to process single file
	echo "Processing $file"
	cd $ScriptPath

	# There is too much text output in these jobs. So for now, don't pipe it
	# to a log file...
	# nohup ./ProcessSingleFile.sh $RunDir $file $LPDir > $RunDir/$FILE.log 2>&1 &
	nohup ./ProcessSingleFileEVD.sh $RunDir $file $LPDir >> /dev/null 2>&1 &

    } done

echo ""

rm -fv /tmp/Batch-35t-NearlineEVD.LOCK

echo ""
echo ""
