#!/bin/bash

args=("$@")
argssize=${#args[*]}
if [ $argssize -ne 4 ];then
    echo ""
    echo "Usage:   ./ProcessNewFiles.sh {maximum number of jobs to run} {version} {qualifier} {compiler}"
    echo ""
    echo "Example:  ./ProcessNewFiles.sh 4 v04_30_03 prof e9"
    echo ""
    echo ""
    exit
fi

maxjobs=${args[0]}
version=${args[1]}
qual=${args[2]}
comp=${args[3]}

echo ""
echo "Processing 35t nearline with maxjobs: $maxjobs"
echo ""
echo ""



# Check Disk space before making a mess
export DiskUsage=`df /home -h |  awk 'NR==2 && match($5,"%"){print substr($5,RSTART-3,1)}'`
if [ $DiskUsage -gt 95 ];then
    echo ""
    echo "Disk too full..."
    echo "nearline-data disk for 35t machine is $DiskUsage % full."
    echo "All processing will be stopped until you clear out some space."
    echo ""
    echo ""
    exit
fi



# Test to see if this job is already running.
if [ -e /tmp/Batch-35t-Nearline.LOCK ];then
    if [ $(( $(date +%s) - $(stat -c %Y /tmp/Batch-35t-Nearline.LOCK) )) -gt 14400 ];then
      echo ""
      echo ".LOCK file is older than 4 hours..."
      echo "There is a stale lock file that needs attention."
      echo "Machine: $(hostname)"
      echo "File: /tmp/Batch-35t-Nearline.LOCK"
      echo "ProcessNewFiles.sh skipped this file at $(date)."
    fi
    echo "lock file /tmp/Batch-35t-Nearline.LOCK exists. Exiting..."
    echo ""
    echo ""
    exit
fi

touch /tmp/Batch-35t-Nearline.LOCK




# Setup necessary pathways:
export ScriptPath=/home/mbaird42/nearline_test_release/srcs/dunetpc/dune/NearlineMonitor/scripts
export RelDir=/home/mbaird42/nearline_test_release/localProducts_larsoft_${version}_${qual}_${comp}
export OutputPath=/home/mbaird42/35t_nearline/output
export BaseFileName=lbne_r
export FileSearch='/home/mbaird42/tickler_data/lbne_r*.root'
export filepos=5



# Only look for new files between 10 min and 10 days old and sort them so that files
# with the largest run/subrun numbers are first in the list.
for file in $( find ${FileSearch} -mtime -10 -mmin +0 | sort -t "/" -k$filepos -r )
do {

	if [ `ps aux | grep ProcessSingleFile | wc -l` -gt $maxjobs ];then
	    echo ""
	    echo "Reached the max number of running jobs.  Exiting..."
	    echo ""
	    echo ""
	    rm -fv /tmp/Batch-35t-Nearline.LOCK
	    exit
	fi

	FILE=`basename $file`

	export bigrun=${FILE:6:3}
	export run=${FILE:6:6}
	export subrun=${FILE:15:2}
        export RunDir=$OutputPath/$bigrun/$run

	mkdir -p $RunDir

	# Check if this file has already been processed or is curently being processed.
	if [ -e $RunDir/$FILE.DONE ];then
	    #echo "SKIP: $FILE has already been processed."
	    continue
	fi
	if [ -e $RunDir/$FILE.LOCK ];then
	    #echo "SKIP: $FILE is currently being processed."
	    continue
	fi

      	# execute script to process single file
	echo "Processing $file"
	cd $ScriptPath
	nohup ./ProcessSingleFile.sh $RunDir $file $RelDir $version $qual $comp > $RunDir/$FILE.log 2>&1 &

    } done

echo ""

rm -fv /tmp/Batch-35t-Nearline.LOCK

echo ""
echo ""
