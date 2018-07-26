#!/bin/sh

# Tom Junk, August 2017
# import data from wa105 raw .dat files to larsoft root files.
# since the CPU usage is so slight, import all of the files in one job.  If the job does not finish, it can
# be resubmitted and it will pick up where it left off.  If we need to parallelize this, then the logic will have
# to be redone a bit.

# set testmode to 1 to test interactively
testmode=0

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunetpc v06_49_00 -q e14:prof

# change this when we install pedestals in dune_pardata:
data_path_to_pedestal=/pnfs/dune/persistent/users/trj/wa105pedestals/pedestal_run729_1.ped
path_to_dropbox=/pnfs/dune/scratch/dunepro/dropbox/data

# Tom's interactive testing

if [[ $testmode = 1 ]]; then
 _CONDOR_SCRATCH_DIR=/dune/data2/users/trj/wa105testprod/
 path_to_dropbox=/pnfs/dune/scratch/users/trj/batchtest/
fi
path_to_working_dir=$_CONDOR_SCRATCH_DIR

mkdir -p $path_to_working_dir
cd $path_to_working_dir
ifdh cp -D ${data_path_to_pedestal} .
path_to_pedestal=`basename ${data_path_to_pedestal}`

touch dropboxlist.txt
rm dropboxlist.txt
ifdh ls ${path_to_dropbox} > dropboxlist.txt

touch alreadyinsamlist.txt
rm alreadyinsamlist.txt
samweb list-files file_name=wa105%.root and data_tier=raw and lbne_data.name=wa105_testdata_2017 and version=v06_49_00 > alreadyinsamlist.txt

for file in `samweb list-files file_name="wa105%.dat" | grep -v recotask | grep -v symlink`
do
  echo Processing: $file

# check to see if we already have this file, either in SAM or waiting in the dropbox

  filebase=`basename $file .dat`
  outputfile=${filebase}_importpass2.root
  if [[ x`grep ${outputfile} dropboxlist.txt` != x ]]; then
    continue
  fi
# efficiency speedup -- use the grep of just one samweb command output instead of querying SAM for every file
  if [[ x`grep ${outputfile} alreadyinsamlist.txt` != x ]]; then
    continue
  fi
#  if [[ x`samweb list-files file_name=${outputfile}`  != x ]]; then
#    continue
#  fi


  fileurl=`samweb get-file-access-url $file`
  ifdh cp -D $fileurl .

  tmpfilename=`echo $file | sed -e "s/wa105\_/wa105/" | sed -e "s/\_s/-/" | sed -e "s/\_.*.dat/.dat/" | sed -e "s/wa105/wa105\_/"`
  ln -s $file $tmpfilename

  tmp2=${tmpfilename#*r}
  runnum=${tmp2%-*}
  echo run number $runnum

  tmp1=${tmpfilename#*-}
  subrun=${tmp1%.*}
  echo subrun $subrun

  touch tmpfcl.fcl
  rm tmpfcl.fcl
cat > tmpfcl.fcl <<EOF
#include "ImportFull311File.fcl"
source.PedestalFile: "$path_to_pedestal"
outputs.out1.compressionLevel: 1
EOF

  lar -c tmpfcl.fcl $tmpfilename -o $outputfile

  jsonfile=${outputfile}.json
  touch $jsonfile
  rm $jsonfile
  echo Making $jsonfile

  size=`stat -c %s $outputfile`
  nev=`echo "Events->GetEntriesFast()" | root -b -l $outputfile 2>&1 | tail -1 | cut -d')' -f2`
  nevm1=$(($nev-1))
  echo Number of events: $nev and minus 1: $nevm1

  dunetpcversion=`ups active | grep dunetpc | awk '{print $2}'`

cat >> $jsonfile <<EOF
{
 "file_name": "${outputfile}",
 "file_size": ${size},
 "file_type": "test-data",
 "event_count": ${nev},
 "first_event": 0,
 "last_event": ${nevm1},
 "runs": [ [ $runnum, "test" ] ],
 "Online.Run": ${runnum},
 "Online.RunNumber": ${runnum},
 "Online.Subrun": ${subrun},
 "file_format": "root",
 "data_tier": "raw",
 "group": "dune",
 "application": {
  "family": "art",
  "name": "ImportFull311File",
  "version": "${dunetpcversion}"
  },
  "lbne_data.name": "wa105_testdata_2017",
  "lbne_data.detector_type": "3x1x1dp",
  "parents": [
    {
     "file_name": "$file"
    }
   ]
}
EOF

  if [[ $testmode = 0 ]]; then
    samweb declare-file $jsonfile
  fi

  if [[ x`ifdh ls ${path_to_dropbox}/${outputfile}` = x ]]; then
    ifdh cp -D ${outputfile} ${path_to_dropbox}
  fi
  if [[ x`ifdh ls ${path_to_dropbox}/${jsonfile}` = x ]]; then
    ifdh cp -D ${jsonfile} ${path_to_dropbox}
  fi

# clean up
  rm $file
  rm $tmpfilename
  rm $outputfile  
  rm $jsonfile
  rm tmpfcl.fcl

  if [[ $testmode = 1 ]]; then
    echo "Test mode. Just running one file"
    exit
  fi

done

