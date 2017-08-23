#!/bin/bash

#Model for a bash script to import .dat files into larsoft rawdata
#Intended use is to set the correct path_to_data and to pedestal (for the moment standard paths on eos @ CERN are selected)
#The script runs iteratively on all the files in the folder and change the import file accordingly
#the imput path points at the correct import file in the dunetpc repository. A file with the same model could be copied elsewhere and a different path cofigured
# A larsoft root file will be created as outout for any input file created
#Runs are arranged in subruns of 335 events each. The last subrun might have less events
#To report any problem: andrea.scarpelli@cern.ch

run=840

path_to_data="/eos/experiment/wa105/data/311/rawdata/$run"
path_to_pedestal="/eos/experiment/wa105/data/311/datafiles/pedestals/pedestals" #can be made an iterative subtraction of the pedestal if is took a pedestal for every run

importfile=$MRB_SOURCE"/dunetpc/dune/3x1x1dp/DataImport/ImportFullFile/ImportFull311File.fcl"

path_to_output="/eos/user/a/ascarpel/3x1x1/data"


for file in $path_to_data/*.dat;
do
        tmp=${file#*-}
        num=${tmp%.*}
        sed -i -e 12c'fileNames: [ "'$file'" ]' $importfile
  	sed -i -e 32c'source.PedestalFile: “ ‘$path_to_pedestal’ “ '  $importfile
        lar -c $importfile -o $path_to_output/raw_$run-$num.root

done
echo "All done"
