#!/bin/bash

echo ""
echo ""
echo "Available runs with event display data:"
echo ""



export FileSearch='/dune/data2/users/lbnedaq/nearline_evd/*/*/*/sliced_pedestal.root'



# Find the available runs
for file in $( find ${FileSearch} -mmin +1 | sort)
do {
	run=${file:53:6}
	echo $run
    } done



echo ""
echo ""
echo "NOTE:   There are 3 types of event display files. This list means that at least one event display file (i.e. - root file) exists for these runs. It is possible that for a given run, not all 3 types of event display file exist (although that would be odd...)"
echo ""
echo ""