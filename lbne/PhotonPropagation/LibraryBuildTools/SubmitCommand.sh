#!/bin/bash

script=OpticalLibraryBuild_Grid_lbne.sh
outdir=/lbne/data/users/ahimmel/OpticalLibv3

clientargs="--resource-provides=usage_model=OPPORTUNISTIC --OS=SL5,SL6 --group=lbne "
toolsargs="-q -g --opportunistic --OS=SL6 "
fileargs="-dROOT $outdir/root -dFCL $outdir/fcl -dLOG $outdir/log "



#Test job - jobsub_tools
#njobs=216000
#thisjob="-T $PWD/$script $njobs"
#echo "jobsub $toolsargs $fileargs $thisjob"
#jobsub $fileargs $thisjob

#real job - jobsub_tools
njobs=500
thisjob="-N $njobs $PWD/$script $njobs"
echo "jobsub $toolsargs $fileargs $thisjob"
jobsub $fileargs $thisjob

#Test job - jobsub_client
#njobs=216000
#thisjob="-T file://$PWD/$script $njobs"
#echo "jobsub_submit $clientargs $fileargs $thisjob"
#jobsub_submit $clientargs $fileargs $thisjob

#Real job - jobsub_client
#njobs=200
#thisjob="-N $njobs file://$PWD/$script $njobs"
#echo "jobsub_submit $clientargs $fileargs $thisjob"
#jobsub_submit $clientargs $fileargs $thisjob

