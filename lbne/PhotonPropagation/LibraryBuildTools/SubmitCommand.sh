#!/bin/bash

script=OpticalLibraryBuild_Grid_lbne.sh
#outdir=/lbne/data/users/ahimmel/OpticalLibv4
outdir=/lbne/data/users/ahimmel/OpticalLibTest
#fcl=/lbne/app/users/ahimmel/testRel_35t_opsim/workdir/srcs/lbnecode/lbne/PhotonPropagation/LibraryBuildTools/lbne35t_buildopticallibrary_grid.fcl
fcl=/lbne/app/users/ahimmel/testRel_35t_opsim/workdir/srcs/lbnecode/lbne/PhotonPropagation/LibraryBuildTools/lbne35t_buildopticallibrary_grid_test.fcl

clientargs="--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --OS=SL5,SL6 --group=lbne -f $fcl --role=Analysis"
toolsargs="-q -g --opportunistic --OS=SL6 "
fileargs="-dROOT $outdir/root -dFCL $outdir/fcl -dLOG $outdir/log "



#Test job - jobsub_tools
#njobs=216000
#thisjob="-T $PWD/$script $njobs"
#echo "jobsub $toolsargs $fileargs $thisjob"
#jobsub $fileargs $thisjob

#real job - jobsub_tools
#njobs=500
#thisjob="-N $njobs $PWD/$script $njobs"
#echo "jobsub $toolsargs $fileargs $thisjob"
#jobsub $fileargs $thisjob

#Test job - jobsub_client
#njobs=216000
#thisjob="-T file://$PWD/$script $njobs"
#echo "jobsub_submit $clientargs $fileargs $thisjob"
#jobsub_submit $clientargs $fileargs $thisjob

#Test job 1 - jobsub_client
njobs=7200
nphotons=10
thisjob="-M -N 1 file://$PWD/$script $njobs $nphotons"

#Test job 2 - jobsub_client
#njobs=500
#nphotons=1000
#thisjob="-M -N 1 file://$PWD/$script $njobs $nphotons"

#Real job - jobsub_client
#njobs=500
#nphotons=30000
#thisjob="-N $njobs file://$PWD/$script $njobs $nphotons"

echo "jobsub_submit $clientargs $fileargs $thisjob"
jobsub_submit $clientargs $fileargs $thisjob
