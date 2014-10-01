#!/bin/bash

script=OpticalLibraryBuild_Grid_lbne.sh
outdir=/lbne/data/users/ahimmel/OpticalTest
njobs=200
args=$njobs

#Real job
#jobsub --opportunistic --X509_USER_PROXY /scratch/bjpjones/grid/bjpjones.uboone.proxy -g -N 9375 -dOUT /uboone/data/users/bjpjones/OpticalProduction  -q OpticalLibraryBuild_Grid.sh `whoami` `pwd`

#Partial job
#jobsub --opportunistic --X509_USER_PROXY /scratch/bjpjones/grid/bjpjones.uboone.proxy -g -N 1500 -dOUT /uboone/data/users/bjpjones/OpticalProduction  -q OpticalLibraryBuild_Grid.sh `whoami` `pwd`


echo "jobsub --opportunistic -g -N $njobs -q -dROOT $outdir/root -dFCL $outdir/fcl -dLOG $outdir/log $script $args"
jobsub --opportunistic -g -N $njobs -q -dROOT $outdir/root -dFCL $outdir/fcl -dLOG $outdir/log $script $args

#Test job
#echo "jobsub -T -q -dROOT $outdir/root -dFCL $outdir/fcl -dLOG $outdir/log $script $args"
#jobsub -T -q -dROOT $outdir/root -dFCL $outdir/fcl -dLOG $outdir/log $script $args

