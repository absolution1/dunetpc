#!/bin/bash
#
# A script to run the optical library building job
#
#
# To run this job:
#
# jobsub -N [NoOfJobs] -M -dROOT /out/dir/root -dFCL /out/dir/fcl -dLOG /out/dir/log OpticalLibraryBuild_Grid_dune.sh
#
# You will get outputs in the area specified by the "outstage" variable 
# which is specified below.
#
# The form of the output is one file for each few voxels. These then need 
# stitching together, which is done after all jobs are done, with a
# dedicated stitching script.
#

#missing="4070 4117"
#
#
##
#if ! echo $missing | grep -w "$PROCESS" > /dev/null; then
#  echo "Exiting since $PROCESS is already complete."
#  exit 0
#fi


#
# Set up our environment
#

njobs=$1
label=${CLUSTER}_$(printf '%04d' $PROCESS)
## This sets all the needed FW and SRT and LD_LIBRARY_PATH envt variables. 
## Then we cd back to our TMP area. EC, 23-Nov-2010.

# In each voxel, run this many photons:
NPhotonsPerVoxel=$2

umask 0002

export GROUP=dune
export HOME=$CONDOR_DIR_ROOT
export CONDOR_SCRATCH_DIR=$TEMP

LOG=${CONDOR_DIR_LOG}/pd_library_gen_${label}.log
FCL=${CONDOR_DIR_FCL}/pd_library_gen_${label}.fcl
RTF=${CONDOR_DIR_ROOT}/pd_library_gen_${label}.root

#
# Prepare to run
#
cd $CONDOR_DIR_ROOT
          echo "Writing log to ${LOG}"
          touch ${LOG}
          date 1>> ${LOG} 2>&1
          


#
# Library building parameters
#

# Copy fcl file and configure for this PROCESS 
          echo "Create this job's fhicl file" 1>> ${LOG} 2>&1
          mv -v ${CONDOR_DIR_INPUT}/*.fcl $FCL 1>> ${LOG} 2>&1
          
NX=`awk '/NX/{ print $2 }' $FCL`
          NY=`awk '/NY/{ print $2 }' $FCL`
          NZ=`awk '/NZ/{ print $2 }' $FCL`
          
# Total number of voxels
          #NTopVoxel=216000
          NTopVoxel=`echo "$NX*$NY*$NZ" | bc`
          
echo "Voxels: $NTopVoxel = $NX * $NY * $NZ" 1>> ${LOG} 2>&1
          




# In each grid job, do this many voxels:
NVoxelsPerJob=`echo "$NTopVoxel/$njobs" | bc`
          
# This works out which voxels this job should focus on: 
FirstVoxel=`echo "($NVoxelsPerJob * $PROCESS ) % $NTopVoxel" | bc` 
          LastVoxel=`echo "(($NVoxelsPerJob * $PROCESS ) + $NVoxelsPerJob - 1 ) % $NTopVoxel" | bc`
          
# Construct the run-time configuration file for this job;
# Make the random number seeds a function of the PROCESS number.
generatorSeed=$(( $PROCESS *23 + 31))
          g4Seed=$(( $PROCESS *41 + 37))
          

echo "physics.producers.generator.FirstVoxel: $FirstVoxel" >> $FCL
          echo "physics.producers.generator.LastVoxel: $LastVoxel"   >> $FCL
          echo "physics.producers.generator.N: $NPhotonsPerVoxel"    >> $FCL
          




# And then tell the user about it:
          echo "This job will run from voxel $FirstVoxel to $LastVoxel, generating $NPhotonsPerVoxel in each" 1>> ${LOG} 2>&1
          echo "CLUSTER:    " $CLUSTER 1>> ${LOG} 2>&1
          echo "PROCESS:    " $PROCESS 1>> ${LOG} 2>&1
          echo "PWD:        " $PWD     1>> ${LOG} 2>&1
          

# No need to set random seeds - generated from machine state
#echo "physics.producers.generator.RandomSeed: $generatorSeed">> $FCL
#echo "physics.producers.largeant.RandomSeed: $g4Seed">> $FCL




echo "> Setup $GROUP environment"                                                              1>> ${LOG} 2>&1
          echo "> source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh"                    1>> ${LOG} 2>&1
          source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh                             1>> ${LOG} 2>&1
          
echo "> INPUT_TAR_FILE=$INPUT_TAR_FILE"                                                        1>> ${LOG} 2>&1
          echo ">x\$INPUT_TAR_FILE={x$INPUT_TAR_FILE}"                                                      1>> ${LOG} 2>&1
          echo "Bool {x$INPUT_TAR_FILE != x}"                                                              1>>${LOG} 2>&1
          
if [ x$INPUT_TAR_FILE != x ]; then
     echo "CONDOR_SCRATCH_DIR=$CONDOR_SCRATCH_DIR"   1>>${LOG} 2>&1
     echo "TMP=$TMP"  1>>${LOG} 2>&1
     echo "TEMP=$TEMP"  1>>${LOG} 2>&1
     mkdir $CONDOR_SCRATCH_DIR/local 1>>${LOG} 2>&1
     cd $CONDOR_SCRATCH_DIR/local  1>>${LOG} 2>&1
     pwd                         1>>${LOG} 2>&1
     echo "Extracting TAR"      1>>${LOG} 2>&1
     cd ${TMP}/local
     tar -xf $INPUT_TAR_FILE   1>>${LOG} 2>&1

     echo "Extracted TAR"     1>>${LOG} 2>&1
     echo "Files extracted are:"     1>>${LOG} 2>&1
     ls                       1>>${LOG} 2>&1
     echo "Initializing localProducts from tarball ${INPUT_TAR_FILE}."     1>>${LOG} 2>&1
      sed "s@setenv MRB_INSTALL.*@setenv MRB_INSTALL ${TMP}/local@" $TMP/     local/setup | \
      sed "s@setenv MRB_TOP.*@setenv MRB_TOP ${TMP}@" > $TMP/local/setup.local


     echo "Setting up products"     1>>${LOG} 2>&1
     echo "ls of local dir: "       1>>${LOG} 2>&1
     #. ${TMP}/local/setup             1>>${LOG} 2>&1
     . ${TMP}/local/setup.local             1>>${LOG} 2>&1
     echo "mrbslp next"     1>>${LOG} 2>&1
     mrbslp 1>>${LOG} 2>&1
     echo "Setup tarbal done" 1>>${LOG} 2>&1
else
  echo "> ups list -aK+ $mrb_project"                                                            1>> ${LOG} 2>&1
  ups list -aK+ $mrb_project                                                                     1>> ${LOG} 2>&1
  echo "> setup $mrb_project $mrb_version -q$mrb_quals"                                          1>> ${LOG} 2>&1
  setup $mrb_project $mrb_version -q $mrb_quals                                                   1>> ${LOG} 2>&1
fi
          
          echo "Setting up environment log" 1>>${LOG} 2>&1
          ENVLOG=$CONDOR_DIR_LOG/environment_${label}.log
          touch $ENVLOG
          uname -a 1>> $LOG     2>&1
          uname -a 1>> $ENVLOG  2>&1
          cat /etc/redhat-release 1>> $LOG     2>&1
          cat /etc/redhat-release 1>> $ENVLOG  2>&1
          env >> $ENVLOG
          
          
          # Run the job
          echo ""                                                                 1>> ${LOG} 2>&1
          echo "***** Starting job"                                               1>> ${LOG} 2>&1
          echo "lar -c $FCL -n $NVoxelsPerJob -T ${RTF}" 1>> ${LOG} 2>&1
          lar -c $FCL -n $NVoxelsPerJob -T ${RTF}        1>> ${LOG} 2>&1
          ret=$?
          echo   "***** Job completed ($ret)"                                     1>> ${LOG} 2>&1
          echo                                                                    1>> ${LOG} 2>&1
          date                                                                    1>> ${LOG} 2>&1
          
exit $ret
