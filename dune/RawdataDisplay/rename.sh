#!/bin/sh

# Script to change the AnalysisExample name to something else.
# Usage: ./rename.sh <new-name>

newname=$1

# Rename the files.
mv AnalysisExample_module.cc  ${newname}_module.cc
mv AnalysisExample.fcl        ${newname}.fcl

# Search through all the files, replacing the text "AnalysisExample"
# with the new name. 
sed -i -e "s/AnalysisExample/${newname}/g" ${newname}.fcl ${newname}_module.cc GNUmakefile

# Because of the way SRT handles compilations, it's necessary to
# create a link in $SRT_PRIVATE_CONTEXT/include for any package
# that you do not check out with addpkg_svn. 

# Get the name of the directory we're in (just in case it's not ${newname}).
fullpath=`pwd`
dir=`basename ${fullpath}`
# Create the link
ln -sf ../${dir} ../include/${dir}
