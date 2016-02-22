#!/bin/bash

outputDir=/lbne/data2/users/lbnedaq/nearline

# (optional): Remove nearline output files.
# find ${outputDir}/*/*/*/ -mtime +90 -name '*.DONE' -exec rm -v {} \;
# find ${outputDir}/*/*/*/ -mtime +90 -name '*.LOCK' -exec rm -v {} \;

# Delete stale output files...
find ${outputDir}/*/*/*/ -mtime +1 -name 'TFileService*.root' -exec rm -v {} \;
find /home/lbnedaq/nearline/temp/ -mtime +1 -name '*.root*' -exec rm -v {} \:

# Delete files from the EVD processing older than 12 hours...
find ${outputDir}_evd/*/*/*/ -mtime +0.5 -name '*.root' -exec rm -v {} \;
find ${outputDir}_evd/*/*/*/ -mtime +0.5 -name '*.LOCK' -exec rm -v {} \;
find ${outputDir}_evd/*/*/*/ -mtime +0.5 -name '*.txt' -exec rm -v {} \;
