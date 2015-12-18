#!/bin/bash

outputDir=/data/nearline

# No need to keep done or lock files after 90 days...
find ${outputDir}/*/*/ -mtime +90 -name '*.DONE' -exec rm -v {} \;
find ${outputDir}/*/*/ -mtime +90 -name '*.LOCK' -exec rm -v {} \;

# Delete stale output files...
find ${outputDir}/*/*/ -mtime +1  -name 'TFileService*.root' -exec rm -v {} \;
