#!/bin/csh

if( $#argv < 1) then
    echo "must supply geometry base file name"
    exit 0
endif

set basename = $argv[1];

echo "testing geometry for $basename"

sed s/detector/$basename/ geotest.fcl > ${basename}test.fcl

lar -c ./${basename}test.fcl

echo "done testing geometry"

rm ${basename}test.fcl
rm *.log
