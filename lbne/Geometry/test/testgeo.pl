#!/usr/bin/perl

if ( !$ARGV[0] ) {
    print "must supply geometry base file name\n";
    exit;
}

$basename=$ARGV[0];

print "testing geometry for $basename\n";

`sed s/detector/${basename}/ geotest.fcl > ${basename}test.fcl`;

$fclname="${basename}test.fcl";
system("lar -c $fclname");

print "done testing geometry\n";

`rm ${basename}test.fcl`;
