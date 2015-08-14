#!/usr/bin/perl
#
# Build a description of the detectors from gdml fragments
#

#
# List all the files that need to be glomed together
#


@defFile = qw(global-defs.gdml dune/dune-defs.gdml);
@gdmlFile = qw(dune/dune-rotations.gdml
	       materials.gdml
	       dune/dune-vertplane.gdml
	       dune/dune-plane.gdml
	       dune/dune-tpc.gdml
	       dune/dune-cryostat.gdml
	       dune/dune-enclosure.gdml
	       dune/dune-world.gdml
	      );


#
# Build the table of variable replacements
#

$ikey=0;
for ($ifile=0; $ifile<@defFile; ++$ifile) {
    
    $DEFINES = @defFile[$ifile];
    open(DEFINES) or die("Could not open file $DEFINES");
    foreach $line (<DEFINES>) {
	chomp($line);
	($key,$name,$value) = split(' ',$line);
	if ($key eq "<constant") {
	    $name  =~ tr/ //d;
	    $value =~ tr/ //d;
	    ($tmp,$name)  = split('=',$name);
	    ($tmp,$value) = split('=',$value);
	    ($tmp,$name, $tmp) = split('"',$name);
	    ($tmp,$value,$tmp) = split('"',$value);
	    $value =~ tr/\"//d;
	    @defNAME[$ikey]  = $name;
	    @defVALUE[$ikey] = "($value)";
	    print @defName[$ikey];
	    ++$ikey;
	}
    }
    close(DEFINES);
}

#
# Open each file in turn, read line-by-line making variable
# substitutions, and write to standard output
#

print "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n";
print "<gdml xmlns:gdml=\"http://cern.ch/2001/Schemas/GDML\"\n";
print "      xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n";
print "      xsi:noNamespaceSchemaLocation=\"GDMLSchema/gdml.xsd\">\n";
for ($i=0; $i<@gdmlFile; ++$i) 
{ 
    $FILE = @gdmlFile[$i];
    open(FILE) or die("Could not open file @gdmlFile[$i] for read.");
    foreach $line (<FILE>) {
	# 
	# Work backwards to keep dependency order correct
        #
	for ($j=@defNAME-1; $j>=0; --$j) {
	    $line =~ s/@defNAME[$j]/@defVALUE[$j]/g;
	}
	print $line;
    }
    close(FILE);	      
}
print "\n<setup name=\"Default\" version=\"1.0\">\n";
print "  <world ref=\"volWorld\" />\n";
print "</setup>\n\n";
print "</gdml>\n";
