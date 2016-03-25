#!/usr/bin/perl

# Much of this program is taken straight from generate_gdml.pl that 
# generates MicroBooNE fragment files (Thank you.)

# Each subroutine generates a fragment GDML file, and the last subroutine
# creates an XML file that make_gdml.pl will use to appropriately arrange
# the fragment GDML files to create the final desired DUNE GDML file, 
# to be named by make_gdml output command

# If you are playing with different geometries, you can use the
# suffix command to help organize your work.

use vars;

use Math::Trig;
use Getopt::Long;
use Math::BigFloat;
Math::BigFloat->precision(-16);

GetOptions( "help|h" => \$help,
	    "suffix|s:s" => \$suffix,
	    "output|o:s" => \$output,
	    "wires|w:s" => \$wires,
            "helpcube|c" => \$helpcube);

if ( defined $help )
{
    # If the user requested help, print the usage notes and exit.
    usage();
    exit;
}

if ( ! defined $suffix )
{
    # The user didn't supply a suffix, so append nothing to the file
    # names.
    $suffix = "";
}
else
{
    # Otherwise, stick a "-" before the suffix, so that a suffix of
    # "test" applied to filename.gdml becomes "filename-test.gdml".
    $suffix = "-" . $suffix;
}




#++++++++++++++++++++++++ Begin defining variables +++++++++++++++++++++++++

# Define detector geometry variables - later to be put in a parameters
# XML file to be parsed as an input?

# set wires on to be the default, unless given an input by the user
$wires_on = 1; # 1=on, 0=off
if (defined $wires)
{
$wires_on = $wires
}

$tpc_on=1;
$inch = 2.54;



#################################################
####   4APA 35t parameters from DocDb 7550   ####
#################################################



##################################################################
##################### wire plane parameters ######################

$UWirePitch = .4888;
$VWirePitch = .5012;
$XWirePitch = .45;

$UAng[0] = 45.707;
$VAng[0] = 44.275;
$UAng[1] = 45.707;
$VAng[1] = 44.275;
$UAng[2] = 45.707;
$VAng[2] = 44.275;
$UAng[3] = 45.707;
$VAng[3] = 44.275;

#$UAng[0] = 36;
#$VAng[0] = 36;
#$UAng[3] = 36;
#$VAng[3] = 36;




##################################################################
######################## TPC parameters ##########################

$LongDrift               = 227;
$ShortDrift              = 27;
$APAFrame_x              = 2*$inch; #this does not include the wire spacing
$APAWirePlaneSpacing     = 0.476; # center to center spacing between all of the wire planes (g, u, v, and x)

$TPCWireThickness        = 0.015;
$TPCWirePlaneThickness   = $TPCWireThickness;

$APA_z = 50.175;


# Let APAs be numbered as follows
#  0 & 3 - Largest
#  1 - Middle
#  2 - Smallest

# Here are the APA heights indexed by APA number (lengths, z, are constant)

$APAHeight[0] = 196.0;
$APAHeight[1] = 112.0;
$APAHeight[2] = 84.0;
$APAHeight[3] = $APAHeight[0];


  # include APA spacing in y and z so volTPCs touch in y and z directions with correct APA
  # spacing - this makes for smoother event generation. 

$TPCSmallestLongDrift_x  = $LongDrift + 3*$APAWirePlaneSpacing + $TPCWirePlaneThickness;
$TPCSmallestLongDrift_y  = $APAHeight[2] + $APAVerticalGap;
$TPCSmallestLongDrift_z  = $APA_z + $APALongGap;

$TPCMidLongDrift_x       = $LongDrift + 3*$APAWirePlaneSpacing + $TPCWirePlaneThickness;
$TPCMidLongDrift_y       = $APAHeight[1] + $APAVerticalGap;
$TPCMidLongDrift_z       = $APA_z + $APALongGap;

$TPCLargestLongDrift_x   = $LongDrift + 3*$APAWirePlaneSpacing + $TPCWirePlaneThickness;
$TPCLargestLongDrift_y   = $APAHeight[0];
$TPCLargestLongDrift_z   = $APA_z + $APALongGap;

$TPCSmallestShortDrift_x = $ShortDrift + 3*$APAWirePlaneSpacing + $TPCWirePlaneThickness;
$TPCSmallestShortDrift_y = $APAHeight[2] + $APAVerticalGap;
$TPCSmallestShortDrift_z = $APA_z + $APALongGap;

$TPCMidShortDrift_x      = $ShortDrift + 3*$APAWirePlaneSpacing + $TPCWirePlaneThickness;
$TPCMidShortDrift_y      = $APAHeight[1] + $APAVerticalGap;
$TPCMidShortDrift_z      = $APA_z + $APALongGap;

$TPCLargestShortDrift_x  = $ShortDrift + 3*$APAWirePlaneSpacing + $TPCWirePlaneThickness;
$TPCLargestShortDrift_y  = $APAHeight[0];
$TPCLargestShortDrift_z  = $APA_z + $APALongGap;



############################################################
############### Optical Detector parameters ################

# TODO: while the structure is exactly what we need, the parameters for
# paddle height and positioning are all placeholders

$nPaddlesInAPA[0] = 3;
$nPaddlesInAPA[1] = 1;
$nPaddlesInAPA[2] = 1;
$nPaddlesInAPA[3] = $nPaddlesInAPA[0]; # for now, this one is the same as 0

$SiPM_y = 0;
$LightPaddle_x = 0.476;
$LightPaddle_y = 56; #  in cm from docDb 7803
$LightPaddle_z = 4*$inch;

# z and x are given by APA frame center.
#  Hardcode y distance of each paddle from the  
#  bottom of the APA to the paddle y-center.
# To be used in make_APA like [apa#][paddle#]

$PaddleYPositions[0][0] = $APAHeight[0]/2; # this puts it in the y center 
$PaddleYPositions[1][0] = $APAHeight[1] - 4*$inch - $LightPaddle_y/2; 
$PaddleYPositions[2][0] = $APAHeight[2] - 4*$inch - $LightPaddle_y/2;
$PaddleYPositions[3][0] = $PaddleYPositions[0][0];

$PaddleYPositions[0][1] = $APAHeight[0] - 4*$inch - $LightPaddle_y/2;
$PaddleYPositions[3][1] = $PaddleYPositions[0][1];

$PaddleYPositions[0][2] = 4*$inch + $LightPaddle_y/2;
$PaddleYPositions[3][2] = $PaddleYPositions[0][2];



##################################################################
###################### Cryostat parameters #######################

# "Dead" LAr outside of the fiducial volume but inside the cryostat
$SpaceCPAToCryoWall    =    10; #$SideLArPadding
$SpaceAPAToFloor       =    10; #$BottomLArPadding
$SpaceAPAToTopLAr      =    3.53; #$TopLArPadding
$HeightGaseousAr       =    13.5; #$Height of Gaseous Ar region
$UpstreamLArPadding    =    10;
$DownstreamLArPadding  =    10;

$APALongGap            =    1.5; #separation between APAs along the incident beam axis
$APAVerticalGap        =    2.5; #separation between APAs along the vertical axis

$CPA_x                 =    5.08;
$CPA_y                 =  ( $APAHeight[1] 
                          + $APAVerticalGap 
                          + $APAHeight[2] ) / 2;
$CPA_z                 =    77;    

# Liquid and Gaseous Argon dimensions 
$Argon_x	       =    $TPCLargestShortDrift_x 
                          + $APAFrame_x 
                          + $TPCLargestLongDrift_x 
                          + 2*$CPA_x
                          + 2*$SpaceCPAToCryoWall;
$Argon_y	       =    $APAHeight[1] 
                          + $APAVerticalGap 
                          + $APAHeight[2]
                          + $SpaceAPAToFloor 
                          + $SpaceAPAToTopLAr
                          + $HeightGaseousAr;  # assuming mid_y+smallest_y > largest_y
$Argon_z	       =    3*$APA_z + 2*$APALongGap
                          + $UpstreamLArPadding 
                          + $DownstreamLArPadding;

# Cryostat Dimensions
$SteelThickness	       =    0.5*$inch;
$Cryostat_x	       =    $Argon_x+2*$SteelThickness;
$Cryostat_y	       =    $Argon_y+2*$SteelThickness;
$Cryostat_z	       =    $Argon_z+2*$SteelThickness;



##################################################################
################# Detector Enclosure parameters ##################

$ConcretePadding       =    50;
$FoamPadding           =    80;
$TotalPadding	       =    $ConcretePadding+$FoamPadding;
$DetEnc_x	       =    $Cryostat_x+2*$TotalPadding;
$DetEnc_y	       =    $Cryostat_y+2*$TotalPadding-$FoamPadding; # no foam on bottom
$CryoInDetEnc_ypos     =    -$DetEnc_y/2 + $ConcretePadding + $Cryostat_y/2;
$DetEnc_z              =    $Cryostat_z+2*$TotalPadding;

  # We want the world origin to be at the very front of the fiducial volume.
  # move it to the front of the enclosure, then back it up through the concrete/foam, 
  # then through the Cryostat shell, then through the upstream dead LAr (including the
  # dead LAr on the edge of the TPC, but this is covered in $UpstreamLArPadding).
  # This is to be added to the z position of every volume in volWorld

$OriginZSet             =       $DetEnc_z/2 
                              - $TotalPadding 
                              - $SteelThickness 
                              - $UpstreamLArPadding;

  # We want the world origin to be vertically centered between the two stacked APAs.
  # (for now, that is, so the sorting works. this is quite asymetric, but then again
  #  so is the entire 35t geometry. this may be kept.)
  # the cryostat sits on top of concrete padding, move the detector enclosure back
  # and then move the world origin to the bottom of the smallest/lowest TPC, then
  # and then up through the TPC, then back up to being centered between the stacked APAs.
  # This is to be added to the y/x position of every volume in volWorld

$OriginYSet             =       $DetEnc_y/2 
                              - $ConcretePadding
                              - $SteelThickness
                              - $SpaceAPAToFloor
                              - $TPCSmallestShortDrift_y
                              - $APAVerticalGap/2;



$OriginXSet             =     - $DetEnc_x/2
                              + $TotalPadding
                              + $SteelThickness
                              + $SpaceCPAToCryoWall
                              + $CPA_x
                              + $LongDrift;


$PosDirCubeSide = 0;
if (defined $helpcube)
{
$PosDirCubeSide = $ArToAr; #seems to be a good proportion
}


$World_x = 2*$DetEnc_x;
$World_y = 2*$DetEnc_y;
$World_z = 2*$DetEnc_z;




##################################################################
######################### CPA positions ##########################

$posCPA0_x  =  - $Argon_x/2 + $SpaceCPAToCryoWall   + $CPA_x/2;
$posCPA0_y  =  - $Argon_y/2 + $SpaceAPAToFloor      + $CPA_y/2;
$posCPA0_z  =  - $Argon_z/2 + $UpstreamLArPadding   + $CPA_z/2;  


$posCPA1_x  =    $Argon_x/2 - $SpaceCPAToCryoWall   - $CPA_x/2;
$posCPA1_y  =  - $Argon_y/2 + $SpaceAPAToFloor      + $CPA_y/2;
$posCPA1_z  =  - $Argon_z/2 + $UpstreamLArPadding   + $CPA_z/2;


$posCPA2_x  =  - $Argon_x/2 + $SpaceCPAToCryoWall   + $CPA_x/2;
$posCPA2_y  =    $Argon_y/2 - $SpaceAPAToTopLAr   
                            - $HeightGaseousAr      - $CPA_y/2;
$posCPA2_z  =  - $Argon_z/2 + $UpstreamLArPadding   + $CPA_z/2;


$posCPA3_x  =    $Argon_x/2 - $SpaceCPAToCryoWall   - $CPA_x/2;
$posCPA3_y  =    $Argon_y/2 - $SpaceAPAToTopLAr   
                            - $HeightGaseousAr      - $CPA_y/2;
$posCPA3_z  =  - $Argon_z/2 + $UpstreamLArPadding   + $CPA_z/2;


$posCPA4_x  =  - $Argon_x/2 + $SpaceCPAToCryoWall   + $CPA_x/2;
$posCPA4_y  =  - $Argon_y/2 + $SpaceAPAToFloor      + $CPA_y/2;
$posCPA4_z  =    $Argon_z/2 - $DownstreamLArPadding - $CPA_z/2;


$posCPA5_x  =    $Argon_x/2 - $SpaceCPAToCryoWall   - $CPA_x/2;
$posCPA5_y  =  - $Argon_y/2 + $SpaceAPAToFloor      + $CPA_y/2;
$posCPA5_z  =    $Argon_z/2 - $DownstreamLArPadding - $CPA_z/2;


$posCPA6_x  =  - $Argon_x/2 + $SpaceCPAToCryoWall   + $CPA_x/2;
$posCPA6_y  =    $Argon_y/2 - $SpaceAPAToTopLAr   
                            - $HeightGaseousAr      - $CPA_y/2;
$posCPA6_z  =    $Argon_z/2 - $DownstreamLArPadding - $CPA_z/2;


$posCPA7_x  =    $Argon_x/2 - $SpaceCPAToCryoWall   - $CPA_x/2;
$posCPA7_y  =    $Argon_y/2 - $SpaceAPAToTopLAr   
                            - $HeightGaseousAr      - $CPA_y/2;
$posCPA7_z  =    $Argon_z/2 - $DownstreamLArPadding - $CPA_z/2;




##################################################################
######################### TPC positions ##########################



# Largest Short Drift
$posTPClsd_x  =  - $Argon_x/2 
                 + $SpaceCPAToCryoWall
                 + $CPA_x
                 + $TPCLargestShortDrift_x/2;
$posTPClsd_y  =    $Argon_y/2
                 - $HeightGaseousAr
                 - $SpaceAPAToTopLAr
                 - $TPCLargestShortDrift_y/2;

# Largest Long Drift
$posTPClld_x  =    $Argon_x/2 
                 - $SpaceCPAToCryoWall
                 - $CPA_x
                 - $TPCLargestLongDrift_x/2;
$posTPClld_y  =    $Argon_y/2
                 - $HeightGaseousAr
                 - $SpaceAPAToTopLAr
                 - $TPCLargestLongDrift_y/2;

$posTPCl_z_upstream    =    - $Argon_z/2
                            + $UpstreamLArPadding ## this steps into the APALongGap/2 space on outside of TPC...
                            + $APA_z/2;
$posTPCl_z_downstream  =      $Argon_z/2
                            - $DownstreamLArPadding ## ...as well as this
                            - $APA_z/2;



# Smallest Long Drift
$posTPCsld_x  =    $Argon_x/2 
                 - $SpaceCPAToCryoWall
                 - $CPA_x
                 - $TPCSmallestLongDrift_x/2; 
$posTPCsld_y  =  - $Argon_y/2
                 + $SpaceAPAToFloor
                 + $APAHeight[2]/2;  
$posTPCsld_z  =    $Argon_z/2
                 - $DownstreamLArPadding
                 - $APA_z
                 - $APALongGap
                 - $APA_z/2;

# Mid Long Drift
$posTPCmld_x  =    $Argon_x/2 
                 - $SpaceCPAToCryoWall
                 - $CPA_x
                 - $TPCMidLongDrift_x/2; 
$posTPCmld_y  =    $Argon_y/2
                 - $HeightGaseousAr
                 - $SpaceAPAToTopLAr
                 - $APAHeight[1]/2;
$posTPCmld_z  =    $Argon_z/2
                 - $DownstreamLArPadding
                 - $APA_z
                 - $APALongGap
                 - $APA_z/2;



# Smallest Short Drift
$posTPCssd_x  =  - $Argon_x/2 
                 + $SpaceCPAToCryoWall
                 + $CPA_x
                 + $TPCSmallestShortDrift_x/2; 
$posTPCssd_y  =  - $Argon_y/2
                 + $SpaceAPAToFloor
                 + $APAHeight[2]/2; 
$posTPCssd_z  =    $Argon_z/2
                 - $DownstreamLArPadding
                 - $APA_z
                 - $APALongGap
                 - $APA_z/2; 

# Mid Short Drift
$posTPCmsd_x  =  - $Argon_x/2 
                 + $SpaceCPAToCryoWall
                 + $CPA_x
                 + $TPCMidShortDrift_x/2; 
$posTPCmsd_y  =    $Argon_y/2
                 - $HeightGaseousAr
                 - $SpaceAPAToTopLAr
                 - $APAHeight[1]/2; 
$posTPCmsd_z  =    $Argon_z/2
                 - $DownstreamLArPadding
                 - $APA_z
                 - $APALongGap
                 - $APA_z/2; 



#TODO: define APA centers more fundamentally and then place TPC centers based off of this
$APA_Xcenter     =   $Argon_x/2 
                   - $SpaceCPAToCryoWall
                   - $CPA_x
                   - $TPCLargestLongDrift_x
                   - $APAFrame_x/2;

$APACenter[0][0] = $APA_Xcenter;
$APACenter[0][1] = $posTPClld_y;
$APACenter[0][2] = $posTPCl_z_upstream;

$APACenter[1][0] = $APA_Xcenter;
$APACenter[1][1] = $posTPCmld_y;
$APACenter[1][2] = $posTPCmld_z;

$APACenter[2][0] = $APA_Xcenter;
$APACenter[2][1] = $posTPCsld_y;
$APACenter[2][2] = $posTPCsld_z;

$APACenter[3][0] = $APACenter[0][0];
$APACenter[3][1] = $APACenter[0][1];
$APACenter[3][2] = $posTPCl_z_downstream;


################

#+++++++++++++++++++++++++ End defining variables ++++++++++++++++++++++++++



# run the sub routines that generate the fragments

gen_Define(); 	 # generates definitions at beginning of GDML
gen_Materials(); # generates materials to be used

     # pass a name to gen_TPC that begins with TPC
    gen_TPC( $TPCSmallestLongDrift_x, $TPCSmallestLongDrift_y, $TPCSmallestLongDrift_z, 'SmallestLongDrift', 1);
    gen_TPC( $TPCMidLongDrift_x, $TPCMidLongDrift_y, $TPCMidLongDrift_z, 'MidLongDrift', 2);
    gen_TPC( $TPCLargestLongDrift_x, $TPCLargestLongDrift_y, $TPCLargestLongDrift_z, 'LargestLongDrift', 0);

    gen_TPC( $TPCSmallestShortDrift_x, $TPCSmallestShortDrift_y, $TPCSmallestShortDrift_z, 'SmallestShortDrift', 1);
    gen_TPC( $TPCMidShortDrift_x, $TPCMidShortDrift_y, $TPCMidShortDrift_z, 'MidShortDrift', 2) ;
    gen_TPC( $TPCLargestShortDrift_x, $TPCLargestShortDrift_y, $TPCLargestShortDrift_z, 'LargestShortDrift', 0);


gen_Cryostat();
gen_Enclosure();
gen_World();	 


write_fragments(); # writes the XML input for make_gdml.pl
			# which zips together the final GDML
exit;


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++ usage +++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub usage()
{
    print "Usage: $0 [-h|--help] [-o|--output <fragments-file>] [-s|--suffix <string>]\n";
    print "       if -o is omitted, output goes to STDOUT; <fragments-file> is input to make_gdml.pl\n";
    print "       -s <string> appends the string to the file names; useful for multiple detector versions\n";
    print "       -h prints this message, then quits\n";
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_Define +++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Define()
{

# Create the <define> fragment file name, 
# add file to list of fragments,
# and open it
    $DEF = "dune_10kT_Def" . $suffix . ".gdml";
    push (@gdmlFiles, $DEF);
    $DEF = ">" . $DEF;
    open(DEF) or die("Could not open file $DEF for writing");


print DEF <<EOF;
<?xml version='1.0'?>
<gdml>
<define>

<!-- 


-->

   <position name="posOriginSet"        unit="cm" x="$OriginXSet" y="$OriginYSet" z="$OriginZSet"/>

   <position name="posTPCLargestShortDrift_Pos"  unit="cm" x="$posTPClsd_x" y="$posTPClsd_y" z="$posTPCl_z_downstream"/>
   <position name="posTPCLargestLongDrift_Pos"   unit="cm" x="$posTPClld_x" y="$posTPClld_y" z="$posTPCl_z_downstream"/>
   <position name="posTPCLargestShortDrift_Neg"  unit="cm" x="$posTPClsd_x" y="$posTPClsd_y" z="$posTPCl_z_upstream"/>
   <position name="posTPCLargestLongDrift_Neg"   unit="cm" x="$posTPClld_x" y="$posTPClld_y" z="$posTPCl_z_upstream"/>
   <position name="posTPCSmallestShortDrift"     unit="cm" x="$posTPCssd_x" y="$posTPCssd_y" z="$posTPCssd_z"/>
   <position name="posTPCSmallestLongDrift"      unit="cm" x="$posTPCsld_x" y="$posTPCsld_y" z="$posTPCsld_z"/>
   <position name="posTPCMidShortDrift"          unit="cm" x="$posTPCmsd_x" y="$posTPCmsd_y" z="$posTPCmsd_z"/>
   <position name="posTPCMidLongDrift"           unit="cm" x="$posTPCmld_x" y="$posTPCmld_y" z="$posTPCmld_z"/>


   <position name="posCathode0" unit="cm" x="$posCPA0_x" y="$posCPA0_y" z="$posCPA0_z"/>
   <position name="posCathode1" unit="cm" x="$posCPA1_x" y="$posCPA1_y" z="$posCPA1_z"/>
   <position name="posCathode2" unit="cm" x="$posCPA2_x" y="$posCPA2_y" z="$posCPA2_z"/>
   <position name="posCathode3" unit="cm" x="$posCPA3_x" y="$posCPA3_y" z="$posCPA3_z"/>
   <position name="posCathode4" unit="cm" x="$posCPA4_x" y="$posCPA4_y" z="$posCPA4_z"/>
   <position name="posCathode5" unit="cm" x="$posCPA5_x" y="$posCPA5_y" z="$posCPA5_z"/>
   <position name="posCathode6" unit="cm" x="$posCPA6_x" y="$posCPA6_y" z="$posCPA6_z"/>
   <position name="posCathode7" unit="cm" x="$posCPA7_x" y="$posCPA7_y" z="$posCPA7_z"/>

   <position name="posCenter"           unit="cm" x="0" y="0" z="0"/>
   <rotation name="rPlus90AboutX"       unit="deg" x="90" y="0" z="0"/>
   <rotation name="rMinus90AboutY"      unit="deg" x="0" y="270" z="0"/>
   <rotation name="rMinus90AboutYMinus90AboutX"       unit="deg" x="270" y="270" z="0"/>
   <rotation name="rPlus180AboutY"	unit="deg" x="0" y="180"   z="0"/>
   <rotation name="rIdentity"		unit="deg" x="0" y="0"   z="0"/>
</define>
</gdml>
EOF
    close (DEF);
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++ gen_Materials +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Materials() 
{

# Create the <materials> fragment file name,
# add file to list of output GDML fragments,
# and open it
    $MAT = "dune_10kT_Materials" . $suffix . ".gdml";
    push (@gdmlFiles, $MAT);
    $MAT = ">" . $MAT;
    open(MAT) or die("Could not open file $MAT for writing");


  print MAT <<EOF;
<materials>
  <element name="videRef" formula="VACUUM" Z="1">  <atom value="1"/> </element>
  <element name="bromine" formula="Br" Z="35"> <atom value="79.904"/> </element>
  <element name="hydrogen" formula="H" Z="1">  <atom value="1.0079"/> </element>
  <element name="nitrogen" formula="N" Z="7">  <atom value="14.0067"/> </element>
  <element name="oxygen" formula="O" Z="8">  <atom value="15.999"/> </element>
  <element name="aluminum" formula="Al" Z="13"> <atom value="26.9815"/>  </element>
  <element name="silicon" formula="Si" Z="14"> <atom value="28.0855"/>  </element>
  <element name="carbon" formula="C" Z="6">  <atom value="12.0107"/>  </element>
  <element name="potassium" formula="K" Z="19"> <atom value="39.0983"/>  </element>
  <element name="chromium" formula="Cr" Z="24"> <atom value="51.9961"/>  </element>
  <element name="iron" formula="Fe" Z="26"> <atom value="55.8450"/>  </element>
  <element name="nickel" formula="Ni" Z="28"> <atom value="58.6934"/>  </element>
  <element name="calcium" formula="Ca" Z="20"> <atom value="40.078"/>   </element>
  <element name="magnesium" formula="Mg" Z="12"> <atom value="24.305"/>   </element>
  <element name="sodium" formula="Na" Z="11"> <atom value="22.99"/>    </element>
  <element name="titanium" formula="Ti" Z="22"> <atom value="47.867"/>   </element>
  <element name="argon" formula="Ar" Z="18"> <atom value="39.9480"/>  </element>
  <element name="sulphur" formula="S" Z="16"> <atom value="32.065"/>  </element>
  <element name="phosphorus" formula="P" Z="15"> <atom value="30.973"/>  </element>

  <material name="Vacuum" formula="Vacuum">
   <D value="1.e-25" unit="g/cm3"/>
   <fraction n="1.0" ref="videRef"/>
  </material>

  <material name="ALUMINUM_Al" formula="ALUMINUM_Al">
   <D value="2.6990" unit="g/cm3"/>
   <fraction n="1.0000" ref="aluminum"/>
  </material>

  <material name="SILICON_Si" formula="SILICON_Si">
   <D value="2.3300" unit="g/cm3"/>
   <fraction n="1.0000" ref="silicon"/>
  </material>

  <material name="epoxy_resin" formula="C38H40O6Br4">
   <D value="1.1250" unit="g/cm3"/>
   <composite n="38" ref="carbon"/>
   <composite n="40" ref="hydrogen"/>
   <composite n="6" ref="oxygen"/>
   <composite n="4" ref="bromine"/>
  </material>

  <material name="SiO2" formula="SiO2">
   <D value="2.2" unit="g/cm3"/>
   <composite n="1" ref="silicon"/>
   <composite n="2" ref="oxygen"/>
  </material>

  <material name="Al2O3" formula="Al2O3">
   <D value="3.97" unit="g/cm3"/>
   <composite n="2" ref="aluminum"/>
   <composite n="3" ref="oxygen"/>
  </material>

  <material name="Fe2O3" formula="Fe2O3">
   <D value="5.24" unit="g/cm3"/>
   <composite n="2" ref="iron"/>
   <composite n="3" ref="oxygen"/>
  </material>

  <material name="CaO" formula="CaO">
   <D value="3.35" unit="g/cm3"/>
   <composite n="1" ref="calcium"/>
   <composite n="1" ref="oxygen"/>
  </material>

  <material name="MgO" formula="MgO">
   <D value="3.58" unit="g/cm3"/>
   <composite n="1" ref="magnesium"/>
   <composite n="1" ref="oxygen"/>
  </material>

  <material name="Na2O" formula="Na2O">
   <D value="2.27" unit="g/cm3"/>
   <composite n="2" ref="sodium"/>
   <composite n="1" ref="oxygen"/>
  </material>

  <material name="TiO2" formula="TiO2">
   <D value="4.23" unit="g/cm3"/>
   <composite n="1" ref="titanium"/>
   <composite n="2" ref="oxygen"/>
  </material>

  <material name="FeO" formula="FeO">
   <D value="5.745" unit="g/cm3"/>
   <composite n="1" ref="iron"/>
   <composite n="1" ref="oxygen"/>
  </material>

  <material name="CO2" formula="CO2">
   <D value="1.562" unit="g/cm3"/>
   <composite n="1" ref="iron"/>
   <composite n="2" ref="oxygen"/>
  </material>

  <material name="P2O5" formula="P2O5">
   <D value="1.562" unit="g/cm3"/>
   <composite n="2" ref="phosphorus"/>
   <composite n="5" ref="oxygen"/>
  </material>

  <material formula=" " name="DUSEL_Rock">
    <D value="2.82" unit="g/cm3"/>
    <fraction n="0.5267" ref="SiO2"/>
    <fraction n="0.1174" ref="FeO"/>
    <fraction n="0.1025" ref="Al2O3"/>
    <fraction n="0.0473" ref="MgO"/>
    <fraction n="0.0422" ref="CO2"/>
    <fraction n="0.0382" ref="CaO"/>
    <fraction n="0.0240" ref="carbon"/>
    <fraction n="0.0186" ref="sulphur"/>
    <fraction n="0.0053" ref="Na2O"/>
    <fraction n="0.00070" ref="P2O5"/>
    <fraction n="0.0771" ref="oxygen"/>
  </material> 

  <material name="fibrous_glass">
   <D value="2.74351" unit="g/cm3"/>
   <fraction n="0.600" ref="SiO2"/>
   <fraction n="0.118" ref="Al2O3"/>
   <fraction n="0.001" ref="Fe2O3"/>
   <fraction n="0.224" ref="CaO"/>
   <fraction n="0.034" ref="MgO"/>
   <fraction n="0.010" ref="Na2O"/>
   <fraction n="0.013" ref="TiO2"/>
  </material>

  <material name="FR4">
   <D value="1.98281" unit="g/cm3"/>
   <fraction n="0.47" ref="epoxy_resin"/>
   <fraction n="0.53" ref="fibrous_glass"/>
  </material>

  <material name="STEEL_STAINLESS_Fe7Cr2Ni" formula="STEEL_STAINLESS_Fe7Cr2Ni">
   <D value="7.9300" unit="g/cm3"/>
   <fraction n="0.0010" ref="carbon"/>
   <fraction n="0.1792" ref="chromium"/>
   <fraction n="0.7298" ref="iron"/>
   <fraction n="0.0900" ref="nickel"/>
  </material>

  <material name="LAr" formula="LAr">
   <D value="1.40" unit="g/cm3"/>
   <fraction n="1.0000" ref="argon"/>
  </material>

  <material name="ArGas" formula="ArGas">
   <D value="0.00166" unit="g/cm3"/>
   <fraction n="1.0" ref="argon"/>
  </material>

  <material formula=" " name="Air">
   <D value="0.001205" unit="g/cm3"/>
   <fraction n="0.781154" ref="nitrogen"/>
   <fraction n="0.209476" ref="oxygen"/>
   <fraction n="0.00934" ref="argon"/>
  </material>

  <material formula=" " name="G10">
   <D value="1.7" unit="g/cm3"/>
   <fraction n="0.2805" ref="silicon"/>
   <fraction n="0.3954" ref="oxygen"/>
   <fraction n="0.2990" ref="carbon"/>
   <fraction n="0.0251" ref="hydrogen"/>
  </material>

  <material formula=" " name="Granite">
   <D value="2.7" unit="g/cm3"/>
   <fraction n="0.438" ref="oxygen"/>
   <fraction n="0.257" ref="silicon"/>
   <fraction n="0.222" ref="sodium"/>
   <fraction n="0.049" ref="aluminum"/>
   <fraction n="0.019" ref="iron"/>
   <fraction n="0.015" ref="potassium"/>
  </material>

  <material formula=" " name="ShotRock">
   <D value="2.7*0.6" unit="g/cm3"/>
   <fraction n="0.438" ref="oxygen"/>
   <fraction n="0.257" ref="silicon"/>
   <fraction n="0.222" ref="sodium"/>
   <fraction n="0.049" ref="aluminum"/>
   <fraction n="0.019" ref="iron"/>
   <fraction n="0.015" ref="potassium"/>
  </material>

  <material formula=" " name="Dirt">
   <D value="1.7" unit="g/cm3"/>
   <fraction n="0.438" ref="oxygen"/>
   <fraction n="0.257" ref="silicon"/>
   <fraction n="0.222" ref="sodium"/>
   <fraction n="0.049" ref="aluminum"/>
   <fraction n="0.019" ref="iron"/>
   <fraction n="0.015" ref="potassium"/>
  </material>

  <material formula=" " name="Concrete">
   <D value="2.3" unit="g/cm3"/>
   <fraction n="0.530" ref="oxygen"/>
   <fraction n="0.335" ref="silicon"/>
   <fraction n="0.060" ref="calcium"/>
   <fraction n="0.015" ref="sodium"/>
   <fraction n="0.020" ref="iron"/>
   <fraction n="0.040" ref="aluminum"/>
  </material>

  <material formula="H2O" name="Water">
   <D value="1.0" unit="g/cm3"/>
   <fraction n="0.1119" ref="hydrogen"/>
   <fraction n="0.8881" ref="oxygen"/>
  </material>

  <material formula="Ti" name="Titanium">
   <D value="4.506" unit="g/cm3"/>
   <fraction n="1." ref="titanium"/>
  </material>

  <material name="TPB" formula="TPB">
   <D value="1.40" unit="g/cm3"/>
   <fraction n="1.0000" ref="argon"/>
  </material>

  <material name="Glass">
   <D value="2.74351" unit="g/cm3"/>
   <fraction n="0.600" ref="SiO2"/>
   <fraction n="0.118" ref="Al2O3"/>
   <fraction n="0.001" ref="Fe2O3"/>
   <fraction n="0.224" ref="CaO"/>
   <fraction n="0.034" ref="MgO"/>
   <fraction n="0.010" ref="Na2O"/>
   <fraction n="0.013" ref="TiO2"/>
  </material>

  <material name="Acrylic">
   <D value="1.19" unit="g/cm3"/>
   <fraction n="0.600" ref="carbon"/>
   <fraction n="0.320" ref="oxygen"/>
   <fraction n="0.080" ref="hydrogen"/>
  </material>

</materials>
EOF

close(MAT);
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++ gen_TPC ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


sub gen_TPC()
{

# $_[0] = $TPC_x
# $_[1] = $TPC_y
# $_[2] = $TPC_z
# $_[3] = 'name'

    my $TPCWirePlaneHeight = $_[1];
    my $TPCWirePlaneLength = $_[2];

    my $TPCActive_x   =  $_[0]-(4*$APAWirePlaneSpacing);
    my $TPCActive_y   =  $TPCWirePlaneHeight;
    my $TPCActive_z   =  $TPCWirePlaneLength; 

    my $UAngle = $UAng[$_[4]];
    my $VAngle = $VAng[$_[4]];

    my $SinUAngle = sin( deg2rad($UAngle) );
    my $CosUAngle = cos( deg2rad($UAngle) );
    my $TanUAngle = tan( deg2rad($UAngle) );
    
    my $SinVAngle = sin( deg2rad($VAngle) );
    my $CosVAngle = cos( deg2rad($VAngle) );
    my $TanVAngle = tan( deg2rad($VAngle) );
    
    my $UWire_yint = $UWirePitch/$SinUAngle;
    my $UWire_zint = $UWirePitch/$CosUAngle;
    
    my $VWire_yint = $VWirePitch/$SinVAngle;
    my $VWire_zint = $VWirePitch/$CosVAngle;

#constructs everything inside volTPC, namely
# (moving from left to right, or from +x to -x)
#  -volCPA
#  -volTPCPlaneU: with wires angled from vertical slightly different than in V
#  -volTPCPlaneV: with wires angled from vertical slightly differently than in U
#  -volTPCPlaneX: with vertical wires


# Create the TPC fragment file name,
# add file to list of output GDML fragments,
# and open it
    $TPC = "35t_TPC_${_[3]}" . $suffix . ".gdml";
    push (@gdmlFiles, $TPC);
    $TPC = ">" . $TPC;
    open(TPC) or die("Could not open file $TPC for writing");


# The standard XML prefix and starting the gdml
    print TPC <<EOF;
<?xml version='1.0'?>
<gdml>
EOF


# All the TPC solids save the wires.
print TPC <<EOF;
<solids>
    <box name="$_[3]" lunit="cm" 
      x="$_[0]" 
      y="$_[1]" 
      z="$_[2]"/>
    <box name="${_[3]}Plane" lunit="cm" 
      x="$TPCWirePlaneThickness" 
      y="$TPCWirePlaneHeight" 
      z="$TPCWirePlaneLength"/>
    <box name="${_[3]}Active" lunit="cm"
      x="$TPCActive_x"
      y="$TPCActive_y"
      z="$TPCActive_z"/>
EOF


#++++++++++++++++++++++++++++ Wire Solids ++++++++++++++++++++++++++++++

print TPC <<EOF;

    <tube name="${_[3]}WireVert"
      rmax="0.5*$TPCWireThickness"
      z="$TPCWirePlaneHeight"               
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

# Set number of wires to default to zero, when $wires_on = 0, for a low memory 
# version. But if $wires_on = 1, calculate the number of wires on each side of each
# plane to be used in the for loops

my $NumberCornerUWires = 0;
my $NumberSideUWires = 0;
my $NumberCommonUWires = 0;
my $NumberCornerVWires = 0;
my $NumberSideVWires = 0;
my $NumberCommonVWires = 0;
my $NumberVerticalWires = 0;

if ($wires_on == 1)
{
   # Number of wires in one corner
$NumberCornerUWires = int( $TPCWirePlaneLength/($UWirePitch/$CosUAngle) );

$NumberCornerVWires = int( $TPCWirePlaneLength/($VWirePitch/$CosVAngle) );

   # Total number of wires touching one vertical (longer) side
   # Note that the total number of wires per plane is this + another set of corner wires
$NumberSideUWires = int( $TPCWirePlaneHeight/($UWirePitch/$SinUAngle) );

$NumberSideVWires = int( $TPCWirePlaneHeight/($VWirePitch/$SinVAngle) );

   # Number of wires per side that aren't cut off by the corner
$NumberCommonUWires = $NumberSideUWires - $NumberCornerUWires;

$NumberCommonVWires = $NumberSideVWires - $NumberCornerVWires;

   # number of wires on the vertical plane
# TODO: 1 is a temporary safety!!
$NumberVerticalWires = int( ($TPCWirePlaneLength - $TPCWireThickness)/$XWirePitch );
}

# The corner wires for the U plane
if ($wires_on==1) 
{
    for ($i = 0; $i < $NumberCornerUWires; ++$i)
    {
  # Subtraction to avoid corners of wires overlapping the TPCPlane sides,
  # equal to 0.5*TCPWireThickness*($TanUAngle+1/$TanUAngle),
  # allowing for 30deg<UAngle

   print TPC <<EOF;
    <tube name="${_[3]}WireU$i"
      rmax="0.5*$TPCWireThickness"
      z="$UWirePitch*($TanUAngle+1/$TanUAngle)*($i+1)-0.01732"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

    }
  # Next, the wire used many times in the middle of the U plane.
  # Subtraction again to avoid wire corners overlapping, equal to 
  # 0.5*TCPWireThickness*2/$TanVAngle, allowing for 30deg<VAngle

   print TPC <<EOF;
    <tube name="${_[3]}WireUCommon"
      rmax="0.5*$TPCWireThickness"
      z="$TPCWirePlaneLength/$SinUAngle-0.02598"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

}


# The corner wires for the V plane
if ($wires_on==1) 
{
    for ($i = 0; $i < $NumberCornerVWires; ++$i)
    {
    # Same subtraction to avoid corners of wires overlapping 
    # the TPCPlane sides

	print TPC <<EOF;

    <tube name="${_[3]}WireV$i"
      rmax="0.5*$TPCWireThickness"
      z="$VWirePitch*($TanVAngle+1/$TanVAngle)*($i+1)-0.01732"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>

EOF

    }

    # The wire used many times in the middle of the V plane
    # Same subtraction as U common 

   print TPC <<EOF;
    <tube name="${_[3]}WireVCommon"
      rmax="0.5*$TPCWireThickness"
      z="$TPCWirePlaneLength/$SinVAngle-0.02598"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

}



# Begin structure and create the vertical wire logical volume
print TPC <<EOF;
</solids>
<structure>
    <volume name="volTPCActive${_[3]}">
      <materialref ref="LAr"/>
      <solidref ref="${_[3]}Active"/>
    </volume>

EOF


if ($wires_on==1) 
{ 
  print TPC <<EOF;
    <volume name="volTPCWireVert${_[3]}">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="${_[3]}WireVert"/>
    </volume>
EOF

  # Corner U wires logical volumes
  for ($i = 0; $i < $NumberCornerUWires; ++$i)
  {
  print TPC <<EOF;
    <volume name="volTPCWireU$i${_[3]}">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="${_[3]}WireU$i"/>
    </volume>
EOF

  }

  # Common U wire logical volume, referenced many times
  print TPC <<EOF;
    <volume name="volTPCWireUCommon${_[3]}">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="${_[3]}WireUCommon"/>
    </volume>
EOF

  # Corner V wires logical volumes
  for ($i = 0; $i < $NumberCornerVWires; ++$i)
  {
  print TPC <<EOF;
    <volume name="volTPCWireV$i${_[3]}">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="${_[3]}WireV$i"/>
    </volume>
EOF

  }

  # Common V wire logical volume, referenced many times
  print TPC <<EOF;
    <volume name="volTPCWireVCommon${_[3]}">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="${_[3]}WireVCommon"/>
    </volume>
EOF

}




#+++++++++++++++++++++++++ Position physical wires ++++++++++++++++++++++++++

#            ++++++++++++++++++++++  U Plane  +++++++++++++++++++++++

# Create U plane logical volume
print TPC <<EOF;
    <volume name="volTPCPlaneU${_[3]}">
      <materialref ref="LAr"/>
      <solidref ref="${_[3]}Plane"/>
EOF

if ($wires_on==1)
{

# Starting with the bottom left corner wires:
   # x=0 to center the wires in the plane
   # y positioning: (-0.5*$TPCWirePlaneHeight) starts the incremental increase
        # from the bottom of the plane, and trigonometry gives the increment
   # z positioning: Looking at the plane from the positive x direction,
        # (0.5*$TPCWirePlaneLength) starts the incremental increase from
        # the lower left corner.
   # rotation: same as common wire in code below

for ($i = 0; $i < $NumberCornerUWires; ++$i)
{
my $ypos = (-0.5*$TPCWirePlaneHeight)+0.5*($i+1)*$UWire_yint;
my $zpos = (0.5*$TPCWirePlaneLength)-0.5*($i+1)*$UWire_zint;

my $diff=(0.5*$TPCWirePlaneLength)-0.5*($NumberCornerUWires)*$UWire_zint;
my $zpos=$zpos-$diff;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireU$i${_[3]}"/>
        <position name="pos${_[3]}WireU$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotation name="rUAngle$i"   unit="deg" x="90-$UAngle" y="0"   z="0"/>
      </physvol>
EOF

}


# Moving upwards to the common wires:
   # x and z are zero to center the wires along a vertical axis
   # y positioning: The trick is positioning the lowest common wire so that the pitch
        # is consistent, then the increment is double the increment of
        # the corner wires since there is no z incriment.
   # rotation: wires in  \\\\  direction, so +90deg to bring them to vertical and 
        # +UAngle counterclockwise to arrive at proper orientation
# Note that the counter maintains wire number (in pos. name) counting bottom to top

for ($i = $NumberCornerUWires; $i < $NumberSideUWires; ++$i)
{
my $ypos = (-0.5*$TPCWirePlaneHeight)+0.5*($NumberCornerUWires)*$UWire_yint+($i+1-$NumberCornerUWires)*$UWire_yint;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireUCommon${_[3]}"/>
        <position name="pos${_[3]}WireU$i" unit="cm" x="0" y="$ypos " z="0"/>
	<rotation name="rUAngle$i"   unit="deg" x="90-$UAngle" y="0"   z="0"/>
      </physvol>
EOF

}


# Finally moving to the corner wires on the top right:
   # x=0 to center the wires in the plane
   # y positioning: plug wire number into same equation
   # z positioning: start at z=0 and go negatively at the same z increment
   # rotation: same as common wire in code above
# note that the counter maintains wire number shown in the position name

for ($i = $NumberSideUWires; $i < $NumberSideUWires+$NumberCornerUWires-1; ++$i)
{
   # Make a counter to recall the right logical volume reference:
   # We want the last U wire in this loop (the highest wire) to be the 
   # first wire in the logical volume loop for U wires. 

$j = $NumberSideUWires+$NumberCornerUWires - $i - 2;

   # Note that since we are referencing the same logical volumes/same solids for
   # the top wires as well as the bottom, the pattern of "stacking" wire on top of wire
   # with an incremental separation is likely to cause the top corner wires to be a 
   # a little shorter than they can be, but never any longer. There is no immediately
   # elegant way to fix this, but at 5mm pitch and around 45deg wire orientation, the
   # wires can be at most 1cm shorter than possible which is negligible until the top
   # 20 wires or so where 1cm is >5% of their length. This also means that there
   # could be one more space for a wire left over, but that is highly unlikely.

my $ypos = (-0.5*$TPCWirePlaneHeight)+0.5*($NumberCornerUWires)*$UWire_yint+($NumberCommonUWires)*$UWire_yint+0.5*($i+1-$NumberSideUWires)*$UWire_yint;
my $zpos = -0.5*($i+1-$NumberSideUWires)*$UWire_zint;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireU$j${_[3]}"/>
        <position name="pos${_[3]}WireU$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotation name="rUAngle$i"   unit="deg" x="90-$UAngle" y="0"   z="0"/>
      </physvol>
EOF

}

} #ends if wires on


#            ++++++++++++++++++++++  V Plane  +++++++++++++++++++++++

# End U plane and create V plane logical volume
print TPC <<EOF;
    </volume>

    <volume name="volTPCPlaneV${_[3]}">
      <materialref ref="LAr"/>
      <solidref ref="${_[3]}Plane"/>
EOF

if ($wires_on==1)
{


# Starting with the bottom right corner wires:
   # x=0 to center the wires in the plane
   # y positioning: (-0.5*$TPCWirePlaneHeight) starts the incremental increase
        # from the bottom of the plane, and trigonometry gives the increment
   # z positioning: Looking at the plane from the positive x direction,
        # (-0.5*$TPCWirePlaneLength) starts the incremental increase from 
        # the lower right corner.
   # rotation: same as common wire in code below

for ($i = 0; $i < $NumberCornerVWires; ++$i)
{
my $ypos = (-0.5*$TPCWirePlaneHeight)+0.5*($i+1)*$VWire_yint;
my $zpos = (-0.5*$TPCWirePlaneLength)+0.5*($i+1)*$VWire_zint;

my $diff=(-0.5*$TPCWirePlaneLength)+0.5*($NumberCornerVWires)*$VWire_zint;
my $zpos=$zpos-$diff;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireV$i${_[3]}"/>
        <position name="pos${_[3]}WireV$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotation name="rVAngle$i"   unit="deg" x="90+$VAngle" y="0"   z="0"/>
      </physvol>
EOF

}


# Moving upwards to the common wires:
   # x and z are zero to center the wires along a vertical axis
   # y positioning: Plug wire number into the same corner ypos equation
   # rotation: wires in  ////  direction, so +90deg to bring them to vertical and 
        # --VAngle counterclockwise to arrive at proper orientation
# Note that the counter maintains wire number in the position name

for ($i = $NumberCornerVWires; $i < $NumberSideVWires; ++$i)
{
my $ypos = (-0.5*$TPCWirePlaneHeight)-0.5*($NumberCornerVWires)*$VWire_yint+($i+1)*$VWire_yint;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireVCommon${_[3]}"/>
        <position name="pos${_[3]}WireV$i" unit="cm" x="0" y="$ypos " z="0"/>
        <rotation name="rVAngle$i"   unit="deg" x="90+$VAngle" y="0"   z="0"/>
      </physvol>
EOF

}


# Finally moving to the corner wires on the top right:
   # x=0 to center the wires in the plane
   # y positioning: plug wire number into same equation
   # z positioning: start at z=0 and go positively at the same z increment
   # rotation: same as common wire in code above
# note that the counter maintains wire number shown in the position name

for ($i = $NumberSideVWires; $i < $NumberSideVWires+$NumberCornerVWires-1; ++$i)
{
   # Make a counter to recall the right logical volume reference where the last
   # wire in this loop is the smallest, first wire in the logical volume loop, just as in U

$j = $NumberSideVWires+$NumberCornerVWires - $i - 2;

   # Note that since we are referencing the same logical volumes/same solids for
   # the top wires as well as the bottom, the pattern of "stacking" wire on top of wire
   # with an incremental separation is likely to cause the top corner wires to be a 
   # a little shorter than they can be, but never any longer. Just as in U

my $ypos = (-0.5*$TPCWirePlaneHeight)+0.5*($NumberCornerVWires)*$VWire_yint+($NumberCommonVWires)*$VWire_yint+0.5*($i+1-$NumberSideVWires)*$VWire_yint;
my $zpos = 0.5*($i+1-$NumberSideVWires)*$VWire_zint;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireV$j${_[3]}"/>
        <position name="pos${_[3]}WireV$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotation name="rVAngle$i"   unit="deg" x="90+$VAngle" y="0"   z="0"/>
      </physvol>
EOF
}

} #ends if wires on



#            ++++++++++++++++++++++  X Plane  +++++++++++++++++++++++

# End V plane and create X plane logical volume
print TPC <<EOF;
    </volume>

    <volume name="volTPCPlaneX${_[3]}">
      <materialref ref="LAr"/>
      <solidref ref="${_[3]}Plane"/>
EOF

if ($wires_on==1)
{

# This is the simplest plane, one loop creates all of the wires
   # x and y position at zero to center the wires
   # z position: moving from front of detector to back, in the positive z direction,
        # starting at (-0.5*$TPCWirePlaneLength), the right side looking from 
        # the +x direction

for ($i=0; $i<$NumberVerticalWires; ++$i)
{
my $zpos = (-0.5*$TPCWirePlaneLength)+$XWirePitch*($i+0.5)+$TPCWireThickness/2;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireVert${_[3]}"/>
        <position name="pos${_[3]}WireX$i" unit="cm" x="0" y="0 " z="$zpos"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
EOF

}

} #ends if wires on

print TPC <<EOF;
    </volume>
EOF

#+++++++++++++++++++++ Position physical wires Above +++++++++++++++++++++

#wrap up the TPC file
print TPC <<EOF;
    <volume name="volTPC${_[3]}">
      <materialref ref="LAr"/>
      <solidref ref="${_[3]}"/>
     <physvol>
       <volumeref ref="volTPCPlaneU${_[3]}"/>
       <position name="pos${_[3]}PlaneU" unit="cm" x="-($_[0]/2)+3*$APAWirePlaneSpacing" y="0" z="0"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneV${_[3]}"/>
       <position name="pos${_[3]}PlaneV" unit="cm" x="-($_[0]/2)+2*$APAWirePlaneSpacing" y="0" z="0"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneX${_[3]}"/>
       <position name="pos${_[3]}PlaneX" unit="cm" x="-($_[0]/2)+$APAWirePlaneSpacing" y="0" z="0"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCActive${_[3]}"/>
       <position name="pos${_[3]}Active" unit="cm" x="-($_[0]/2)+$APAWirePlaneSpacing + $TPCActive_x/2" y="0" z="0"/>
     </physvol>
    </volume>
</structure>
</gdml>
EOF

    close(GDML);

} #end of sub gen_TPC





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_Cryostat +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Cryostat()
{

# Create the cryostat fragment file name,
# add file to list of output GDML fragments,
# and open it
    $CRYO = "35t_Cryostat" . $suffix . ".gdml";
    push (@gdmlFiles, $CRYO);
    $CRYO = ">" . $CRYO;
    open(CRYO) or die("Could not open file $CRYO for writing");


# The standard XML prefix and starting the gdml
    print CRYO <<EOF;
<?xml version='1.0'?>
<gdml>
EOF


# All the cryostat solids.
print CRYO <<EOF;
<solids>
    <box name="Cryostat" lunit="cm" 
      x="$Cryostat_x" 
      y="$Cryostat_y" 
      z="$Cryostat_z"/>

    <box name="ArgonInterior" lunit="cm" 
      x="$Argon_x"
      y="$Argon_y"
      z="$Argon_z"/>

    <subtraction name="SteelShell">
      <first ref="Cryostat"/>
      <second ref="ArgonInterior"/>
    </subtraction>

    <box name="Cathode" lunit="cm"
      x="$CPA_x"
      y="$CPA_y"
      z="$CPA_z"/>

    <box name="LightPaddle" lunit="cm"
      x="$LightPaddle_x"
      y="$LightPaddle_y + $SiPM_y"
      z="$LightPaddle_z"/>

EOF


    # Add APA frames and optical detectors (paddles)
    solid_APA( 0 );
    solid_APA( 1 );
    solid_APA( 2 );
    solid_APA( 3 );


# Cryostat structure
print CRYO <<EOF;
</solids>
<structure>
    <volume name="volSteelShell">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="SteelShell"/>
    </volume>

    <volume name="volCathode">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="Cathode"/>
    </volume>

    <volume name="volOpDetSensitive">
      <materialref ref="Acrylic"/>
      <solidref ref="LightPaddle"/>
    </volume>

EOF

    vol_APA( 0 );
    vol_APA( 1 );
    vol_APA( 2 );
    vol_APA( 3 );

print CRYO <<EOF;

    <volume name="volCryostat">
      <materialref ref="LAr"/>
      <solidref ref="Cryostat"/>

      <physvol>
        <volumeref ref="volSteelShell"/>
        <position name="posSteelShell" unit="cm" x="0" y="0" z="0"/> 
      </physvol>


     <physvol>
       <volumeref ref="volCathode"/>
       <positionref ref="posCathode0"/>
     </physvol>
     <physvol>
       <volumeref ref="volCathode"/>
       <positionref ref="posCathode1"/>
     </physvol>
     <physvol>
       <volumeref ref="volCathode"/>
       <positionref ref="posCathode2"/>
     </physvol>
     <physvol>
       <volumeref ref="volCathode"/>
       <positionref ref="posCathode3"/>
     </physvol>

     <physvol>
       <volumeref ref="volCathode"/>
       <positionref ref="posCathode4"/>
     </physvol>
     <physvol>
       <volumeref ref="volCathode"/>
       <positionref ref="posCathode5"/>
     </physvol>
     <physvol>
       <volumeref ref="volCathode"/>
       <positionref ref="posCathode6"/>
     </physvol>
     <physvol>
       <volumeref ref="volCathode"/>
       <positionref ref="posCathode7"/>
     </physvol>


     <!-- The Smallest APA -->
     <physvol>
        <volumeref ref="volTPCSmallestLongDrift"/>
        <positionref ref="posTPCSmallestLongDrift"/>
	<rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPCSmallestShortDrift"/>
        <positionref ref="posTPCSmallestShortDrift"/>
	<rotationref ref="rPlus180AboutY"/>
      </physvol>

EOF

place_APA($APACenter[2][0], $APACenter[2][1], $APACenter[2][2], 2);



print CRYO <<EOF;
      <!-- The Middle-sized APA -->
      <physvol>
        <volumeref ref="volTPCMidLongDrift"/>
        <positionref ref="posTPCMidLongDrift"/>
	<rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPCMidShortDrift"/>
        <positionref ref="posTPCMidShortDrift"/>
	<rotationref ref="rPlus180AboutY"/>
      </physvol>

EOF

place_APA($APACenter[1][0], $APACenter[1][1], $APACenter[1][2], 1);



print CRYO <<EOF;
      <!-- The Largest APAs, Upstream and Downstream -->
      <physvol>
        <volumeref ref="volTPCLargestLongDrift"/>
        <positionref ref="posTPCLargestLongDrift_Neg"/>
	<rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPCLargestShortDrift"/>
        <positionref ref="posTPCLargestShortDrift_Neg"/>
	<rotationref ref="rPlus180AboutY"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPCLargestLongDrift"/>
        <positionref ref="posTPCLargestLongDrift_Pos"/>
	<rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPCLargestShortDrift"/>
        <positionref ref="posTPCLargestShortDrift_Pos"/>
	<rotationref ref="rPlus180AboutY"/>
      </physvol>

EOF

place_APA($APACenter[0][0], $APACenter[0][1], $APACenter[0][2], 0);
place_APA($APACenter[3][0], $APACenter[3][1], $APACenter[3][2], 3);



print CRYO <<EOF;
    </volume>
</structure>
</gdml>
EOF

close(CRYO);
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ solid_APA ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Must be called only within gen_Cryostat(), 


# $_[0] = APA number
#   convention: 0 is Largest, 1 is Mid, 2 is Smallest

sub solid_APA()
{

    $APA_i = $_[0];

####################################################################
################# APA Frame and Paddle Dimensions ##################

    $APA_y = $APAHeight[$APA_i];

$APAFrameZSide_x = $APAFrame_x;
$APAFrameZSide_y = 4*$inch;
$APAFrameZSide_z = $APA_z;

$APAFrameYSide_x = $APAFrame_x;
$APAFrameYSide_y = $APA_y-2*$APAFrameZSide_y;
$APAFrameYSide_z = 4*$inch;

# Two outer Y supports will sandwich the light paddles
$APAFrameYOuterSupport_x = ($APAFrame_x-$LightPaddle_x)/2;
$APAFrameYOuterSupport_y = $APA_y-2*$APAFrameZSide_y;
$APAFrameYOuterSupport_z = 4*$inch;

$EdgeFrameSteelThickness = 0.12*$inch;
$InnerFrameSteelThickness = 0.062*$inch;


  print CRYO <<EOF;
     <box name="APAFrameYSideHollow\-$APA_i" lunit="cm"
      x="$APAFrameYSide_x-2*$EdgeFrameSteelThickness"
      y="$APAFrameYSide_y-2*$EdgeFrameSteelThickness"
      z="$APAFrameYSide_z"/>
     <box name="APAFrameYSideShell\-$APA_i" lunit="cm"
      x="$APAFrameYSide_x"
      y="$APAFrameYSide_y"
      z="$APAFrameYSide_z"/>
     <subtraction name="APAFrameYSide\-$APA_i">
      <first  ref="APAFrameYSideShell\-$APA_i"/>
      <second ref="APAFrameYSideHollow\-$APA_i"/>
      <positionref ref="posCenter"/>
      <rotationref ref="rIdentity"/>
      </subtraction>

     <box name="APAFrameZSideHollow\-$APA_i" lunit="cm"
      x="$APAFrameZSide_x-2*$EdgeFrameSteelThickness"
      y="$APAFrameZSide_y-2*$EdgeFrameSteelThickness"
      z="$APAFrameZSide_z"/>
     <box name="APAFrameZSideShell\-$APA_i" lunit="cm"
      x="$APAFrameZSide_x"
      y="$APAFrameZSide_y"
      z="$APAFrameZSide_z"/>
     <subtraction name="APAFrameZSide\-$APA_i">
      <first  ref="APAFrameZSideShell\-$APA_i"/>
      <second ref="APAFrameZSideHollow\-$APA_i"/>
      <positionref ref="posCenter"/>
      <rotationref ref="rIdentity"/>
      </subtraction>

     <box name="APAFrameYOuterSupport\-$APA_i" lunit="cm"
      x="$EdgeFrameSteelThickness"
      y="$APAFrameYOuterSupport_y"
      z="$APAFrameYOuterSupport_z"/>


EOF

}





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ vol_APA ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Must be called only within gen_Cryostat(), 


# $_[0] = x APA center
# $_[1] = y APA center
# $_[2] = z APA center
# $_[3] = APA number
#   convention: 0 and 3 are Largest, 1 is Mid, 2 is Smallest

sub vol_APA()
{

    $APA_i = $_[0];

print CRYO <<EOF;

    <volume name="volAPAFrameYSide\-$APA_i">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="APAFrameYSide\-$APA_i"/>
    </volume>

    <volume name="volAPAFrameZSide\-$APA_i">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="APAFrameZSide\-$APA_i"/>
    </volume>

    <volume name="volAPAFrameYOuterSupport\-$APA_i">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="APAFrameYOuterSupport\-$APA_i"/>
    </volume>

EOF

}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ place_APA ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Must be called only within gen_Cryostat(), 


# $_[0] = x APA center
# $_[1] = y APA center
# $_[2] = z APA center
# $_[3] = APA number
#   convention: 0 and 3 are Largest, 1 is Mid, 2 is Smallest

sub place_APA()
{

    $APA_i = $_[3];

####################################################################
################# APA Frame and Paddle Dimensions ##################

    $APA_y = $APAHeight[$APA_i];

$APAFrameZSide_x = $APAFrame_x;
$APAFrameZSide_y = 4*$inch;
$APAFrameZSide_z = $APA_z;

$APAFrameYSide_x = $APAFrame_x;
$APAFrameYSide_y = $APA_y-2*$APAFrameZSide_y;
$APAFrameYSide_z = 4*$inch;

# Two outer Y supports will sandwich the light paddles
$APAFrameYOuterSupport_x = ($APAFrame_x-$LightPaddle_x)/2;
$APAFrameYOuterSupport_y = $APA_y-2*$APAFrameZSide_y;
$APAFrameYOuterSupport_z = 4*$inch;

# if there were an inner support to fill the hole
$APAFrameYInnerSupport_x = $LightPaddle_x;

$EdgeFrameSteelThickness = 0.12*$inch;
$InnerFrameSteelThickness = 0.062*$inch;

#TODO: Place TPCs in here as well, but for now just do frames and paddles



   # First put in the frame
  print CRYO <<EOF;
      <!-- <physvol>
        <volumeref ref="volAPAFrameYOuterSupport\-$APA_i"/>
        <position name="posAPAFrameYOuterSupportNeg\-$APA_i" unit="cm" 
	x="$_[0] - ($APAFrameYOuterSupport_x + $APAFrameYInnerSupport_x/2 - $EdgeFrameSteelThickness/2)" 
	y="$_[1]" 
	z="$_[2]"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameYOuterSupport\-$APA_i"/>
        <position name="posAPAFrameYOuterSupportPos\-$APA_i" unit="cm" 
	x="$_[0] + ($APAFrameYOuterSupport_x + $APAFrameYInnerSupport_x/2 - $EdgeFrameSteelThickness/2)" 
	y="$_[1]" 
	z="$_[2]"/> 
        <rotationref ref="rIdentity"/>
      </physvol>  -->
      <physvol>
        <volumeref ref="volAPAFrameYSide\-$APA_i"/>
        <position name="posAPAFrameYSideNeg\-$_[3]" unit="cm" 
	x="$_[0]" 
	y="$_[1]" 
	z="$_[2] - $APA_z/2 + $APAFrameYSide_z/2"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameYSide\-$APA_i"/>
        <position name="posAPAFrameYSidePos\-$_[3]" unit="cm" 
	x="$_[0]" 
	y="$_[1]" 
	z="$_[2] + $APA_z/2 - $APAFrameYSide_z/2"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameZSide\-$APA_i"/>
        <position name="posAPAFrameZSidePos\-$_[3]" unit="cm" 
	x="$_[0]" 
	y="$_[1] + $APA_y/2 - $APAFrameZSide_y/2" 
	z="$_[2]"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameZSide\-$APA_i"/>
        <position name="posAPAFrameZSideNeg\-$_[3]" unit="cm" 
	x="$_[0]" 
	y="$_[1]  - $APA_y/2 + $APAFrameZSide_y/2" 
	z="$_[2]"/> 
        <rotationref ref="rIdentity"/>
      </physvol>

EOF


   # Now loop through paddle y positions and place them
  for( $p=0; $p<$nPaddlesInAPA[$APA_i]; $p++ ){

print CRYO <<EOF;

     <physvol>
       <volumeref ref="volOpDetSensitive"/>
       <position name="posPaddle\-$p\-APA\-$APA_i" unit="cm" 
         x="$_[0]" 
	 y="$_[1] - $APAHeight[$APA_i]/2 + $PaddleYPositions[$APA_i][$p]" 
         z="$_[2]"/>
       <rotationref ref="rIdentity"/>
     </physvol>

EOF

  }

}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++ gen_Enclosure +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Enclosure()
{

# Create the detector enclosure fragment file name,
# add file to list of output GDML fragments,
# and open it
    $ENCL = "35t_DetEnclosure" . $suffix . ".gdml";
    push (@gdmlFiles, $ENCL);
    $ENCL = ">" . $ENCL;
    open(ENCL) or die("Could not open file $ENCL for writing");


# The standard XML prefix and starting the gdml
    print ENCL <<EOF;
<?xml version='1.0'?>
<gdml>
EOF


# All the detector enclosure solids.
print ENCL <<EOF;
<solids>

    <box name="DetEnclosure" lunit="cm" 
      x="$DetEnc_x"
      y="$DetEnc_y"
      z="$DetEnc_z"/>

</solids>
EOF



# Detector enclosure structure
    print ENCL <<EOF;
<structure>

    <volume name="volDetEnclosure">
      <materialref ref="Concrete"/>
      <solidref ref="DetEnclosure"/>
      <physvol>
        <volumeref ref="volCryostat"/>
        <position name="posCryo" unit="cm" x="0" y="$CryoInDetEnc_ypos" z="0"/>
      </physvol>
    </volume>
EOF

print ENCL <<EOF;
</structure>
</gdml>
EOF

close(ENCL);
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++ gen_World +++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_World()
{

# Create the WORLD fragment file name,
# add file to list of output GDML fragments,
# and open it
    $WORLD = "35t_World" . $suffix . ".gdml";
    push (@gdmlFiles, $WORLD);
    $WORLD = ">" . $WORLD;
    open(WORLD) or die("Could not open file $WORLD for writing");


# The standard XML prefix and starting the gdml
    print WORLD <<EOF;
<?xml version='1.0'?>
<gdml>
EOF


# All the World solids.
print WORLD <<EOF;
<solids>
    <box name="World" lunit="cm" 
      x="$World_x" 
      y="$World_y" 
      z="$World_z"/>
</solids>
EOF

# World structure
print WORLD <<EOF;
<structure>
    <volume name="volWorld" >
      <materialref ref="Air"/>
      <solidref ref="World"/>
      <physvol>
        <volumeref ref="volDetEnclosure"/>
        <position name="posDetEnclosure" unit="cm" x="$OriginXSet" y="$OriginYSet" z="$OriginZSet"/>
      </physvol>
    </volume>
</structure>
</gdml>
EOF

# make_gdml.pl will take care of <setup/>

close(WORLD);
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++ write_fragments ++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub write_fragments()
{
   # This subroutine creates an XML file that summarizes the the subfiles output
   # by the other sub routines - it is the input file for make_gdml.pl which will
   # give the final desired GDML file. Specify its name with the output option.
   # (you can change the name when running make_gdml)

    if ( ! defined $output )
    {
	$output = "-"; # write to STDOUT 
    }

    # Set up the output file.
    $OUTPUT = ">" . $output;
    open(OUTPUT) or die("Could not open file $OUTPUT");

    print OUTPUT <<EOF;
<?xml version='1.0'?>

<!-- Input to Geometry/gdml/make_gdml.pl; define the GDML fragments
     that will be zipped together to create a detector description. 
     -->

<config>

   <constantfiles>

      <!-- These files contain GDML <constant></constant>
           blocks. They are read in separately, so they can be
           interpreted into the remaining GDML. See make_gdml.pl for
           more information. 
	   -->
	   
EOF

    foreach $filename (@defFiles)
    {
	print OUTPUT <<EOF;
      <filename> $filename </filename>
EOF
    }

    print OUTPUT <<EOF;

   </constantfiles>

   <gdmlfiles>

      <!-- The GDML file fragments to be zipped together. -->

EOF

    foreach $filename (@gdmlFiles)
    {
	print OUTPUT <<EOF;
      <filename> $filename </filename>
EOF
    }

    print OUTPUT <<EOF;

   </gdmlfiles>

</config>
EOF

    close(OUTPUT);
}
