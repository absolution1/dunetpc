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
#use strict;
#use warnings;
#use Switch;
use Math::Trig;
use Getopt::Long;
use Math::BigFloat;
Math::BigFloat->precision(-15);

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


$UVReadoutBoardPitch = 0.698516;
#$UVReadoutBoardPitch = 0.7;

$UAng[0] = 45.705;
$VAng[0] = 44.274;
$UAng[1] = 45.705;
$VAng[1] = 44.274;
$UAng[2] = 45.705;
$VAng[2] = 44.274;
$UAng[3] = 45.705;
$VAng[3] = 44.274;
# WARM
# $UAng[0] = 45.707;
# $VAng[0] = 44.275;
# $UAng[1] = 45.707;
# $VAng[1] = 44.275;
# $UAng[2] = 45.707;
# $VAng[2] = 44.275;
# $UAng[3] = 45.707;
# $VAng[3] = 44.275;

$UWirePitch = $UVReadoutBoardPitch*cos(deg2rad($UAng[0]));
$VWirePitch = $UVReadoutBoardPitch*cos(deg2rad($VAng[0]));
$XWirePitch = 0.449055;


##################################################################
######################## TPC parameters ##########################

$G10thickness = 0.3155;
$WrapCover    = 0.1578;
#$G10thickness = $inch/8;
#$WrapCover    = $inch/16;

$LongDrift               = 222.46;
$ShortDrift              = 27.16;
$APAFrame_x              = 5.0661; # ~2in -- this does not include the wire spacing

$TPCWireThickness        = 0.015;
$TPCWirePlaneThickness   = $TPCWireThickness;
$APAWirePlaneSpacing     = 0.4730488 + $TPCWirePlaneThickness; # center to center spacing between all of the wire planes (g, u, v, and x)

# At creation of the plane volumes, the y and z boundaries will be increased
# by this much at each of the 4 edges. this is so the corners of the wire 
# tubes don't extrude. For all other purposes, the plane dimensions stay as originally defined
$UVPlaneBoundNudge = $TPCWireThickness;

# The following are all widths about the same z center,
# namely the center of the corresponding APA
#$Zactive_z    = 49.8441 + $TPCWireThickness;
$Zactive_z    = 112*$XWirePitch + $TPCWireThickness;
$APAFrame_z   = 50.2619;
#$Vactive_z    = 50.8929 - 2*$G10thickness;
#$Uactive_z    = 51.5240 - 2*$G10thickness;
$Vactive_z    = $APAFrame_z;
$Uactive_z    = $APAFrame_z + 2*$G10thickness;
$APAphys_z    = 51.8395; 


 # NUMBER VERTICAL WIRES = (Zview_z / pitch) + 1
   #   Since Zview_z is defined in docdb 7550 to be distance 
   #   between outer vertical wires, + 1 since the floor of this
   #   division will be one under, giveing the amt of spaces, not wires
 # POSITIONING:  plane centered between +/- APAActive_z/2

# Let APAs be numbered as follows
#  0 & 3 - Largest
#  1 - Middle (Middle Top)
#  2 - Smallest (Middle Bottom)
#  APA heights and positions will be indexed by APA number

$APAFrame_y[0] = 203.06;
$APAFrame_y[1] = 119.29;
$APAFrame_y[2] = 91.37;
$APAFrame_y[3] = $APAFrame_y[0];

# the physical assembply includes 4 boards to anchor Grid, U, V, and Z, then a cover 
for($apa = 0; $apa < 4; ++$apa){
    $APAphys_y[$apa]    = $APAFrame_y[$apa] + 4*$G10thickness + $WrapCover; 
}

# The field cage dictates the boundaries of the active volumes.
$ReadoutBoardOverlap = 7.61; #board overlaps wires, chop this off of their active height
$OuterWireToCage      = 2.33;
$TopAPAToCage         = -5.6647; # negative because the cage is actually lower than the APA frame
$BottomTallAPAToCage  = 0.673;
$BottomShortAPAToCage =    $BottomTallAPAToCage
                        - ($APAphys_y[1] + $APAGap_y + $APAphys_y[2] - $APAphys_y[0]); # also negative

for($apa = 0; $apa < 4; ++$apa){

    $overlap = $ReadoutBoardOverlap;

      # the smallest bottom APA wire planes extrude ~1cm past the field cage
      # before being overlapped by the readout board. For the purposes of the 
      # LArSoft Geometry interface, we want to end the wires at the edge of the
      # active volume, or rather the edge of the cage.
    if($apa==2){ $overlap = -$BottomShortAPAToCage + $UVPlaneBoundNudge; } 
                            # negative sign so overlap is positive

      # each view has its own G10 board to wrap around at the bottom
      # and is covered by the readout board at the top
    $Zactive_y[$apa]    = $APAFrame_y[$apa] + 0*$G10thickness - $overlap;
    $Vactive_y[$apa]    = $APAFrame_y[$apa] + 1*$G10thickness - $overlap;
    $Uactive_y[$apa]    = $APAFrame_y[$apa] + 2*$G10thickness - $overlap;

}


$APAGap_y     =    0.0845;  #separation between APAs along the incident beam axis
$APAGap_z     =    0.0845;  #separation between APAs along the vertical axis


  # include APA spacing in y and z so volTPCs touch in y and z directions with correct APA
  # spacing - this makes for smoother event generation. 

$TPCLongDrift_x  = $LongDrift  + 3*$APAWirePlaneSpacing + $TPCWirePlaneThickness;
$TPCShortDrift_x = $ShortDrift + 3*$APAWirePlaneSpacing + $TPCWirePlaneThickness;



############################################################
#################### AuxDet parameters #####################

$AuxDetHousingHeight = 27.72*$inch;
$AuxDetHousingLongSide = 13.07*$inch;
$AuxDetHousingShortSide = 10.64*$inch;
$AuxDetHousingThickness = 3;
$AuxDetScintillatorHeight = 24.80*$inch;
$AuxDetScintillatorLongSide = 12.82*$inch;
$AuxDetScintillatorShortSide = 10.65*$inch;
$AuxDetScintillatorThickness = 0.95;

# email from Michelle:
# the housing (AL box) is 13.07/10.64 x 27.72 (inches)
# the actual scintillator is 12.82/10.65 x 24.80 (inches)


# AuxDet configuration:
#
#  /  \ _______ /  \ _______   } This spacing is only in the housing,
# /    \\     //    \\     /      so it is not implimented in AuDet
#/______\\   //______\\   /       which only represents the scintillator
#         \_/          \_/


# email from Michelle:
# Assume that the telescope scintillator is a rectangle, 64 3/4" x 6 5/8" x 3/4".

$AuxDetBSUScintWidth = (6+5/8)*$inch;
$AuxDetBSUScintHeight = 0.75*$inch;
$AuxDetBSUScintLength = 64.75*$inch;


# these positions pulled from a coordinate system with positive x and z are reversed
# positions are from a particular corner of the scintillators
############################
############################
@BSULayer1_xpos = (  16.21,  9.46,   2.71,  -4.04,
		    -17.33, -24.08, -30.83, -37.58,
		    -44.33, -51.08, -57.83, -64.58, -71.33 );
for($i = 0; $i < 13; ++$i){ 
    $BSULayer1_xpos[$i] *= -$inch; # negative reverses coordinate system back
    $BSULayer1_xpos[$i] -= $AuxDetBSUScintWidth/2; # move to center of volume
}
$BSULayer1_ypos = 75.14*$inch;
$BSULayer1_zpos = -1.21*$inch;
$BSULayer1_zpos = -1*$BSULayer1_zpos + $AuxDetBSUScintLength/2; # correct coordinates and move to center of scint
############################
############################
$BSULayer2_xpos = -67.24*$inch;
$BSULayer2_xpos = -1*$BSULayer2_xpos - $AuxDetBSUScintLength/2;
$BSULayer2_ypos = 89.49*$inch;
#@BSULayer2_zpos = ( -0.56, -7.31,  -18.84, -25.59,
#		    -32.34, -39.09, -45.84, -52.59, -59.34 );
# assume pitch of 6.75:
@BSULayer2_zpos = (  6.19,  -0.56,  -7.31, -18.84, -25.59,
		    -32.34, -39.09, -45.84, -52.59, -59.34 );
for($i = 0; $i < 10; ++$i){ 
    $BSULayer2_zpos[$i] *= -$inch; 
    $BSULayer2_zpos[$i] -= $AuxDetBSUScintWidth/2; 
}
############################
############################
@BSULayer3_xpos = (  16.21,  9.46,   2.71,  -4.04,
		    -10.79, -17.54, -24.29, -31.04,
		    -37.79, -44.54, -51.29, -58.04, 
		    -64.79, -71.54, -78.29, -85.04  );
for($i = 0; $i < 16; ++$i){ 
    $BSULayer3_xpos[$i] *= -$inch; 
    $BSULayer3_xpos[$i] -= $AuxDetBSUScintWidth/2; 
}
$BSULayer3_ypos = 306.55*$inch;
$BSULayer3_zpos = -1.21*$inch;
$BSULayer3_zpos = -1*$BSULayer3_zpos + $AuxDetBSUScintLength/2;
############################
############################
$BSULayer4_xpos = -67.24*$inch;
$BSULayer4_xpos = -1*$BSULayer4_xpos - $AuxDetBSUScintLength/2;
$BSULayer4_ypos = 310.80*$inch;
@BSULayer4_zpos = (  1.42,  -5.34,  -12.09, -18.84,
		    -25.59, -32.34, -39.09, -45.84,
		    -52.59, -59.34 );
for($i = 0; $i < 10; ++$i){ 
    $BSULayer4_zpos[$i] *= -$inch; 
    $BSULayer4_zpos[$i] -= $AuxDetBSUScintWidth/2; 
}
############################
## Layer 5 is not included in the final gdml
############################
@BSULayer5_xpos = (  16.21,  9.46,   2.71,  -4.04,
		    -10.79, -17.54, -24.29, -31.04,
		    -37.79, -44.54, -51.29, -58.04,
		    -64.79, -71.54, -78.29, -85.04 );
for($i = 0; $i < 16; ++$i){ 
    $BSULayer5_xpos[$i] *= -$inch; 
    $BSULayer5_xpos[$i] -= $AuxDetBSUScintWidth/2; 
}
$BSULayer5_ypos = 314.55*$inch;
$BSULayer5_zpos = -1.21*$inch;
$BSULayer5_zpos = -1*$BSULayer5_zpos + $AuxDetBSUScintLength/2;



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

$PaddleYPositions[0][0] = $APAFrame_y[0]/2; # this puts it in the y center 
$PaddleYPositions[1][0] = $APAFrame_y[1] - 4*$inch - $LightPaddle_y/2; 
$PaddleYPositions[2][0] = $APAFrame_y[2] - 4*$inch - $LightPaddle_y/2;
$PaddleYPositions[3][0] = $PaddleYPositions[0][0];

$PaddleYPositions[0][1] = $APAFrame_y[0] - 4*$inch - $LightPaddle_y/2;
$PaddleYPositions[3][1] = $PaddleYPositions[0][1];

$PaddleYPositions[0][2] = 4*$inch + $LightPaddle_y/2;
$PaddleYPositions[3][2] = $PaddleYPositions[0][2];



##################################################################
###################### Cryostat parameters #######################


$HeightGaseousAr       =    14;
$FloorToSmallAPAFrame  =    22.432;
$APAToTopCryo          =    34.441; 
$ClosestAPAToEastWall  =    33.225; # To G10 cover, places all of the APAs in z

#  Long       West      Short
#  Drift   __________   Drift
#         |          |
#         |       |  |
#  South  |       |  |  North Wall
#         |__________|
#
#             East

$CPAToFloor     = 26.955;
$CPAToCeiling   = 36.010;
$CPAToEastWall  = 27.994;
$CPAToWestWall  = 76.234;
$CPAToSouthWall = 121.945;
$CPAToNorthWall = 22.907;




$CPATube_OD = 5.066;
$CPATube_ID = 4.747;

$CPATubeYSide_CenterToCenter = 166.269; # length in *z* direction
$CPATubeZSide_CenterToCenter = 207.532; # height in *y* direction

$Cathode_x                 =    0.016;
$Cathode_y                 =    $CPATubeZSide_CenterToCenter - $CPATube_OD;
$Cathode_z                 =    $CPATubeYSide_CenterToCenter - $CPATube_OD;    

# Liquid and Gaseous Argon dimensions 
$Argon_x	       =    $TPCShortDrift_x 
                          + $APAFrame_x 
                          + $TPCLongDrift_x 
                          + $CPAToNorthWall 
                          + $CPAToSouthWall + $CPATube_OD;
# note that the height y includes liquid and gaseous argon
$Argon_y	       =    $APAphys_y[1] 
                          + $APAGap_y
                          + $APAphys_y[2]
                          + $APAToTopCryo
                          + $FloorToSmallAPAFrame;  # assuming mid_y+smallest_y > largest_y

$Argon_z	       =    $CPATubeYSide_CenterToCenter
                          + $CPAToEastWall 
                          + $CPAToWestWall;


# try hardcoding parameters now that placements/dimensions are more accurate
# inside cryostat parameters
#$Argon_x  =  400.497;
#$Argon_y  =  270.497;
#$Argon_z  =  270.497;
$NeckInside_x   =  100;
$NeckInside_y   =  75;

# Cryostat Dimensions
$SteelShellThickness	       =    0.5*$inch;

$Cryostat_x	       =    $Argon_x; # move the steel shell out of volCryostat
$Cryostat_y	       =    $Argon_y;
$Cryostat_z	       =    $Argon_z;



##################################################################
################# Detector Enclosure parameters ##################


# Around the cryostat
$ConcretePadding       =    30;
$FoamPadding           =    39.75148;
$TotalPadding	       =    $ConcretePadding + $FoamPadding;
$CryoWithPadding_x         =   $Argon_x + 2*$SteelShellThickness + 2*$TotalPadding;
$CryoWithPadding_y_noneck  =   $Argon_y + 2*$SteelShellThickness + 2*$FoamPadding + $ConcretePadding;
$CryoWithPadding_y_neck    =   ($Argon_y + $NeckInside_y) + 2*$SteelShellThickness + $TotalPadding;
$CryoWithPadding_z         =    $Argon_z + 2*$SteelShellThickness + 2*$TotalPadding;

print "CryoWithPadding: x=$CryoWithPadding_x, y=$CryoWithPadding_y_noneck ($CryoWithPadding_y_neck), z=$CryoWithPadding_z\n";
print "Cryo: x=$Argon_x,  y=$Argon_y (+$NeckInside_y),  z=$Argon_z\n";


# The actual enclosure
$TrenchLength = 9*$CryoWithPadding_x; # make this a reasonable length, still guessing
     # quick measurments by Michelle Stancari, definitely good enough for now
$TrenchWallThickness = 14*$inch;
$WalkwayWidth = 34.5*$inch;
$EastWallToTrench = 10*$inch;
$PlateAToGroundLevel = 50*$inch;


# Bird's eye view of trench, walkway is on the west side
#    |          |
#    |          |<-- trench wall
#    |          |
#    |   ___N__ |
#    |  | ____ ||
#    |  |      ||
#    | W|______||E    +z <---|  (+y out of screen)
#    |  |______||            |
#    |      ^   |            v
#    |     Neck |           +x
#    |          |
#    |          |
#    |          |

# the container volume dimensions
$DetEnc_x	       =    $TrenchLength;
$DetEnc_y	       =    $TrenchWallThickness  # use this for thickness of floor too
                           + $CryoWithPadding_y_noneck 
                           + $PlateAToGroundLevel; # Make ground level the top of the DetEnc
$DetEnc_z              =    $CryoWithPadding_z 
                           + $WalkwayWidth + $EastWallToTrench
                           + 2*$TrenchWallThickness;


$posCryoInDetEnc_x     =  0;
$posCryoInDetEnc_y     =  - $DetEnc_y/2 + $TrenchWallThickness + $TotalPadding + $SteelShellThickness + $Argon_y/2;
$posCryoInDetEnc_z     =  - $DetEnc_z/2 + $TrenchWallThickness + $EastWallToTrench + $CryoWithPadding_z/2;




$PosDirCubeSide = 0;
if (defined $helpcube)
{
$PosDirCubeSide = $ArToAr; #seems to be a good proportion
}


# The world dimensions are critical in the CRY cosmics generator
#   following uboone's lead, make world much larger
#   the cry helper needs a lot of room


$World_x = 5*$DetEnc_x;
$World_y = 25*$DetEnc_y;
$World_z = 25*$DetEnc_z;

$BermRadius = 5*12*$inch;





##################################################################
######################### TPC positions ##########################


$APA_Xcenter     =   $Argon_x/2 
                   - $CPAToSouthWall # to center CPATube
                   - $CPATube_OD
                   - $TPCLongDrift_x
                   - $APAFrame_x/2;

# 0: One of the 2 identical tall APAs (Largest), call it the "upstream" one
$APACenter[0][0] =   $APA_Xcenter;
$APACenter[0][1] =   $Argon_y/2
                   - $APAToTopCryo         # This subsumes the half vertical gap on the top...
                   - $APAphys_y[0]/2;      # ... so use APAphys_y instead of TPC_y
$APACenter[0][2] = - $Argon_z/2
                   + $ClosestAPAToEastWall  # ..Similarly, this already steps into the APAGap_z/2
                   + $APAphys_z/2;


# 1: The top middle APA (Mid)
$APACenter[1][0] =   $APA_Xcenter;
$APACenter[1][1] =   $Argon_y/2
                   - $APAToTopCryo
                   - $APAphys_y[1]/2;
$APACenter[1][2] =   $APACenter[0][2]
                   + $APAphys_z + $APAGap_z;

# 2: The bottom middle APA (Smallest)
$APACenter[2][0] =   $APA_Xcenter;
$APACenter[2][1] =   $APACenter[1][1] # place relative to APA above it
                   - $APAphys_y[1]/2
                   - $APAGap_y
                   - $APAphys_y[2]/2;
$APACenter[2][2] =   $APACenter[0][2]
                   + $APAphys_z + $APAGap_z;

# 3: The other tall APA, call it the "downstream" one
$APACenter[3][0] =   $APACenter[0][0];
$APACenter[3][1] =   $APACenter[0][1];
$APACenter[3][2] =   $APACenter[1][2]
                   + $APAphys_z + $APAGap_z;

$posTPCShortDrift_x  =    $APACenter[0][0]
                        - $APAFrame_x/2
                        - $TPCShortDrift_x/2;

$posTPCLongDrift_x   =    $APACenter[0][0]
                        + $APAFrame_x/2
                        + $TPCLongDrift_x/2;




# We want the active volumes to extend all the way to the field cage
# so that LArG4 tracks everything inside the cage.
$TPCActive_z[0] =   $APAphys_z
                  - ($APAphys_z-$Zactive_z)/2
                  + $APAGap_z/2
                  + $OuterWireToCage;
$TPCActive_z[1] =   $APAphys_z + $APAGap_z;
$TPCActive_z[2] =   $APAphys_z + $APAGap_z;
$TPCActive_z[3] =   $TPCActive_z[0];


$TPCActive_y[0] =   $APAphys_y[0]
                  + $TopAPAToCage
                  + $BottomTallAPAToCage;
$TPCActive_y[1] =   $APAphys_y[1]
                  + $APAGap_y/2
                  + $TopAPAToCage;
$TPCActive_y[2] =   $APAphys_y[2]
                  + $APAGap_y/2
                  + $BottomShortAPAToCage
                  + $UVPlaneBoundNudge; # Without this relatively negligible addition, the nudge
                                        # further downstream on the plane causes the planes to extrude 
                                        # the TPC volume, whose height is TPCActive_y.
                                        # The plane nudge is to avoid many overlaps between the wire 
                                        # corners and the plane boundaries.
$TPCActive_y[3] =   $TPCActive_y[0];




$TPCCenter[0][2] =   $APACenter[0][2] 
                   - $Zactive_z/2
                   - $OuterWireToCage
                   + $TPCActive_z[0]/2;
$TPCCenter[1][2] =   $APACenter[1][2]; 
$TPCCenter[2][2] =   $APACenter[2][2]; 
$TPCCenter[3][2] =   $APACenter[3][2] 
                   + $Zactive_z/2
                   + $OuterWireToCage
                   - $TPCActive_z[0]/2;



$TPCCenter[0][1] =   $APACenter[0][1] 
                   - $APAphys_y[0]/2
                   - $BottomTallAPAToCage
                   + $TPCActive_y[0]/2;
$TPCCenter[1][1] =   $APACenter[1][1]
                   + $APAphys_y[1]/2
                   + $TopAPAToCage
                   - $TPCActive_y[1]/2;
$TPCCenter[2][1] =   $APACenter[2][1]
                   - $APAphys_y[2]/2
                   - $BottomShortAPAToCage
                   + $TPCActive_y[2]/2;
$TPCCenter[3][1] =   $TPCCenter[0][1];

    



  # We want the world origin to be at the very front of the fiducial volume.
  # move it to the front of the enclosure, then back it up through the concrete/foam, 
  # then through the Cryostat shell, then through the upstream dead LAr (including the
  # dead LAr on the edge of the TPC, but this is covered in $UpstreamLArPadding).
  # This is to be added to the z position of every volume in volWorld

$OriginZSet =       
    - $DetEnc_z/2
    + $TrenchWallThickness
    + $WalkwayWidth
    + $TotalPadding 
    + $SteelShellThickness
    + $Argon_z/2    # at this point, we are at the center of the cryostat...
    - $APACenter[0][2] # ... and now at the center of the East-most APA
    + $Uactive_z/2;

  # We want the world origin to be vertically centered between the two stacked APAs.
  # (for now, that is, so the sorting works. this is quite asymetric, but then again
  #  so is the entire 35t geometry. this may be kept.)
  # the cryostat sits on top of concrete padding, move the detector enclosure back
  # and then move the world origin to the bottom of the smallest/lowest TPC, then
  # and then up through the TPC, then back up to being centered between the stacked APAs.
  # This is to be added to the y/x position of every volume in volWorld

$OriginYSet =       
      $DetEnc_y/2
    - $TrenchWallThickness
    - $TotalPadding
    - $SteelShellThickness
    - $Argon_y
    + $APAToTopCryo
    + $APAphys_y[1]
    + $APAGap_y/2;




#$OriginXSet             =       $DetEnc_x/2 
#                              - $TotalPadding
#                              - $SteelShellThickness
#                              - $CPAToNorthWall
#                              - $CPATube_OD
#                              - $ShortDrift ... through APA frame

$OriginXSet  =    
    - $CryoWithPadding_x/2 
    + $TotalPadding
    + $SteelShellThickness
    + $Argon_x/2   # at this point, we are at the center of the cryostat...
    - $APA_Xcenter    # ... and now at the APA's center x coordinate
    - $APAFrame_x/2
    - 3*$APAWirePlaneSpacing
    - $TPCWirePlaneThickness;



##
## Correct Telescope Aux Det vertical position if necessary
##
if(   $posCryoInDetEnc_y 
    + $Argon_y/2 
    + $SteelShellThickness 
    + $FoamPadding 
    + $AuxDetBSUScintHeight/2   >   $BSULayer1_ypos - $OriginYSet ){

    $nudge_y = (  $posCryoInDetEnc_y 
		+ $Argon_y/2  
		+ $SteelShellThickness 
		+ $FoamPadding 
		+ $AuxDetBSUScintHeight/2 
		+ $OriginYSet              ) - $BSULayer1_ypos;

    print "--- NOTE: nudging telescope counters $nudge_y cm up to sit on top of GDML cryostat\n";

    $BSULayer1_ypos += $nudge_y;
    $BSULayer2_ypos += $nudge_y;
    $BSULayer3_ypos += $nudge_y;
    $BSULayer4_ypos += $nudge_y;
    $BSULayer5_ypos += $nudge_y;

}



##################################################################
######################### CPA positions ##########################

#  Long       West      Short
#  Drift   __________   Drift
#         |          |
#         |       |  |
#  South  |       |  |  North Wall
#         |__________|
#
#             East

#$posCPAShortDrift_x  =  - $Argon_x/2 + $CPAToNorthWall + $CPATube_OD/2;
#$posCPALongDrift_x   =    $Argon_x/2 - $CPAToSouthWall - $CPATube_OD/2;

# ^^^ these would be ideal, but quick fix to check in
# v3 without overlaps:

#$posCPAShortDrift_x  =  - $Argon_x/2 + $CPAToNorthWall + $CPATube_OD/2;
#$posCPALongDrift_x   =    $Argon_x/2 - $CPAToSouthWall - $CPATube_OD/2;
$posCPAShortDrift_x  =    $posTPCShortDrift_x - $TPCShortDrift_x/2 
                                - $CPATube_OD; #<-- temp overlap fix
$posCPALongDrift_x   =    $posTPCLongDrift_x + $TPCLongDrift_x/2
                                + $CPATube_OD; #<-- temp overlap fix


$posCPAShortDrift_y  =  - $Argon_y/2 + $CPAToFloor 
                        + ($CPATubeZSide_CenterToCenter + $CPATube_OD)/2;
$posCPAShortDrift_z  =  - $Argon_z/2 + $CPAToEastWall
                        + ($CPATubeYSide_CenterToCenter + $CPATube_OD)/2;

$posCPALongDrift_y  =    $posCPAShortDrift_y;
$posCPALongDrift_z  =    $posCPAShortDrift_z;





##################################################################
#################### Bar Fiber Module numbers ####################
$numberofbarmodules=4;
$numberoffibermodules=3;
$numberofplankmodules=1;

$PaddleCenterX=$APA_Xcenter;
$PaddleCenterY[0][0]=$APACenter[0][1]-$APAFrame_y[0]/2+ $PaddleYPositions[0][0];
$PaddleCenterY[0][1]=$APACenter[0][1]-$APAFrame_y[0]/2+ $PaddleYPositions[0][1];
$PaddleCenterY[0][2]=$APACenter[0][1]-$APAFrame_y[0]/2+ $PaddleYPositions[0][2];
$PaddleCenterY[1][0]=$APACenter[1][1]-$APAFrame_y[1]/2+ $PaddleYPositions[1][0];
$PaddleCenterY[2][0]=$APACenter[2][1]-$APAFrame_y[2]/2+ $PaddleYPositions[2][0];
$PaddleCenterY[3][0]=$PaddleCenterY[0][0];
$PaddleCenterY[3][1]=$PaddleCenterY[0][1];
$PaddleCenterY[3][2]=$PaddleCenterY[0][2];
$PaddleCenterZ[0]=$APACenter[0][2];
$PaddleCenterZ[1]=$APACenter[1][2];
$PaddleCenterZ[2]=$APACenter[2][2];
$PaddleCenterZ[3]=$APACenter[3][2];
################

#+++++++++++++++++++++++++ End defining variables ++++++++++++++++++++++++++



# run the sub routines that generate the fragments

gen_Define(); 	 # generates definitions at beginning of GDML
gen_Materials(); # generates materials to be used


# gen_TPC( x dimension, y, z, string appended to TPC for name, APA number)
# generate the short drift and long drift sides of the APA as seperate TPCs

open(my $wout, '>', 'gdmlWireCenters.txt');

    gen_TPC( $TPCLongDrift_x,  $TPCActive_y[0], $TPCActive_z[0], 'LargestLongDriftUpstream', 0);
    gen_TPC( $TPCShortDrift_x, $TPCActive_y[0], $TPCActive_z[0], 'LargestShortDriftUpstream', 0);

    gen_TPC( $TPCLongDrift_x,  $TPCActive_y[2], $TPCActive_z[2], 'SmallestLongDrift', 2);
    gen_TPC( $TPCShortDrift_x, $TPCActive_y[2], $TPCActive_z[2], 'SmallestShortDrift', 2);

    gen_TPC( $TPCLongDrift_x,  $TPCActive_y[1], $TPCActive_z[1], 'MidLongDrift', 1);
    gen_TPC( $TPCShortDrift_x, $TPCActive_y[1], $TPCActive_z[1], 'MidShortDrift', 1);

    gen_TPC( $TPCLongDrift_x,  $TPCActive_y[0], $TPCActive_z[3], 'LargestLongDriftDownstream', 3);
    gen_TPC( $TPCShortDrift_x, $TPCActive_y[0], $TPCActive_z[3], 'LargestShortDriftDownstream', 3);

close $wout;

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

   <position name="posTPCLargestShortDrift_Pos"  unit="cm" x="$posTPCShortDrift_x" y="$TPCCenter[3][1]"  z="$TPCCenter[3][2]"/>
   <position name="posTPCLargestLongDrift_Pos"   unit="cm" x="$posTPCLongDrift_x"  y="$TPCCenter[3][1]"  z="$TPCCenter[3][2]"/>
   <position name="posTPCLargestShortDrift_Neg"  unit="cm" x="$posTPCShortDrift_x" y="$TPCCenter[0][1]"  z="$TPCCenter[0][2]"/>
   <position name="posTPCLargestLongDrift_Neg"   unit="cm" x="$posTPCLongDrift_x"  y="$TPCCenter[0][1]"  z="$TPCCenter[0][2]"/>
   <position name="posTPCSmallestShortDrift"     unit="cm" x="$posTPCShortDrift_x" y="$TPCCenter[2][1]"  z="$TPCCenter[2][2]"/>
   <position name="posTPCSmallestLongDrift"      unit="cm" x="$posTPCLongDrift_x"  y="$TPCCenter[2][1]"  z="$TPCCenter[2][2]"/>
   <position name="posTPCMidShortDrift"          unit="cm" x="$posTPCShortDrift_x" y="$TPCCenter[1][1]"  z="$TPCCenter[1][2]"/>
   <position name="posTPCMidLongDrift"           unit="cm" x="$posTPCLongDrift_x"  y="$TPCCenter[1][1]"  z="$TPCCenter[1][2]"/>


   <position name="posCathodeLongDrift"  unit="cm" x="$posCPAShortDrift_x" y="$posCPAShortDrift_y" z="$posCPAShortDrift_z"/>
   <position name="posCathodeShortDrift" unit="cm" x="$posCPALongDrift_x" y="$posCPALongDrift_y" z="$posCPALongDrift_z"/>

   <rotation name="rAuxDetNSWallUp"    unit="deg" x="90"  y="0"  z="90"/>
   <rotation name="rAuxDetNSWallDown"  unit="deg" x="270" y="0"  z="90"/>
   <rotation name="rAuxDetEWWallUp"    unit="deg" x="90"  y="0"  z="0"/>
   <rotation name="rAuxDetEWWallDown"  unit="deg" x="270" y="0"  z="0"/>

   <position name="posCenter"           unit="cm" x="0" y="0" z="0"/>
   <rotation name="rPlus90AboutX"       unit="deg" x="90" y="0"  z="0"/>
   <rotation name="rPlus90AboutY"       unit="deg" x="0"  y="90" z="0"/>
   <rotation name="rMinus90AboutY"      unit="deg" x="0" y="270" z="0"/>
   <rotation name="rMinus90AboutYMinus90AboutX"       unit="deg" x="270" y="270" z="0"/>
   <rotation name="rPlus180AboutY"	unit="deg" x="0" y="180"   z="0"/>
   <rotation name="rPlus180AboutX"	unit="deg" x="180" y="0"   z="0"/>
   <rotation name="rPlus180AboutXandY"	unit="deg" x="180" y="180"   z="0"/>
   <rotation name="rIdentity"		unit="deg" x="0" y="0"   z="0"/>
EOF

##################################################################
###################### Bar Module Position #######################
for ($k=1; $k<$numberofbarmodules+1; ++$k)
{
if($k==1) {$APA_i=1;$p=0;$PD_zoffset=0;}
elsif($k==2) {$APA_i=0;$p=1;$PD_zoffset=0;}
elsif($k==3) {$APA_i=0;$p=2;$PD_zoffset=-0.00051;}
elsif($k==4) {$APA_i=3;$p=0;$PD_zoffset=0;}

#for ($j=1; $j<$numberofbars+1; ++$j)
#{
#$bar_z=-4.11 + 2.74*($j-1)+$PaddleCenterZ[$APA_i];
$bar_z=$PaddleCenterZ[$APA_i]+2.794*1.5+$PD_zoffset;
$bar_name="Bar" . $k . "Pos";
print DEF <<EOF;
	<position name="$bar_name" x="$PaddleCenterX" y="$PaddleCenterY[$APA_i][$p]" z="$bar_z" unit="cm"/>
EOF
#}
}

##################################################################
##################### Plank Module Position ######################
for ($k=1; $k<$numberofplankmodules+1; ++$k)
{
if($k==1) {$APA_i=3;$p=2;}

$plank_name="PlankPos";
print DEF <<EOF;
	<position name="$plank_name" x="$PaddleCenterX" y="$PaddleCenterY[$APA_i][$p]" z="$PaddleCenterZ[$APA_i]" unit="cm"/>
EOF

}
##################################################################
##################### Fiber Module Position ######################

for ($k=1; $k<$numberoffibermodules+1; ++$k)
{
if($k==1) {$APA_i=2;$p=0;$PD_zoffset=0;}
elsif($k==2) {$APA_i=0;$p=0;$PD_zoffset=0;}
elsif($k==3) {$APA_i=3;$p=1;$PD_zoffset=0.00051;}

$fiber_x=$PaddleCenterX+0.1445;
$fiber_y=$PaddleCenterY[$APA_i][$p];
$fiber_z=$PaddleCenterZ[$APA_i]+4.1905+$PD_zoffset;

$fiber_name="Fiber" . $k. "Pos";
print DEF <<EOF;
	<position name="$fiber_name" x="$fiber_x" y="$fiber_y" z="$fiber_z" unit="cm"/>
EOF

}

print DEF <<EOF;
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
  <element name="hydrogen" formula="H" Z="1">  <atom value="1.0079"/> </element>
  <element name="beryllium" formula="Be" Z="4">  <atom value="9.0121831"/>  </element>
  <element name="carbon" formula="C" Z="6">  <atom value="12.0107"/>  </element>
  <element name="nitrogen" formula="N" Z="7">  <atom value="14.0067"/> </element>
  <element name="oxygen" formula="O" Z="8">  <atom value="15.999"/> </element>
  <element name="sodium" formula="Na" Z="11"> <atom value="22.99"/>    </element>
  <element name="magnesium" formula="Mg" Z="12"> <atom value="24.305"/>   </element>
  <element name="aluminum" formula="Al" Z="13"> <atom value="26.9815"/>  </element>
  <element name="silicon" formula="Si" Z="14"> <atom value="28.0855"/>  </element>
  <element name="phosphorus" formula="P" Z="15"> <atom value="30.973"/>  </element>
  <element name="sulphur" formula="S" Z="16"> <atom value="32.065"/>  </element>
  <element name="argon" formula="Ar" Z="18"> <atom value="39.9480"/>  </element>
  <element name="potassium" formula="K" Z="19"> <atom value="39.0983"/>  </element>
  <element name="calcium" formula="Ca" Z="20"> <atom value="40.078"/>   </element>
  <element name="titanium" formula="Ti" Z="22"> <atom value="47.867"/>   </element>
  <element name="chromium" formula="Cr" Z="24"> <atom value="51.9961"/>  </element>
  <element name="iron" formula="Fe" Z="26"> <atom value="55.8450"/>  </element>
  <element name="nickel" formula="Ni" Z="28"> <atom value="58.6934"/>  </element>
  <element name="copper" formula="Cu" Z="29">  <atom value="63.546"/>  </element>
  <element name="bromine" formula="Br" Z="35"> <atom value="79.904"/> </element>

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

  <material name="Polystyrene" formula="C8H8">
   <D value="1.05" unit="g/cm3"/>
   <composite n="8" ref="carbon"/>
   <composite n="8" ref="hydrogen"/>
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
   <composite n="1" ref="carbon"/>
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
   <D value="2.58" unit="g/cm3"/>
   <fraction n="0.600" ref="SiO2"/>
   <fraction n="0.118" ref="Al2O3"/>
   <fraction n="0.001" ref="Fe2O3"/>
   <fraction n="0.224" ref="CaO"/>
   <fraction n="0.034" ref="MgO"/>
   <fraction n="0.010" ref="Na2O"/>
   <fraction n="0.013" ref="TiO2"/>
  </material>

<!-- The following fractional components are placeholders,
     to be fixed (though they have very little effect,
     as long as the density is correct) -->
  <material name="polyurethane_foam">
   <D value=".13" unit="g/cm3"/>
   <fraction n="0.600" ref="SiO2"/>
   <fraction n="0.118" ref="Al2O3"/>
   <fraction n="0.001" ref="Fe2O3"/>
   <fraction n="0.224" ref="CaO"/>
   <fraction n="0.034" ref="MgO"/>
   <fraction n="0.010" ref="Na2O"/>
   <fraction n="0.013" ref="TiO2"/>
  </material>

<!-- for the cryostat foam insulation -->
  <material name="R-PUF">
   <D value=".2525" unit="g/cm3"/>
   <fraction n="0.95" ref="polyurethane_foam"/>
   <fraction n="0.05" ref="fibrous_glass"/>
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

  <material name="Copper_Beryllium_alloy25" formula="Copper_Beryllium_alloy25">
   <D value="8.26" unit="g/cm3"/>
   <fraction n="0.981" ref="copper"/>
   <fraction n="0.019" ref="beryllium"/>
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

  <material name="Plastic" formula="Plastic">
   <D value="1.032" unit="g/cm3"/>
   <fraction n=".474" ref="carbon"/>
   <fraction n=".526" ref="hydrogen"/>
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
# $_[4] = APA number

    my $apa = $_[4];

    my $TPCActive_x   =  $_[0]-(3*$APAWirePlaneSpacing);
    my $TPCActive_y   =  $_[1] - $APAGap_y/2 - $ReadoutBoardOverlap;

    my $UAngle = $UAng[$apa];
    my $VAngle = $VAng[$apa];

    my $SinUAngle = sin( deg2rad($UAngle) );
    my $CosUAngle = cos( deg2rad($UAngle) );
    my $TanUAngle = tan( deg2rad($UAngle) );
    
    my $SinVAngle = sin( deg2rad($VAngle) );
    my $CosVAngle = cos( deg2rad($VAngle) );
    my $TanVAngle = tan( deg2rad($VAngle) );
    
    my $UWire_yint = $UWirePitch/$SinUAngle;
    my $UWire_zint = $UVReadoutBoardPitch  ;
    
    my $VWire_yint = $VWirePitch/$SinVAngle;
    my $VWire_zint = $UVReadoutBoardPitch  ;

#constructs everything inside volTPC, namely
# (moving from left to right, or from +x to -x)
#  -volCPActive
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


print $wout "\n\n\n----- Wires for $_[3] -----\n\n\n";


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
    <box name="${_[3]}UPlane" lunit="cm" 
      x="$TPCWirePlaneThickness" 
      y="$Uactive_y[$apa] + $UVPlaneBoundNudge" 
      z="$Uactive_z + $UVPlaneBoundNudge"/>
    <box name="${_[3]}VPlane" lunit="cm" 
      x="$TPCWirePlaneThickness" 
      y="$Vactive_y[$apa] + $UVPlaneBoundNudge" 
      z="$Vactive_z + $UVPlaneBoundNudge"/>
    <box name="${_[3]}ZPlane" lunit="cm" 
      x="$TPCWirePlaneThickness" 
      y="$Zactive_y[$apa]" 
      z="$Zactive_z"/>
    <box name="${_[3]}Active" lunit="cm"
      x="$TPCActive_x"
      y="$TPCActive_y[$apa]"
      z="$TPCActive_z[$apa]"/>
EOF


#++++++++++++++++++++++++++++ Wire Solids ++++++++++++++++++++++++++++++

print TPC <<EOF;

    <tube name="${_[3]}WireVert"
      rmax="0.5*$TPCWireThickness"
      z="$Zactive_y[$apa]"               
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
#$NumberCornerUWires = int( $APAFrame_z/($UWirePitch/$CosUAngle) );
$NumberCornerUWires = 72;

#$NumberCornerVWires = int( $APAFrame_z/($VWirePitch/$CosVAngle) );
$NumberCornerVWires = 72;

   # Total number of wires touching one vertical (longer) side
   # Note that the total number of wires per plane is this + another set of corner wires
$NumberSideUWires = int( $Uactive_y[$apa]/($UWirePitch/$SinUAngle) );

$NumberSideVWires = int( $Vactive_y[$apa]/($VWirePitch/$SinVAngle) );

   # Number of wires per side that aren't cut off by the corner
$NumberCommonUWires = $NumberSideUWires - $NumberCornerUWires;

$NumberCommonVWires = $NumberSideVWires - $NumberCornerVWires;

   # Number of wires on the vertical plane
   #   Since APA Active z is defined in docdb 7550 to be distance 
   #   between outer vertical wires, + 1 since the floor of this
   #   division will be one under, giveing the amt of spaces, not wires
# $NumberVerticalWires = int( $Zactive_z/$XWirePitch ) + 1;
$NumberVerticalWires = 112; # Zactive now defined in terms of 112, 


$nUchans = 2*$NumberCornerUWires;
$nVchans = 2*$NumberCornerVWires;

print $wout "$nUchans U channels\n";
print $wout "$nVchans V channels\n";
print $wout "$NumberVerticalWires Z channels per side\n";

}


my $FirstUWireOffset = .35 + $G10thickness + 2*$G10thickness*$TanUAngle - $UVReadoutBoardPitch;
my $FirstVWireOffset = .35; # doesnt include a G10 board in width

my $FirstTopUWire_yspan =
    $Uactive_y[$apa]/2
    - ( - $Uactive_y[$apa]/2
        + $FirstUWireOffset/$TanUAngle      # walk us up to the first wire
        + $UWire_yint*($NumberSideUWires-1) # up to the top of the top common wire
        - $Uactive_z/$TanUAngle             # back to the bottom of the top common wire
	+ $UWire_yint);                     # nudge up to bottom of the first top corner wire

my $FirstTopVWire_yspan =
    $Vactive_y[$apa]/2
    - ( - $Vactive_y[$apa]/2
        + $FirstVWireOffset/$TanVAngle      # walk us up to the first wire
        + $VWire_yint*($NumberSideVWires-1) # up to the top of the top common wire
        - $Vactive_z/$TanVAngle             # back to the bottom of the top common wire
	+ $VWire_yint);                     # nudge up to bottom of the first top corner wire

# The corner wires for the U plane
if ($wires_on==1)
{
    for ($i = 0; $i < $NumberCornerUWires; $i++)
    {
	$CornerUWireLength[$i] = ($FirstUWireOffset + $i*$UVReadoutBoardPitch)/$SinUAngle;

   print TPC <<EOF;
    <tube name="${_[3]}WireU$i"
      rmax="0.5*$TPCWireThickness"
      z="$CornerUWireLength[$i]"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

    }

    $CommonUWireLength = $Uactive_z/$SinUAngle;

   print TPC <<EOF;
    <tube name="${_[3]}WireUCommon"
      rmax="0.5*$TPCWireThickness"
      z="$CommonUWireLength"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

    for ($i = 0; $i < $NumberCornerUWires; $i++)
    {

	$TopCornerUWireLength[$i] = ($FirstTopUWire_yspan - $i*$UWire_yint)/$CosUAngle;

	$j = $i + $NumberSideUWires;

   print TPC <<EOF;
    <tube name="${_[3]}WireU$j"
      rmax="0.5*$TPCWireThickness"
      z="$TopCornerUWireLength[$i]"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

    }

}


# The corner wires for the V plane
if ($wires_on==1) 
{
    for ($i = 0; $i < $NumberCornerVWires; ++$i)
    {
	$CornerVWireLength[$i] = ($FirstVWireOffset + $i*$UVReadoutBoardPitch)/$SinVAngle;

	print TPC <<EOF;

    <tube name="${_[3]}WireV$i"
      rmax="0.5*$TPCWireThickness"
      z="$CornerVWireLength[$i]"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>

EOF

    }

    # The wire used many times in the middle of the V plane
    # Same subtraction as U common 

    $CommonVWireLength = $Vactive_z/$SinVAngle;

   print TPC <<EOF;
    <tube name="${_[3]}WireVCommon"
      rmax="0.5*$TPCWireThickness"
      z="$CommonVWireLength"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

    for ($i = 0; $i < $NumberCornerVWires; $i++)
    {

	$TopCornerVWireLength[$i] = ($FirstTopVWire_yspan - $i*$VWire_yint)/$CosVAngle;

	$j = $i + $NumberSideVWires;

   print TPC <<EOF;
    <tube name="${_[3]}WireV$j"
      rmax="0.5*$TPCWireThickness"
      z="$TopCornerVWireLength[$i]"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

    }


}

# make the solids only once per APA
# (here only from the lang drift TPC)
if($_[0]>100){ solid_TPCG10( $_[4],  $_[0],  $_[1],  $_[2]); }

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
      <materialref ref="Copper_Beryllium_alloy25"/>
      <solidref ref="${_[3]}WireVert"/>
    </volume>
EOF

  # Corner U wires logical volumes
  for ($i = 0; $i < $NumberCornerUWires; ++$i)
  {
  print TPC <<EOF;
    <volume name="volTPCWireU$i${_[3]}">
      <materialref ref="Copper_Beryllium_alloy25"/>
      <solidref ref="${_[3]}WireU$i"/>
    </volume>
EOF
  }


  # Top Corner U wires logical volumes
  for ($j = $NumberSideUWires; $j < $NumberSideUWires + $NumberCornerUWires; ++$j)
  {
  print TPC <<EOF;
    <volume name="volTPCWireU$j${_[3]}">
      <materialref ref="Copper_Beryllium_alloy25"/>
      <solidref ref="${_[3]}WireU$j"/>
    </volume>
EOF
  }


  # Common U wire logical volume, referenced many times
  print TPC <<EOF;
    <volume name="volTPCWireUCommon${_[3]}">
      <materialref ref="Copper_Beryllium_alloy25"/>
      <solidref ref="${_[3]}WireUCommon"/>
    </volume>
EOF

  # Corner V wires logical volumes
  for ($i = 0; $i < $NumberCornerVWires; ++$i)
  {
  print TPC <<EOF;
    <volume name="volTPCWireV$i${_[3]}">
      <materialref ref="Copper_Beryllium_alloy25"/>
      <solidref ref="${_[3]}WireV$i"/>
    </volume>
EOF

  }

  # Top Corner V wires logical volumes
  for ($j = $NumberSideVWires; $j < $NumberSideVWires + $NumberCornerVWires; ++$j)
  {
  print TPC <<EOF;
    <volume name="volTPCWireV$j${_[3]}">
      <materialref ref="Copper_Beryllium_alloy25"/>
      <solidref ref="${_[3]}WireV$j"/>
    </volume>
EOF
  }

  # Common V wire logical volume, referenced many times
  print TPC <<EOF;
    <volume name="volTPCWireVCommon${_[3]}">
      <materialref ref="Copper_Beryllium_alloy25"/>
      <solidref ref="${_[3]}WireVCommon"/>
    </volume>
EOF

}

# generate the G10 board solids and logical volumes
# make the volumes only once per APA
# (here only from the lang drift TPC)
if($_[0]>100){ vol_TPCG10( $_[4] ); }

my $lastYpos = 0;
my $lastZpos = 0;


#+++++++++++++++++++++++++ Position physical wires ++++++++++++++++++++++++++

#            ++++++++++++++++++++++  U Plane  +++++++++++++++++++++++

# Create U plane logical volume
print TPC <<EOF;
    <volume name="volTPCPlaneU${_[3]}">
      <materialref ref="LAr"/>
      <solidref ref="${_[3]}UPlane"/>
EOF


print $wout "\n-     Wires for U plane  -\n\n";
print $wout " Uplane_y: $Uactive_y[$apa]\n";
print $wout " Uplane_z: $Uactive_z\n";


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

    $FirstU_ypos = - $Uactive_y[$apa]/2 + $FirstUWireOffset/$TanUAngle/2;
    $FirstU_zpos = - $Uactive_z/2 + $FirstUWireOffset/2;

for ($i = 0; $i < $NumberCornerUWires; ++$i)
{

my $ypos = $FirstU_ypos + ($i)*0.5*$UWire_yint;
my $zpos = $FirstU_zpos + ($i)*0.5*$UVReadoutBoardPitch;

# cant actually define like this:
  #    my $ypos = (-0.5*$Uactive_y[$apa]) + $CornerUWireLength[$i]*$CosUAngle/2;
  #    my $zpos = (+0.5*$Uactive_z) - $CornerUWireLength[$i]*$SinUAngle/2;
# since the wire lengths need to be slightly reduced to avoid the 
# wire corner's overlap with the plane boundary

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireU$i${_[3]}"/>
        <position name="pos${_[3]}WireU$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotation name="rUAngle$i"   unit="deg" x="90+$UAngle" y="0"   z="0"/>
      </physvol>
EOF

$topY = $ypos + ($CosUAngle*$CornerUWireLength[$i]/2);
$bottomY = $ypos - ($CosUAngle*$CornerUWireLength[$i]/2);
$edgeZ_p = $zpos + ($SinUAngle*$CornerUWireLength[$i]/2);
$edgeZ_m = $zpos - ($SinUAngle*$CornerUWireLength[$i]/2);
print $wout "U$i: ( $ypos , $zpos ) \n";
print $wout "  -- Y: $bottomY to $topY -- Z: $edgeZ_m to $edgeZ_p \n";

$lastYpos = $ypos;
$lastZpos = $zpos;

}


# Moving upwards to the common wires:
   # x and z are zero to center the wires along a vertical axis
   # y positioning: The trick is positioning the lowest common wire so that the pitch
        # is consistent, then the increment is double the increment of
        # the corner wires since there is no z incriment.
   # rotation: wires in  \\\\  direction, so +90deg to bring them to vertical and 
        # +UAngle counterclockwise to arrive at proper orientation
# Note that the counter maintains wire number (in pos. name) counting bottom to top


my $StartCommonUWires_ypos = $lastYpos + $UWire_yint - abs( $lastZpos )/$TanUAngle;
#print "$StartCommonUWires_ypos = $lastYpos + $UWire_yint - abs( $lastZpos )/$TanUAngle\n";


for ($i = $NumberCornerUWires; $i < $NumberSideUWires; ++$i)
{

    $j = $i - $NumberCornerUWires;

#my $ypos = (-0.5*$Uactive_y[$apa])
#           + 0.5*($NumberCornerUWires)*$UWire_yint+($i+1-$NumberCornerUWires)*$UWire_yint;
    my $ypos = $StartCommonUWires_ypos + $UWire_yint*($j);

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireUCommon${_[3]}"/>
        <position name="pos${_[3]}WireU$i" unit="cm" x="0" y="$ypos " z="0"/>
	<rotation name="rUAngle$i"   unit="deg" x="90+$UAngle" y="0"   z="0"/>
      </physvol>
EOF

$topY = $ypos + ($CosUAngle*$CommonUWireLength/2);
$bottomY = $ypos - ($CosUAngle*$CommonUWireLength/2);
$edgeZ_p = + ($SinUAngle*$CommonUWireLength/2);
$edgeZ_m = - ($SinUAngle*$CommonUWireLength/2);
print $wout "U$i: ( $ypos , 0 ) \n";
print $wout "  -- Y: $bottomY to $topY -- Z: $edgeZ_m to $edgeZ_p \n";

$lastYpos = $ypos;
#$lastZpos = $zpos; always 0

}



my $FirstTopUWire_zspan = $FirstTopUWire_yspan*$TanUAngle;
my $StartTopUWires_ypos =  + $Uactive_y[$apa]/2 - $FirstTopUWire_yspan/2;
my $StartTopUWires_zpos =  + $Uactive_z/2 - $FirstTopUWire_zspan/2;
 
# Finally moving to the corner wires on the top right:
   # x=0 to center the wires in the plane
   # y positioning: plug wire number into same equation
   # z positioning: start at z=0 and go negatively at the same z increment
   # rotation: same as common wire in code above
# note that the counter maintains wire number shown in the position name

for ($j = $NumberSideUWires; $j < $NumberSideUWires+$NumberCornerUWires; ++$j)
{

$i = $j - $NumberSideUWires;

my $ypos = $StartTopUWires_ypos + ($i)*0.5*$UWire_yint;
my $zpos = $StartTopUWires_zpos + ($i)*0.5*$UVReadoutBoardPitch;


print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireU$j${_[3]}"/>
        <position name="pos${_[3]}WireU$j" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotation name="rUAngle$j"   unit="deg" x="90+$UAngle" y="0"   z="0"/>
      </physvol>
EOF

$topY = $ypos + ($CosUAngle*$TopCornerUWireLength[$i]/2);
$bottomY = $ypos - ($CosUAngle*$TopCornerUWireLength[$i]/2);
$edgeZ_p = $zpos + ($SinUAngle*$TopCornerUWireLength[$i]/2);
$edgeZ_m = $zpos - ($SinUAngle*$TopCornerUWireLength[$i]/2);
print $wout "U$i: ( $ypos , $zpos ) \n";
print $wout "  -- Y: $bottomY to $topY -- Z: $edgeZ_m to $edgeZ_p \n";

}

} #ends if wires on


#            ++++++++++++++++++++++  V Plane  +++++++++++++++++++++++

# End U plane and create V plane logical volume
print TPC <<EOF;
    </volume>

    <volume name="volTPCPlaneV${_[3]}">
      <materialref ref="LAr"/>
      <solidref ref="${_[3]}VPlane"/>
EOF

print $wout "\n-     Wires for V plane  -\n\n";
print $wout " Vplane_y: $Vactive_y[$apa]\n";
print $wout " Vplane_z: $Vactive_z\n";

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

    $FirstV_ypos = - $Vactive_y[$apa]/2 + $FirstVWireOffset/$TanVAngle/2;
    $FirstV_zpos = + $Vactive_z/2 - $FirstVWireOffset/2;

for ($i = 0; $i < $NumberCornerVWires; ++$i)
{

my $ypos = $FirstV_ypos + ($i)*0.5*$VWire_yint;
my $zpos = $FirstV_zpos - ($i)*0.5*$UVReadoutBoardPitch;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireV$i${_[3]}"/>
        <position name="pos${_[3]}WireV$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotation name="rVAngle$i"   unit="deg" x="90-$VAngle" y="0"   z="0"/>
      </physvol>
EOF

$topY = $ypos + ($CosVAngle*$CornerVWireLength[$i]/2);
$bottomY = $ypos - ($CosVAngle*$CornerVWireLength[$i]/2);
$edgeZ_p = $zpos + ($SinVAngle*$CornerVWireLength[$i]/2);
$edgeZ_m = $zpos - ($SinVAngle*$CornerVWireLength[$i]/2);
print $wout "V$i: ( $ypos , $zpos ) \n";
print $wout "  -- Y: $bottomY to $topY -- Z: $edgeZ_m to $edgeZ_p \n";

$lastYpos = $ypos;
$lastZpos = $zpos;

}


# Moving upwards to the common wires:
   # x and z are zero to center the wires along a vertical axis
   # y positioning: Plug wire number into the same corner ypos equation
   # rotation: wires in  ////  direction, so +90deg to bring them to vertical and 
        # --VAngle counterclockwise to arrive at proper orientation
# Note that the counter maintains wire number in the position name

my $StartCommonVWires_ypos = $lastYpos + $VWire_yint - abs( $lastZpos )/$TanVAngle;

for ($i = $NumberCornerVWires; $i < $NumberSideVWires; ++$i)
{

    $j = $i - $NumberCornerVWires;
    my $ypos = $StartCommonVWires_ypos + $VWire_yint*($j);

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireVCommon${_[3]}"/>
        <position name="pos${_[3]}WireV$i" unit="cm" x="0" y="$ypos " z="0"/>
        <rotation name="rVAngle$i"   unit="deg" x="90-$VAngle" y="0"   z="0"/>
      </physvol>
EOF

$topY = $ypos + ($CosVAngle*$CommonVWireLength/2);
$bottomY = $ypos - ($CosVAngle*$CommonVWireLength/2);
$edgeZ_p = + ($SinVAngle*$CommonVWireLength/2);
$edgeZ_m = - ($SinVAngle*$CommonVWireLength/2);
print $wout "V$i: ( $ypos , 0 ) \n";
print $wout "  -- Y: $bottomY to $topY -- Z: $edgeZ_m to $edgeZ_p \n";

$lastYpos = $ypos;
#$lastZpos = $zpos; always 0

}


my $FirstTopVWire_zspan = $FirstTopVWire_yspan*$TanVAngle;
my $StartTopVWires_ypos =  + $Vactive_y[$apa]/2 - $FirstTopVWire_yspan/2;
my $StartTopVWires_zpos =  - $Vactive_z/2 + $FirstTopVWire_zspan/2;

# Finally moving to the corner wires on the top right:
   # x=0 to center the wires in the plane
   # y positioning: plug wire number into same equation
   # z positioning: start at z=0 and go positively at the same z increment
   # rotation: same as common wire in code above
# note that the counter maintains wire number shown in the position name

for ($j = $NumberSideVWires; $j < $NumberSideVWires+$NumberCornerVWires; ++$j)
{

$i = $j - $NumberSideVWires;

my $ypos = $StartTopVWires_ypos + ($i)*0.5*$VWire_yint;
my $zpos = $StartTopVWires_zpos - ($i)*0.5*$UVReadoutBoardPitch;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireV$j${_[3]}"/>
        <position name="pos${_[3]}WireV$j" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotation name="rVAngle$j"   unit="deg" x="90-$VAngle" y="0"   z="0"/>
      </physvol>
EOF

$topY = $ypos + ($CosVAngle*$TopCornerVWireLength[$i]/2);
$bottomY = $ypos - ($CosVAngle*$TopCornerVWireLength[$i]/2);
$edgeZ_p = $zpos + ($SinVAngle*$TopCornerVWireLength[$i]/2);
$edgeZ_m = $zpos - ($SinVAngle*$TopCornerVWireLength[$i]/2);
print $wout "V$i: ( $ypos , $zpos ) \n";
print $wout "  -- Y: $bottomY to $topY -- Z: $edgeZ_m to $edgeZ_p \n";

}






} #ends if wires on



#            ++++++++++++++++++++++  Z Plane  +++++++++++++++++++++++

# End V plane and create Z plane logical volume
print TPC <<EOF;
    </volume>

    <volume name="volTPCPlaneZ${_[3]}">
      <materialref ref="LAr"/>
      <solidref ref="${_[3]}ZPlane"/>
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
my $zpos = (-0.5*$Zactive_z) + $i*$XWirePitch + $TPCWireThickness/2;

print TPC <<EOF;
      <physvol>/
        <volumeref ref="volTPCWireVert${_[3]}"/>
        <position name="pos${_[3]}WireZ$i" unit="cm" x="0" y="0 " z="$zpos"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
EOF

}

} #ends if wires on

print TPC <<EOF;
    </volume>
EOF

#+++++++++++++++++++++ ^^ Position physical wires Above ^^ +++++++++++++++++++++

## make the TPC active volume extend down to the G10 for the grid 

    my $BottomOfAPA = - $TPC_y[$apa]/2 + $APAGap_y/2;


if   ($apa==0){
    $BottomOfAPA = -$TPCActive_y[0]/2 + $BottomTallAPAToCage;
    $zz = -1;
} elsif($apa==1){
    $BottomOfAPA = -$TPCActive_y[1]/2 + $APAGap_y/2;
    $zz = 0;  
} elsif($apa==2){ 
    $BottomOfAPA = -$TPCActive_y[2]/2 + $APAGap_y/2;
    $zz = 0;  
} elsif($apa==3){ 
    $BottomOfAPA = -$TPCActive_y[3]/2 + $BottomTallAPAToCage;
    $zz = 1;  
}

if ($TPCActive_x<100){ $xx = -1; }
else                 { $xx = 1;  }

    $posZplane[0]   = $xx*(-$_[0]/2 + $APAWirePlaneSpacing - $TPCWirePlaneThickness/2);
    $posZplane[1]   = $BottomOfAPA + $WrapCover + 4*$G10thickness + $Zactive_y[$apa]/2;
    $posZplane[2]   = $zz*(-$_[2]/2 + $APAGap_z/2 + $APAphys_z/2);

    $posVplane[0]   = $posZplane[0] + $xx*$APAWirePlaneSpacing;
    $posVplane[1]   = $BottomOfAPA + $WrapCover + 3*$G10thickness + $Vactive_y[$apa]/2;
    $posVplane[2]   = $posZplane[2];

    $posUplane[0]   = $posVplane[0] + $xx*$APAWirePlaneSpacing;
    $posUplane[1]   = $BottomOfAPA + $WrapCover + 2*$G10thickness + $Uactive_y[$apa]/2;
    $posUplane[2]   = $posZplane[2];

    $posTPCActive[0] = $posUplane[0] + $xx*($TPCWirePlaneThickness/2 + $TPCActive_x/2);
    $posTPCActive[1] = 0;
    $posTPCActive[2] = 0;

if ($TPCActive_x<100){ $planeRot = "rIdentity"; }
else                 { $planeRot = "rPlus180AboutY";  }


#wrap up the TPC file
print TPC <<EOF;
    <volume name="volTPC${_[3]}">
      <materialref ref="LAr"/>
      <solidref ref="${_[3]}"/>
     <physvol>
       <volumeref ref="volTPCPlaneZ${_[3]}"/>
       <position name="pos${_[3]}PlaneZ" unit="cm" 
         x="$posZplane[0]" y="$posZplane[1]" z="$posZplane[2]"/>
       <rotationref ref="rIdentity"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneV${_[3]}"/>
       <position name="pos${_[3]}PlaneV" unit="cm" 
         x="$posVplane[0]" y="$posVplane[1]" z="$posVplane[2]"/>
       <rotationref ref="$planeRot"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneU${_[3]}"/>
       <position name="pos${_[3]}PlaneU" unit="cm" 
         x="$posUplane[0]" y="$posUplane[1]" z="$posUplane[2]"/>
       <rotationref ref="$planeRot"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCActive${_[3]}"/>
       <position name="pos${_[3]}Active" unit="cm" 
         x="$posTPCActive[0]" y="$posTPCActive[1]" z="$posTPCActive[2]"/>
     </physvol>
EOF

# place the G10 board extensions to the portions placed directly in volCryostat
#place_TPCG10( $_[4],  $_[0],  $_[1],  $_[2] );

print TPC <<EOF;
    </volume>
</structure>
</gdml>
EOF

    close(GDML);

} #end of sub gen_TPC


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ solid_TPCG10 +++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++ vol_TPCG10 +++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ place_TPCG10 +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Must be called only within gen_TPC(), 


# $_[0] = APA number
#   convention: 0 is Largest, 1 is Mid, 2 is Smallest

sub solid_TPCG10()
{

    $apa = $_[0];


$G10BoardYSide_V_x = 2*$APAWirePlaneSpacing; # The rest of the x-direction extension
                                             #   out to the V plane wires
$G10BoardYSide_V_y = $APAFrame_y[$apa];    # Y and Z are the same as the center
$G10BoardYSide_V_z = $G10thickness;          #   parts placed in volCryostat


$G10BoardYSide_U_x = 3*$APAWirePlaneSpacing; # The rest of the x-direction extension
                                             #   out to the U plane wires
$G10BoardYSide_U_y = $APAFrame_y[$apa];    # Y and Z are the same as the center
$G10BoardYSide_U_z = $G10thickness;          #   parts placed in volCryostat


$G10BoardZSide_V_x = 2*$APAWirePlaneSpacing;
$G10BoardZSide_V_y = $G10thickness;
$G10BoardZSide_V_z = $APAFrame_z;
$G10BoardZSide_U_x = 3*$APAWirePlaneSpacing;
$G10BoardZSide_U_y = $G10thickness;
$G10BoardZSide_U_z = $APAFrame_z;
$G10BoardZSide_Grid_x = 4*$APAWirePlaneSpacing;
$G10BoardZSide_Grid_y = $G10thickness;
$G10BoardZSide_Grid_z = $APAFrame_z;

  print TPC <<EOF;

     <box name="G10BoardYSideVSeg\-$apa" lunit="cm"
      x="$G10BoardYSide_V_x"
      y="$G10BoardYSide_V_y"
      z="$G10BoardYSide_V_z"/>

     <box name="G10BoardYSideUSeg\-$apa" lunit="cm"
      x="$G10BoardYSide_U_x"
      y="$G10BoardYSide_U_y"
      z="$G10BoardYSide_U_z"/>

     <box name="G10BoardZSideVSeg\-$apa" lunit="cm"
      x="$G10BoardZSide_V_x"
      y="$G10BoardZSide_V_y"
      z="$G10BoardZSide_V_z"/>

     <box name="G10BoardZSideUSeg\-$apa" lunit="cm"
      x="$G10BoardZSide_U_x"
      y="$G10BoardZSide_U_y"
      z="$G10BoardZSide_U_z"/>

     <box name="G10BoardZSideGridSeg\-$apa" lunit="cm"
      x="$G10BoardZSide_Grid_x"
      y="$G10BoardZSide_Grid_y"
      z="$G10BoardZSide_Grid_z"/>

EOF

}


sub vol_TPCG10()
{

    $apa = $_[0];

  print TPC <<EOF;

    <volume name="volG10BoardYSideVSeg\-$apa">
      <materialref ref="G10"/>
      <solidref ref="G10BoardYSideVSeg\-$apa"/>
    </volume>
    <volume name="volG10BoardYSideUSeg\-$apa">
      <materialref ref="G10"/>
      <solidref ref="G10BoardYSideUSeg\-$apa"/>
    </volume>

    <volume name="volG10BoardZSideVSeg\-$apa">
      <materialref ref="G10"/>
      <solidref ref="G10BoardZSideVSeg\-$apa"/>
    </volume>
    <volume name="volG10BoardZSideUSeg\-$apa">
      <materialref ref="G10"/>
      <solidref ref="G10BoardZSideUSeg\-$apa"/>
    </volume>
    <volume name="volG10BoardZSideGridSeg\-$apa">
      <materialref ref="G10"/>
      <solidref ref="G10BoardZSideGridSeg\-$apa"/>
    </volume>

EOF

}


sub place_TPCG10()
{

    $apa = $_[0];

    $ThisTPC_x = $_[1];
    $ThisTPC_y = $_[2];
    $ThisTPC_z = $_[3];


$G10BoardYSide_V_x = 1*$APAWirePlaneSpacing; # The rest of the x-direction extension
                                             #   out to the V plane wires
$G10BoardYSide_V_y = $APAFrame_y[$APA_i];    # Y and Z are the same as the center
$G10BoardYSide_V_z = $G10thickness;          #   parts placed in volCryostat


$G10BoardYSide_U_x = 2*$APAWirePlaneSpacing; # The rest of the x-direction extension
                                             #   out to the U plane wires
$G10BoardYSide_U_y = $APAFrame_y[$APA_i];    # Y and Z are the same as the center
$G10BoardYSide_U_z = $G10thickness;          #   parts placed in volCryostat


$G10BoardZSide_V_x = 1*$APAWirePlaneSpacing;
$G10BoardZSide_V_y = $G10thickness;
$G10BoardZSide_V_z = $APAFrame_z;
$G10BoardZSide_U_x = 2*$APAWirePlaneSpacing;
$G10BoardZSide_U_y = $G10thickness;
$G10BoardZSide_U_z = $APAFrame_z;
$G10BoardZSide_Grid_x = 3*$APAWirePlaneSpacing;
$G10BoardZSide_Grid_y = $G10thickness;
$G10BoardZSide_Grid_z = $APAFrame_z;


$posG10ZSideZ_y    = $APAFrameCenter_y - $APAFrame_y[$APA_i]/2 - (0+.5)*($G10BoardZSide_y);
$posG10ZSideV_y    = $APAFrameCenter_y - $APAFrame_y[$APA_i]/2 - (1+.5)*($G10BoardZSide_y);
$posG10ZSideU_y    = $APAFrameCenter_y - $APAFrame_y[$APA_i]/2 - (2+.5)*($G10BoardZSide_y);
$posG10ZSideGrid_y = $APAFrameCenter_y - $APAFrame_y[$APA_i]/2 - (3+.5)*($G10BoardZSide_y);

# ... but the smallest APA is upside-down, so put the G10 boards at the top in that case
    if($apa == 2)
    {
	$posG10ZSideZ_y    = $APAFrameCenter_y + $APAFrame_y[$APA_i]/2 + (0+.5)*($G10BoardZSide_y);
	$posG10ZSideV_y    = $APAFrameCenter_y + $APAFrame_y[$APA_i]/2 + (1+.5)*($G10BoardZSide_y);
	$posG10ZSideU_y    = $APAFrameCenter_y + $APAFrame_y[$APA_i]/2 + (2+.5)*($G10BoardZSide_y);
	$posG10ZSideGrid_y = $APAFrameCenter_y + $APAFrame_y[$APA_i]/2 + (3+.5)*($G10BoardZSide_y);
    }


  print TPC <<EOF;

      <!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        - Add the *parts* of the G10 boards that extend those directly in volCryostat.
	- There are two side boards on each the up and downstream end, 
	    one each to wrap the U and V views around the APA frame
	- There are 4 on the bottom which anchor the U V and Z wires and the grid plane
	- The center parts of the G10 boards must be placed directly in volCryostat

	APA $apa

	-->

      <physvol>
        <volumeref ref="volG10BoardYSideVSeg\-$apa"/>
        <position name="posG10BoardYSideVSeg\-up\-$apa" unit="cm" 
	x=" - $ThisTPC_x/2 + $G10BoardYSide_V_x/2" 
	y=" - $APAphys_y[$apa]/2 + $WrapCover + 4*$G10thickness + $APAFrame_y[$apa]/2 " 
	z=" - $APAFrame_z/2 - (0+.5)*($G10BoardYSide_V_z)"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardYSideUSeg\-$apa"/>
        <position name="posG10BoardYSideUSeg\-up\-$apa" unit="cm" 
	x=" - $ThisTPC_x/2 + $G10BoardYSide_U_x/2" 
	y=" - $APAphys_y[$apa]/2 + $WrapCover + 4*$G10thickness + $APAFrame_y[$apa]/2" 
	z=" - $APAFrame_z/2 - (1+.5)*($G10BoardYSide_U_z)"/> 
        <rotationref ref="rIdentity"/>
      </physvol>

      <physvol>
        <volumeref ref="volG10BoardYSideVSeg\-$apa"/>
        <position name="posG10BoardYSideVSeg\-down\-$apa" unit="cm" 
	x=" - $ThisTPC_x/2 + $G10BoardYSide_V_x/2" 
	y=" - $APAphys_y[$apa]/2 + $WrapCover + 4*$G10thickness + $APAFrame_y[$apa]/2 " 
	z=" + $APAFrame_z/2 + (0+.5)*($G10BoardYSide_V_z)"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardYSideUSeg\-$apa"/>
        <position name="posG10BoardYSideUSeg\-down\-$apa" unit="cm" 
	x=" - $ThisTPC_x/2 + $G10BoardYSide_U_x/2" 
	y=" - $APAphys_y[$apa]/2 + $WrapCover + 4*$G10thickness + $APAFrame_y[$apa]/2" 
	z=" + $APAFrame_z/2 + (1+.5)*($G10BoardYSide_U_z)"/> 
        <rotationref ref="rIdentity"/>
      </physvol>



      <physvol>
        <volumeref ref="volG10BoardZSideVSeg\-$apa"/>
        <position name="posG10BoardZSideVSeg\-$apa" unit="cm" 
	x=" - $ThisTPC_x/2 + $G10BoardZSide_V_x/2" 
	y=" - $APAphys_y[$apa]/2 + $WrapCover + 2.5*$G10thickness" 
	z="0"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardZSideUSeg\-$apa"/>
        <position name="posG10BoardZSideUSeg\-$apa" unit="cm" 
	x=" - $ThisTPC_x/2 + $G10BoardZSide_U_x/2" 
	y=" - $APAphys_y[$apa]/2 + $WrapCover + 1.5*$G10thickness" 
	z="0"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardZSideGridSeg\-$apa"/>
        <position name="posG10BoardZSideGridSeg\-$apa" unit="cm" 
	x=" - $ThisTPC_x/2 + $G10BoardZSide_Grid_x/2" 
	y=" - $APAphys_y[$apa]/2 + $WrapCover + 0.5*$G10thickness" 
	z="0"/> 
        <rotationref ref="rIdentity"/>
      </physvol>



EOF

}




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

    <box name="GaseousArgon" lunit="cm" 
      x="$Argon_x"
      y="$HeightGaseousAr"
      z="$Argon_z"/>

    <box name="LightPaddle" lunit="cm"
      x="$LightPaddle_x"
      y="$LightPaddle_y + $SiPM_y"
      z="$LightPaddle_z"/>

<box name="Bar" lunit="cm"
      x="0.6"
      y="50.8"
      z="2.54"/>

  <union name="TwoBars">
   <first ref="Bar"/> <second ref="Bar"/>
   <position name="bs2" z="-2.794" unit="cm"/>
  </union>

  <union name="FourBars">
   <first ref="TwoBars"/> <second ref="TwoBars"/>
   <position name="bs4" z="-5.588" unit="cm"/>
  </union>

    <box name="Plank" lunit="cm"
      x="0.635"
      y="50.8"
      z="11.0"/>

    <box name="Fiber" lunit="cm"
      x="0.289"
      y="50.8"
      z="0.289"/>

    <box name="FiberBottom" lunit="cm"
      x="0.289"
      y="43.18"
      z="0.289"/>

    <box name="FiberTop" lunit="cm"
      x="0.289"
      y="2.2225"
      z="0.289"/>

    <para name="FiberRight" x="0.289" y="0.289" z="5.3975" alpha="0" theta="3.0649" phi="0" aunit="deg" lunit="cm"/>

    <para name="FiberLeft" x="0.289" y="0.289" z="5.3975" alpha="0" theta="-3.0649" phi="0" aunit="deg" lunit="cm"/>

  <union name="FF">
   <first ref="Fiber"/> <second ref="Fiber"/>
   <position name="ff" x="-0.289" z="-0.289" unit="cm"/>
  </union>

  <union name="FFB">
   <first ref="FF"/> <second ref="FiberBottom"/>
   <position name="ffb" x="-0.289" y="-3.81" z="0.289" unit="cm"/>
  </union>

  <union name="FFBB">
   <first ref="FFB"/> <second ref="FiberBottom"/>
   <position name="ffbb" y="-3.81" z="-0.578" unit="cm"/>
  </union>

  <union name="FFBBT">
   <first ref="FFBB"/> <second ref="FiberTop"/>
   <position name="ffbbt" x="-0.289" y="24.28875" unit="cm"/>
  </union>

  <union name="FFBBTT">
   <first ref="FFBBT"/> <second ref="FiberTop"/>
   <position name="ffbbtt" y="24.28875" z="-0.289" unit="cm"/>
  </union>

  <union name="FFBBTTL">
   <first ref="FFBBTT"/> <second ref="FiberLeft"/>
   <position name="ffbbttl" x="0" y="20.47875" z="-0.4335" unit="cm"/> <rotation name="FiberRotLeft" x="-90" y="90" unit="deg"/>
  </union>

  <union name="FourFibers">
   <first ref="FFBBTTL"/> <second ref="FiberRight"/>
   <position name="fs4" x="-0.289" y="20.47875" z="0.1445" unit="cm"/> <rotation name="FiberRotRight" x="-90" y="90" unit="deg"/>
  </union>

  <union name="EightFibers">
   <first ref="FourFibers"/> <second ref="FourFibers"/>
   <position name="fs8" z="-1.2" unit="cm"/>
  </union>

  <union name="SixteenFibers">
   <first ref="EightFibers"/> <second ref="EightFibers"/>
   <position name="fs16" z="-2.4" unit="cm"/>
  </union>

  <union name="ThirtyTwoFibers">
   <first ref="SixteenFibers"/> <second ref="SixteenFibers"/>
   <position name="fs32" z="-4.8" unit="cm"/>
  </union>

EOF


    # Add APA frames with G10 boards on them
    solid_APA( 0 );
    solid_APA( 1 );
    solid_APA( 2 );
    solid_APA( 3 );

    solid_CPA(); # no need for multiple sets of solids, the only 2 CPAs are identical


# Cryostat structure
print CRYO <<EOF;
</solids>
<structure>

    <volume name="volGaseousArgon">
      <materialref ref="ArGas"/>
      <solidref ref="GaseousArgon"/>
    </volume>

    <volume name="volOpDetSensitive_Bar">
      <materialref ref="LAr"/>
      <solidref ref="FourBars"/>
    </volume>

    <volume name="volOpDetSensitive_Plank">
      <materialref ref="LAr"/>
      <solidref ref="Plank"/>
    </volume>

    <volume name="volOpDetSensitive_Fiber">
      <materialref ref="LAr"/>
      <solidref ref="ThirtyTwoFibers"/>
    </volume>

EOF

    vol_APA( 0 );
    vol_APA( 1 );
    vol_APA( 2 );
    vol_APA( 3 );

    vol_CPA();

print CRYO <<EOF;

    <volume name="volCryostat">
      <materialref ref="LAr"/>
      <solidref ref="Cryostat"/>

      <physvol>
        <volumeref ref="volGaseousArgon"/>
        <position name="posGaseousArgon" unit="cm" x="0" y="$Argon_y/2 - $HeightGaseousAr/2" z="0"/> 
      </physvol>

     <physvol>
       <volumeref ref="volCathode"/>
       <positionref ref="posCathodeShortDrift"/>
     </physvol>
     <physvol>
       <volumeref ref="volCathode"/>
       <positionref ref="posCathodeLongDrift"/>
     </physvol>




     <!-- The Smallest APA -->
     <!-- Note this is the only APA with a readout
          at the bottom. TPCs are constructed with
          readouts at the top, so rotate the two APAs
          180 around X to flip upside down             -->

     <physvol>
        <volumeref ref="volTPCSmallestLongDrift"/>
        <positionref ref="posTPCSmallestLongDrift"/>
	<rotationref ref="rPlus180AboutX"/>
      </physvol>

      <physvol>
        <volumeref ref="volTPCSmallestShortDrift"/>
        <positionref ref="posTPCSmallestShortDrift"/>
	<rotationref ref="rPlus180AboutX"/>
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
	<rotationref ref="rIdentity"/>
      </physvol>


EOF

place_APA($APACenter[1][0], $APACenter[1][1], $APACenter[1][2], 1);



print CRYO <<EOF;
      <!-- The Largest APAs, Upstream and Downstream -->

      <physvol>
        <volumeref ref="volTPCLargestLongDriftUpstream"/>
        <positionref ref="posTPCLargestLongDrift_Neg"/>
	<rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPCLargestShortDriftUpstream"/>
        <positionref ref="posTPCLargestShortDrift_Neg"/>
	<rotationref ref="rIdentity"/>
      </physvol>

      <physvol>
        <volumeref ref="volTPCLargestLongDriftDownstream"/>
        <positionref ref="posTPCLargestLongDrift_Pos"/>
	<rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPCLargestShortDriftDownstream"/>
        <positionref ref="posTPCLargestShortDrift_Pos"/>
	<rotationref ref="rIdentity"/>
      </physvol>

EOF

place_APA($APACenter[0][0], $APACenter[0][1], $APACenter[0][2], 0);
place_APA($APACenter[3][0], $APACenter[3][1], $APACenter[3][2], 3);
    
place_CPA(0, $posCPAShortDrift_x, $posCPAShortDrift_y, $posCPAShortDrift_z); 
place_CPA(1, $posCPALongDrift_x, $posCPALongDrift_y, $posCPALongDrift_z); 


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

    $APA_y = $APAFrame_y[$APA_i];

$APAFrameZSide_x = $APAFrame_x;
$APAFrameZSide_y = 4*$inch;
$APAFrameZSide_z = $APAFrame_z;

$APAFrameYSide_x = $APAFrame_x;
$APAFrameYSide_y = $APA_y-2*$APAFrameZSide_y;
$APAFrameYSide_z = 4*$inch;

# Two outer Y supports will sandwich the light paddles
$APAFrameYOuterSupport_x = ($APAFrame_x-$LightPaddle_x)/2;
$APAFrameYOuterSupport_y = $APA_y-2*$APAFrameZSide_y;
$APAFrameYOuterSupport_z = 4*$inch;

$EdgeFrameSteelThickness = 0.12*$inch;
$InnerFrameSteelThickness = 0.062*$inch;


$G10BoardYSide_x = $APAFrame_x;
$G10BoardYSide_y = $APAFrame_y[$APA_i];
$G10BoardYSide_z = $G10thickness;

$G10BoardZSide_x = $APAFrame_x;
$G10BoardZSide_y = $G10thickness;
$G10BoardZSide_z = $APAFrame_z;


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


     <box name="G10BoardYSideCenterSeg\-$APA_i" lunit="cm"
      x="$G10BoardYSide_x"
      y="$G10BoardYSide_y"
      z="$G10BoardYSide_z"/>

     <box name="G10BoardZSideCenterSeg\-$APA_i" lunit="cm"
      x="$G10BoardZSide_x"
      y="$G10BoardZSide_y"
      z="$G10BoardZSide_z"/>

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

    <volume name="volG10BoardYSideCenterSeg\-$APA_i">
      <materialref ref="G10"/>
      <solidref ref="G10BoardYSideCenterSeg\-$APA_i"/>
    </volume>

    <volume name="volG10BoardZSideCenterSeg\-$APA_i">
      <materialref ref="G10"/>
      <solidref ref="G10BoardZSideCenterSeg\-$APA_i"/>
    </volume>


EOF

}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ place_APA ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Must be called only within gen_Cryostat(), 


# $_[0] = x APA physical center
# $_[1] = y APA physical center
# $_[2] = z APA physical center
# $_[3] = APA number
#   convention: 0 and 3 are Largest, 1 is Mid, 2 is Smallest

sub place_APA()
{

    $APA_i = $_[3];

####################################################################
################# APA Frame and Paddle Dimensions ##################

    $APA_y = $APAFrame_y[$APA_i];


# The center passed to this function is the physical APA center,
# which is not quite the frame's center, since there are more boards
# at the bottom. Transform them:

    $APAFrameCenter_x =   $_[0];
    $APAFrameCenter_y =   $_[1] - $APAphys_y[$APA_i]/2 
                        + $WrapCover + 4*$G10thickness
			+ $APAFrame_y[$APA_i]/2;
    $APAFrameCenter_z =   $_[2];

# the smallest APA is different (upside down)
    if($APA_i == 2)
    {
	$APAFrameCenter_y =   $_[1] + $APAphys_y[$APA_i]/2 
                            - $WrapCover - 4*$G10thickness
		     	    - $APAFrame_y[$APA_i]/2;
    }


$APAFrameZSide_x = $APAFrame_x;
$APAFrameZSide_y = 4*$inch;
$APAFrameZSide_z = $APAFrame_z;

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


$G10BoardYSide_x = $APAFrame_x;
$G10BoardYSide_y = $APAFrame_y[$APA_i];
$G10BoardYSide_z = $G10thickness;

$G10BoardZSide_x = $APAFrame_x;
$G10BoardZSide_y = $G10thickness;
$G10BoardZSide_z = $APAFrame_z;

$posG10ZSideZ_y    = $APAFrameCenter_y - $APAFrame_y[$APA_i]/2 - (0+.5)*($G10BoardZSide_y);
$posG10ZSideV_y    = $APAFrameCenter_y - $APAFrame_y[$APA_i]/2 - (1+.5)*($G10BoardZSide_y);
$posG10ZSideU_y    = $APAFrameCenter_y - $APAFrame_y[$APA_i]/2 - (2+.5)*($G10BoardZSide_y);
$posG10ZSideGrid_y = $APAFrameCenter_y - $APAFrame_y[$APA_i]/2 - (3+.5)*($G10BoardZSide_y);

# ... but the smallest APA is upside-down, so put the G10 boards at the top in that case
    if($APA_i == 2)
    {
	$posG10ZSideZ_y    = $APAFrameCenter_y + $APAFrame_y[$APA_i]/2 + (0+.5)*($G10BoardZSide_y);
	$posG10ZSideV_y    = $APAFrameCenter_y + $APAFrame_y[$APA_i]/2 + (1+.5)*($G10BoardZSide_y);
	$posG10ZSideU_y    = $APAFrameCenter_y + $APAFrame_y[$APA_i]/2 + (2+.5)*($G10BoardZSide_y);
	$posG10ZSideGrid_y = $APAFrameCenter_y + $APAFrame_y[$APA_i]/2 + (3+.5)*($G10BoardZSide_y);
    }

   # First put in the frame
  print CRYO <<EOF;

<!--

-->


<!-- 
      <physvol>
        <volumeref ref="volAPAFrameYOuterSupport\-$APA_i"/>
        <position name="posAPAFrameYOuterSupportNeg\-$APA_i" unit="cm" 
	x="$APAFrameCenter_x - ($APAFrameYOuterSupport_x + $APAFrameYInnerSupport_x/2 - $EdgeFrameSteelThickness/2)" 
	y="$APAFrameCenter_y" 
	z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameYOuterSupport\-$APA_i"/>
        <position name="posAPAFrameYOuterSupportPos\-$APA_i" unit="cm" 
	x="$APAFrameCenter_x + ($APAFrameYOuterSupport_x + $APAFrameYInnerSupport_x/2 - $EdgeFrameSteelThickness/2)" 
	y="$APAFrameCenter_y" 
	z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
-->

      <physvol>
        <volumeref ref="volAPAFrameYSide\-$APA_i"/>
        <position name="posAPAFrameYSideNeg\-$_[3]" unit="cm" 
	x="$APAFrameCenter_x" 
	y="$APAFrameCenter_y" 
	z="$APAFrameCenter_z - $APAFrame_z/2 + $APAFrameYSide_z/2"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameYSide\-$APA_i"/>
        <position name="posAPAFrameYSidePos\-$_[3]" unit="cm" 
	x="$APAFrameCenter_x" 
	y="$APAFrameCenter_y" 
	z="$APAFrameCenter_z + $APAFrame_z/2 - $APAFrameYSide_z/2"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameZSide\-$APA_i"/>
        <position name="posAPAFrameZSidePos\-$_[3]" unit="cm" 
	x="$APAFrameCenter_x" 
	y="$APAFrameCenter_y + $APAFrame_y[$APA_i]/2 - $APAFrameZSide_y/2" 
	z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameZSide\-$APA_i"/>
        <position name="posAPAFrameZSideNeg\-$_[3]" unit="cm" 
	x="$APAFrameCenter_x" 
	y="$APAFrameCenter_y  - $APAFrame_y[$APA_i]/2 + $APAFrameZSide_y/2" 
	z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>


      <!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        - Add the *parts* of the G10 boards that exist directly in volCryostat.
	- There are two boards on each the up and downstream end, 
	    one each to wrap the U and V views around the APA frame
	- There are 4 on the bottom which anchor the U V and Z wires and the grid plane
	- The rest of the parts of the G10 boards must be placed directly in volTPC   -->

      <physvol>
        <volumeref ref="volG10BoardYSideCenterSeg\-$APA_i"/>
        <position name="posG10BoardYSideCenterSeg\-Vup\-$APA_i" unit="cm" 
	x="$APAFrameCenter_x" 
	y="$APAFrameCenter_y" 
	z="$APAFrameCenter_z - $APAFrame_z/2 - (0+.5)*($G10BoardYSide_z)"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardYSideCenterSeg\-$APA_i"/>
        <position name="posG10BoardYSideCenterSeg\-Uup\-$APA_i" unit="cm" 
	x="$APAFrameCenter_x" 
	y="$APAFrameCenter_y" 
	z="$APAFrameCenter_z - $APAFrame_z/2 - (1+.5)*($G10BoardYSide_z)"/> 
        <rotationref ref="rIdentity"/>
      </physvol>

      <physvol>
        <volumeref ref="volG10BoardYSideCenterSeg\-$APA_i"/>
        <position name="posG10BoardYSideCenterSeg\-Vdown\-$APA_i" unit="cm" 
	x="$APAFrameCenter_x" 
	y="$APAFrameCenter_y" 
	z="$APAFrameCenter_z + $APAFrame_z/2 + (0+.5)*($G10BoardYSide_z)"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardYSideCenterSeg\-$APA_i"/>
        <position name="posG10BoardYSideCenterSeg\-Udown\-$APA_i" unit="cm" 
	x="$APAFrameCenter_x" 
	y="$APAFrameCenter_y" 
	z="$APAFrameCenter_z + $APAFrame_z/2 + (1+.5)*($G10BoardYSide_z)"/> 
        <rotationref ref="rIdentity"/>
      </physvol>

      <physvol>
        <volumeref ref="volG10BoardZSideCenterSeg\-$APA_i"/>
        <position name="posG10BoardZSideCenterSeg\-Z\-$APA_i" unit="cm" 
	x="$APAFrameCenter_x" 
	y="$posG10ZSideZ_y" 
	z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardZSideCenterSeg\-$APA_i"/>
        <position name="posG10BoardZSideCenterSeg\-V\-$APA_i" unit="cm" 
	x="$APAFrameCenter_x" 
	y="$posG10ZSideV_y" 
	z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardZSideCenterSeg\-$APA_i"/>
        <position name="posG10BoardZSideCenterSeg\-U\-$APA_i" unit="cm" 
	x="$APAFrameCenter_x" 
	y="$posG10ZSideU_y" 
	z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardZSideCenterSeg\-$APA_i"/>
        <position name="posG10BoardZSideCenterSeg\-Grid\-$APA_i" unit="cm" 
	x="$APAFrameCenter_x" 
	y="$posG10ZSideGrid_y"
	z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>



EOF


   # Now loop through paddle y positions and place them
  #for( $p=0; $p<$nPaddlesInAPA[$APA_i]; $p++ ){

#if($nPaddlesInAPA[$APA_i]!=1)
#{
#print CRYO <<EOF;

#     <physvol>
#       <volumeref ref="volOpDetSensitive"/>
#       <position name="posPaddle\-$p\-APA\-$APA_i" unit="cm" 
#         x="$APAFrameCenter_x" 
#	 y="$APAFrameCenter_y - $APAHeight[$APA_i]/2 + $PaddleYPositions[$APA_i][$p]" 
#         z="$APAFrameCenter_z"/>
#       <rotationref ref="rIdentity"/>
#     </physvol>

#EOF


if($APA_i==0){
    place_bar(2);
    place_fiber(2);
    place_bar(3);
} elsif($APA_i==1) {
    place_bar(1);
} elsif($APA_i==2) {
    place_fiber(1);
} elsif($APA_i==3) {
    place_bar(4);
    place_fiber(3);
    place_plank();
}

#switch($APA_i){
#case 0 {	place_bar(2);
#		place_fiber(2);
#		place_bar(3);	}
#case 1 {	place_bar(1);	}
#case 2 {	place_fiber(1);	}
#case 3 {	place_bar(4);
#		place_fiber(3);
#		place_plank();	}
#}


}






#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ solid_CPA +++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++ vol_CPA +++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ place_CPA +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Must be called only within gen_Cryostat(), 

sub solid_CPA()
{

    $apa = $_[0];

print CRYO <<EOF;

    <box name="Cathode" lunit="cm"
      x="$Cathode_x"
      y="$Cathode_y"
      z="$Cathode_z"/>

    <tube name="CPATubeYSide"
      rmin="$CPATube_ID"
      rmax="$CPATube_OD"
      z="$Cathode_y"               
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>

    <tube name="CPATubeZSide"
      rmin="$CPATube_ID"
      rmax="$CPATube_OD"
      z="$Cathode_z"               
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>

EOF

}


sub vol_CPA()
{

  print CRYO <<EOF;

    <volume name="volCathode">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="Cathode"/>
    </volume>
    <volume name="volCPATubeYSide">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="CPATubeYSide"/>
    </volume>
    <volume name="volCPATubeZSide">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="CPATubeZSide"/>
    </volume>


EOF

}


sub place_CPA()
{

    $cpa = $_[0];
    $cpaCenter_x = $_[1];
    $cpaCenter_y = $_[2];
    $cpaCenter_z = $_[3];

print CRYO <<EOF;


      <physvol>
        <volumeref ref="volCathode"/>
        <position name="posCathode\-$cpa" unit="cm" 
	x="$cpaCenter_x - $Cathode_x/2" 
	y="$cpaCenter_y" 
	z="$cpaCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>

<!--
  Curious... The CPATube_OD used here should be CPATube_OD/2,
  but that causes a baffling overlap. look into later! 
-->

      <physvol>
        <volumeref ref="volCPATubeYSide"/>
        <position name="posCPATubeYSideUp\-$cpa" unit="cm" 
	x="$cpaCenter_x" 
	y="$cpaCenter_y" 
	z="$cpaCenter_z - $Cathode_z/2 - $CPATube_OD"/> 
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volCPATubeYSide"/>
        <position name="posCPATubeYSideDown\-$cpa" unit="cm" 
	x="$cpaCenter_x" 
	y="$cpaCenter_y" 
	z="$cpaCenter_z + $Cathode_z/2 + $CPATube_OD"/> 
        <rotationref ref="rPlus90AboutX"/>
      </physvol>

      <physvol>
        <volumeref ref="volCPATubeZSide"/>
        <position name="posCPATubeZSideBottom\-$cpa" unit="cm" 
	x="$cpaCenter_x" 
	y="$cpaCenter_y - $Cathode_y/2 - $CPATube_OD" 
	z="$cpaCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volCPATubeZSide"/>
        <position name="posCPATubeZSideTop\-$cpa" unit="cm" 
	x="$cpaCenter_x" 
	y="$cpaCenter_y + $Cathode_y/2 + $CPATube_OD" 
	z="$cpaCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>


EOF

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



$TopSteelShell_x = $Argon_x - $NeckInside_x;

    $FoamSouth_y = $Argon_y + $NeckInside_y + 2*$SteelShellThickness;
    $FoamTop_x   = $Argon_x - $NeckInside_x;
    $NeckConcreteShell_x = $TotalPadding + 2*$SteelShellThickness + $NeckInside_x + $FoamPadding;

    $TrenchTopConcrete_x = ($DetEnc_x - $CryoWithPadding_x)/2;

# All the detector enclosure solids.
print ENCL <<EOF;
<solids>

    <box name="DetEnclosure" lunit="cm" 
      x="$DetEnc_x"
      y="$DetEnc_y"
      z="$DetEnc_z"/>

    <trd name="AuxDetTrap" lunit="cm"
      x1="$AuxDetScintillatorLongSide" x2="$AuxDetScintillatorShortSide"
      y1="$AuxDetScintillatorThickness" y2="$AuxDetScintillatorThickness"
      z="$AuxDetScintillatorHeight"/>

    <box name="TrenchBottomConcrete" lunit="cm" 
      x="$DetEnc_x"
      y="$DetEnc_y"
      z="$DetEnc_z"/>
    <box name="TrenchBottomConcreteSubtract" lunit="cm" 
      x="$DetEnc_x"
      y="$DetEnc_y -   $TrenchWallThickness"
      z="$DetEnc_z - 2*$TrenchWallThickness"/>
    <subtraction name="TrenchBottomConcreteShell">
      <first ref="TrenchBottomConcrete"/>
      <second ref="TrenchBottomConcreteSubtract"/>
      <position name="posHoleInTrench" unit="cm" x="0" y="$TrenchWallThickness/2" z="0"/> 
    </subtraction>

    <box name="TrenchTopConcrete" lunit="cm" 
      x="$TrenchTopConcrete_x"
      y="$TrenchWallThickness"
      z="$DetEnc_z - 2*$TrenchWallThickness"/>


    <box name="BottomSteel" lunit="cm" 
      x="$Cryostat_x + 2*$SteelShellThickness"
      y="$Cryostat_y +   $SteelShellThickness"
      z="$Cryostat_z + 2*$SteelShellThickness"/>
    <subtraction name="BottomSteelShell">
      <first ref="BottomSteel"/>
      <second ref="Cryostat"/>
      <position name="posLArInShell" unit="cm" x="0" y="$SteelShellThickness/2" z="0"/> 
    </subtraction>


    <box name="Neck" lunit="cm" 
      x="$NeckInside_x + 2*$SteelShellThickness"
      y="$NeckInside_y +   $SteelShellThickness"
      z="$Cryostat_z   + 2*$SteelShellThickness"/>
    <box name="NeckArgon" lunit="cm" 
      x="$NeckInside_x"
      y="$NeckInside_y"
      z="$Cryostat_z"/>
    <subtraction name="NeckSteelShell">
      <first ref="Neck"/>
      <second ref="NeckArgon"/>
      <position name="posLArInNeckShell" unit="cm" x="0" y="-$SteelShellThickness/2" z="0"/> 
    </subtraction>


    <box name="TopSteelShell" lunit="cm" 
      x="$TopSteelShell_x"
      y="$SteelShellThickness"
      z="$Cryostat_z + 2*$SteelShellThickness"/>


    <box name="FoamSouth" lunit="cm" 
      x="$FoamPadding"
      y="$FoamSouth_y"
      z="$Argon_z + 2*$SteelShellThickness + 2*$FoamPadding"/>

    <box name="FoamNorth" lunit="cm" 
      x="$FoamPadding"
      y="$Argon_y + 2*$SteelShellThickness"
      z="$Argon_z + 2*$SteelShellThickness + 2*$FoamPadding"/>

    <box name="FoamEastWest" lunit="cm" 
      x="$Argon_x + 2*$SteelShellThickness"
      y="$Argon_y + 2*$SteelShellThickness"
      z="$FoamPadding"/>

    <box name="FoamEastWestNeck" lunit="cm" 
      x="$NeckInside_x + 2*$SteelShellThickness"
      y="$NeckInside_y + $SteelShellThickness"
      z="$FoamPadding"/>

    <box name="FoamBottom" lunit="cm" 
      x="$Argon_x + 2*$SteelShellThickness + 2*$FoamPadding"
      y="$FoamPadding"
      z="$Argon_z + 2*$SteelShellThickness + 2*$FoamPadding"/>

    <box name="FoamTop" lunit="cm" 
      x="$FoamTop_x"
      y="$FoamPadding"
      z="$Argon_z + 2*$SteelShellThickness + 2*$FoamPadding"/>

    <box name="FoamNorthNeck" lunit="cm" 
      x="$FoamPadding"
      y="$NeckInside_y"
      z="$Argon_z + 2*$SteelShellThickness + 2*$FoamPadding"/>


    <box name="BottomConcrete" lunit="cm" 
      x="$Argon_x + 2*($SteelShellThickness+$FoamPadding) + 2*$ConcretePadding"
      y="$Argon_y + 2*($SteelShellThickness+$FoamPadding) +   $ConcretePadding"
      z="$Argon_z + 2*($SteelShellThickness+$FoamPadding) + 2*$ConcretePadding"/>
    <box name="BottomConcreteSubtract" lunit="cm" 
      x="$Argon_x + 2*($SteelShellThickness+$FoamPadding)"
      y="$Argon_y + 2*($SteelShellThickness+$FoamPadding)"
      z="$Argon_z + 2*($SteelShellThickness+$FoamPadding)"/>
    <subtraction name="BottomConcreteShell">
      <first ref="BottomConcrete"/>
      <second ref="BottomConcreteSubtract"/>
      <position name="posHoleInConcrete" unit="cm" x="0" y="$ConcretePadding/2" z="0"/> 
    </subtraction>


    <box name="NeckConcrete" lunit="cm" 
      x="$NeckConcreteShell_x"
      y="$NeckInside_y +    $SteelShellThickness - $FoamPadding"
      z="$Cryostat_z   + 2*($SteelShellThickness+$FoamPadding) + 2*$ConcretePadding"/>
    <box name="NeckConcreteSubtract" lunit="cm" 
      x="$NeckInside_x + 2*($SteelShellThickness+$FoamPadding)"
      y="$NeckInside_y +    $SteelShellThickness - $FoamPadding"
      z="$Cryostat_z   + 2*($SteelShellThickness+$FoamPadding)"/>
    <subtraction name="NeckConcreteShell">
      <first ref="NeckConcrete"/>
      <second ref="NeckConcreteSubtract"/>
      <position name="posHoleInNeckConcrete" unit="cm" x="-$ConcretePadding/2" y="0" z="0"/> 
    </subtraction>



</solids>
EOF



# Detector enclosure structure
    print ENCL <<EOF;
<structure>

    <volume name="volBotSteelShell">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="BottomSteelShell"/>
    </volume>

    <volume name="volTopSteelShell">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="TopSteelShell"/>
    </volume>

    <volume name="volNeckSteelShell">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="NeckSteelShell"/>
    </volume>

    <volume name="volNeck">
      <materialref ref="ArGas"/>
      <solidref ref="Neck"/>
      <physvol>
        <volumeref ref="volNeckSteelShell"/>
        <position name="posNeckSteelShell" unit="cm" x="0" y="0" z="0"/> 
      </physvol>
    </volume>


    <volume name="volFoamSouth">
      <materialref ref="polyurethane_foam"/>
      <solidref ref="FoamSouth"/>
    </volume>
    <volume name="volFoamNorth">
      <materialref ref="polyurethane_foam"/>
      <solidref ref="FoamNorth"/>
    </volume>
    <volume name="volFoamEastWest">
      <materialref ref="polyurethane_foam"/>
      <solidref ref="FoamEastWest"/>
    </volume>
    <volume name="volFoamEastWestNeck">
      <materialref ref="polyurethane_foam"/>
      <solidref ref="FoamEastWestNeck"/>
    </volume>
    <volume name="volFoamBottom">
      <materialref ref="polyurethane_foam"/>
      <solidref ref="FoamBottom"/>
    </volume>
    <volume name="volFoamTop">
      <materialref ref="polyurethane_foam"/>
      <solidref ref="FoamTop"/>
    </volume>
    <volume name="volFoamNorthNeck">
      <materialref ref="polyurethane_foam"/>
      <solidref ref="FoamNorthNeck"/>
    </volume>

    <volume name="volBottomConcreteShell">
      <materialref ref="Concrete"/>
      <solidref ref="BottomConcreteShell"/>
    </volume>
    <volume name="volNeckConcreteShell">
      <materialref ref="Concrete"/>
      <solidref ref="NeckConcreteShell"/>
    </volume>

    <volume name="volTrenchBottomConcreteShell">
      <materialref ref="Concrete"/>
      <solidref ref="TrenchBottomConcreteShell"/>
    </volume>
    <volume name="volTrenchTopConcrete">
      <materialref ref="Concrete"/>
      <solidref ref="TrenchTopConcrete"/>
    </volume>
EOF

for($i = 1; $i <= 13; ++$i){ 
    print ENCL <<EOF;
    <volume name="volAuxDetBoxBSU-L1-$i" >
      <materialref ref="Acrylic"/>
      <solidref ref="AuxDetBoxBSU"/>
    </volume>
EOF
}
for($i = 1; $i <= 10; ++$i){ 
    print ENCL <<EOF;
    <volume name="volAuxDetBoxBSU-L2-$i" >
      <materialref ref="Acrylic"/>
      <solidref ref="AuxDetBoxBSU"/>
    </volume>
EOF
}
for($i = 1; $i <= 16; ++$i){ 
    print ENCL <<EOF;
    <volume name="volAuxDetBoxBSU-L3-$i" >
      <materialref ref="Acrylic"/>
      <solidref ref="AuxDetBoxBSU"/>
    </volume>
EOF
}
for($i = 1; $i <= 10; ++$i){ 
    print ENCL <<EOF;
    <volume name="volAuxDetBoxBSU-L4-$i" >
      <materialref ref="Acrylic"/>
      <solidref ref="AuxDetBoxBSU"/>
    </volume>
EOF
}

for ($i=1; $i<=12; ++$i) {
    print ENCL <<EOF;
    <volume name="volAuxDetTrap-South-$i" >
      <materialref ref="Acrylic"/>
      <solidref ref="AuxDetTrap"/>
    </volume>
    <volume name="volAuxDetTrap-North-$i" >
      <materialref ref="Acrylic"/>
      <solidref ref="AuxDetTrap"/>
    </volume>
EOF
}
for ($i=1; $i<=10; ++$i) {
    print ENCL <<EOF;
    <volume name="volAuxDetTrap-East-$i" >
      <materialref ref="Acrylic"/>
      <solidref ref="AuxDetTrap"/>
    </volume>
    <volume name="volAuxDetTrap-West-$i" >
      <materialref ref="Acrylic"/>
      <solidref ref="AuxDetTrap"/>
    </volume>
EOF
}

    print ENCL <<EOF;
    <volume name="volDetEnclosure">
      <materialref ref="Air"/>
      <solidref ref="DetEnclosure"/>

      <physvol>
        <volumeref ref="volCryostat"/>
        <position name="posCryo" unit="cm" x="$posCryoInDetEnc_x" y="$posCryoInDetEnc_y" z="$posCryoInDetEnc_z"/>
      </physvol>

      <physvol>
        <volumeref ref="volNeck"/>
        <position name="posNeck" unit="cm" 
	  x="$posCryoInDetEnc_x + $Argon_x/2 - $NeckInside_x/2" 
	  y="$posCryoInDetEnc_y + $Argon_y/2 + ($NeckInside_y + $SteelShellThickness)/2" 
          z="$posCryoInDetEnc_z"/>
      </physvol>
      <physvol>
        <volumeref ref="volTopSteelShell"/>
        <position name="posTopSteelShell" unit="cm" 
	  x="$posCryoInDetEnc_x - $Argon_x/2 - $SteelShellThickness + $TopSteelShell_x/2" 
	  y="$posCryoInDetEnc_y + $Argon_y/2 + $SteelShellThickness/2" 
          z="$posCryoInDetEnc_z"/>
      </physvol>
      <physvol>
        <volumeref ref="volBotSteelShell"/>
        <position name="posBotSteelShell" unit="cm" 
	  x="$posCryoInDetEnc_x"
	  y="$posCryoInDetEnc_y - $SteelShellThickness/2" 
          z="$posCryoInDetEnc_z"/> 
      </physvol>


      <physvol>
        <volumeref ref="volFoamSouth"/>
        <position name="posFoamSouth" unit="cm" 
	  x="$posCryoInDetEnc_x + $Argon_x/2 + $SteelShellThickness + $FoamPadding/2"
	  y="$posCryoInDetEnc_y - $Argon_y/2 - $SteelShellThickness + $FoamSouth_y/2" 
          z="$posCryoInDetEnc_z"/> 
      </physvol>
      <physvol>
        <volumeref ref="volFoamNorth"/>
        <position name="posFoamNorth" unit="cm" 
	  x="$posCryoInDetEnc_x - $Argon_x/2 - $SteelShellThickness - $FoamPadding/2"
	  y="$posCryoInDetEnc_y" 
          z="$posCryoInDetEnc_z"/> 
      </physvol>
      <physvol>
        <volumeref ref="volFoamEastWest"/>
        <position name="posFoamEast" unit="cm" 
	  x="$posCryoInDetEnc_x"
	  y="$posCryoInDetEnc_y" 
          z="$posCryoInDetEnc_z - $Argon_z/2 - $SteelShellThickness - $FoamPadding/2"/> 
      </physvol> 
      <physvol>
        <volumeref ref="volFoamEastWest"/>
        <position name="posFoamWest" unit="cm" 
	  x="$posCryoInDetEnc_x"
	  y="$posCryoInDetEnc_y" 
          z="$posCryoInDetEnc_z + $Argon_z/2 + $SteelShellThickness + $FoamPadding/2"/> 
      </physvol>
      <physvol>
        <volumeref ref="volFoamEastWestNeck"/>
        <position name="posFoamEastNeck" unit="cm" 
	  x="$posCryoInDetEnc_x + $Argon_x/2 - $NeckInside_x/2"
	  y="$posCryoInDetEnc_y + $Argon_y/2 + ($NeckInside_y + $SteelShellThickness)/2" 
          z="$posCryoInDetEnc_z - $Argon_z/2 - $SteelShellThickness - $FoamPadding/2"/> 
      </physvol>
      <physvol>
        <volumeref ref="volFoamEastWestNeck"/>
        <position name="posFoamWestNeck" unit="cm" 
	  x="$posCryoInDetEnc_x + $Argon_x/2 - $NeckInside_x/2"
	  y="$posCryoInDetEnc_y + $Argon_y/2 + ($NeckInside_y + $SteelShellThickness)/2" 
          z="$posCryoInDetEnc_z + $Argon_z/2 + $SteelShellThickness + $FoamPadding/2"/> 
      </physvol>
      <physvol>
        <volumeref ref="volFoamBottom"/>
        <position name="posFoamBottom" unit="cm" 
	  x="$posCryoInDetEnc_x"
	  y="$posCryoInDetEnc_y - $Argon_y/2 - $SteelShellThickness - $FoamPadding/2" 
          z="$posCryoInDetEnc_z"/> 
      </physvol>
      <physvol>
        <volumeref ref="volFoamTop"/>
        <position name="posFoamTop" unit="cm" 
	  x="$posCryoInDetEnc_x - $Argon_x/2 - $SteelShellThickness - $FoamPadding + $FoamTop_x/2"
	  y="$posCryoInDetEnc_y + $Argon_y/2 + $SteelShellThickness + $FoamPadding/2" 
          z="$posCryoInDetEnc_z"/> 
      </physvol>
      <physvol>
        <volumeref ref="volFoamNorthNeck"/>
        <position name="posFoamNorthNeck" unit="cm" 
	  x="$posCryoInDetEnc_x + $Argon_x/2 - $NeckInside_x - $SteelShellThickness - $FoamPadding/2" 
	  y="$posCryoInDetEnc_y + $Argon_y/2 + $SteelShellThickness + $NeckInside_y/2" 
          z="$posCryoInDetEnc_z"/>
      </physvol>


      <physvol>
        <volumeref ref="volBottomConcreteShell"/>
        <position name="posBottomConcreteShell" unit="cm" 
	  x="$posCryoInDetEnc_x"
	  y="$posCryoInDetEnc_y - $ConcretePadding/2" 
          z="$posCryoInDetEnc_z"/> 
      </physvol>
      <physvol>
        <volumeref ref="volNeckConcreteShell"/>
        <position name="posNeckConcreteShell" unit="cm" 
	  x="$posCryoInDetEnc_x + $Argon_x/2 +$SteelShellThickness + $TotalPadding - $NeckConcreteShell_x/2" 
	  y="$posCryoInDetEnc_y + $Argon_y/2 + $FoamPadding + ($NeckInside_y + $SteelShellThickness - $FoamPadding)/2" 
          z="$posCryoInDetEnc_z"/>
      </physvol>



      <physvol>
        <volumeref ref="volTrenchBottomConcreteShell"/>
        <position name="posTrenchConcrete" unit="cm" x="0" y="0" z="0"/>
      </physvol>
<!-- Concrete is too heavy for the roof
      <physvol>
        <volumeref ref="volTrenchTopConcrete"/>
        <position name="posTrenchConcreteTop_Pos" unit="cm" 
	  x=" + $DetEnc_x/2 - $TrenchTopConcrete_x/2" 
	  y="$DetEnc_y/2 - $TrenchWallThickness/2" 
          z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volTrenchTopConcrete"/>
        <position name="posTrenchConcreteTop_Neg" unit="cm" 
	  x=" - $DetEnc_x/2 + $TrenchTopConcrete_x/2" 
	  y="$DetEnc_y/2 - $TrenchWallThickness/2" 
          z="0"/>
      </physvol>
-->


EOF


# place the 12 North Wall AuxDets
#$NorthWallInt_z = ((95+108)/2 - (35+48)/2)*$inch/5;
$NorthWallInt_z = ( 274.32-(274.32-213.36)/4 - (88.9+121.92)/2 )/5;
 # these numbers from a spreadsheet emailed by Ohana Rodrigues
for ($i=1; $i<=12; ++$i)
{
    #$FromCryoCorner_z = (35+48)/2*$inch + (($i-1) % 6)*$NorthWallInt_z;
    $FromCryoCorner_z = (88.9+121.92)/2 + (($i-1) % 6)*$NorthWallInt_z;
    if($i%2==1){
	$rotation = "rAuxDetNSWallUp";
	if($i>6) { $FromCryoCorner_y = (259.08 + 0.12*$inch) + $AuxDetScintillatorHeight/2; }
        else     { $FromCryoCorner_y = (66.04  + 0.12*$inch) + $AuxDetScintillatorHeight/2; }
    }
    else {
	$rotation = "rAuxDetNSWallDown";
	if($i>6) { $FromCryoCorner_y = (325.12 - 0.12*$inch) - $AuxDetScintillatorHeight/2; }
	else     { $FromCryoCorner_y = (132.08 - 0.12*$inch) - $AuxDetScintillatorHeight/2; }
    }

    $posAuxDet_y = $posCryoInDetEnc_y - 375.92/2 + $FromCryoCorner_y;

print ENCL <<EOF;
      <physvol>
        <volumeref ref="volAuxDetTrap-North-$i"/>
        <position name="posAuxDet-North-$i" unit="cm" 
	  x=" $posCryoInDetEnc_x - $CryoWithPadding_x/2 - $AuxDetHousingThickness/2" 
	  y=" $posAuxDet_y"
          z=" $posCryoInDetEnc_z - 411.99/2 + $FromCryoCorner_z"/>
        <rotationref ref="$rotation"/>
      </physvol>
EOF

}

# place the 12 South Wall AuxDets
#$SouthWallInt_z = ((96.5+109.6)/2 - (36.8+49.5)/2)*$inch/5;
$SouthWallInt_z = ( (245.11+278.38)/2 - (93.47+(153.67-93.47)/4) )/5;
 # these numbers from a spreadsheet emailed by Ohana Rodrigues
for ($i=1; $i<=12; ++$i)
{
    #$FromCryoCorner_z = (96.5+109.6)/2*$inch - (($i-1) % 6)*$SouthWallInt_z;
    $FromCryoCorner_z = (245.11+278.38)/2 - (($i-1) % 6)*$SouthWallInt_z;
    if($i%2==1){
	$rotation = "rAuxDetNSWallUp";
	if($i>6) { $FromCryoCorner_y = (272.11 + 0.12*$inch) + $AuxDetScintillatorHeight/2; }
        else     { $FromCryoCorner_y = (31.5   + 0.12*$inch) + $AuxDetScintillatorHeight/2; }
    }
    else {
	$rotation = "rAuxDetNSWallDown";
	if($i>6) { $FromCryoCorner_y = (335.61 - 0.12*$inch) - $AuxDetScintillatorHeight/2; }
	else     { $FromCryoCorner_y = (93.98  - 0.12*$inch) - $AuxDetScintillatorHeight/2; }
    }

    $posAuxDet_y = $posCryoInDetEnc_y - 375.92/2 + $FromCryoCorner_y;

print ENCL <<EOF;
      <physvol>
        <volumeref ref="volAuxDetTrap-South-$i"/>
        <position name="posAuxDet-South-$i" unit="cm" 
	  x=" $posCryoInDetEnc_x + $CryoWithPadding_x/2 + $AuxDetHousingThickness/2" 
	  y=" $posAuxDet_y"
          z=" $posCryoInDetEnc_z - 411.99/2 + $FromCryoCorner_z"/>
        <rotationref ref="$rotation"/>
      </physvol>
EOF

}

# place the 10 East Wall AuxDets
#$EastWallInt_x = ((136.7+149.9)/2 - (28.4+41.4)/2)*$inch/9;
$EastWallInt_x = ( (347.48+380.75)/2 - (72.14+(133.1-72.14)/4) )/9;
 # these numbers from a spreadsheet emailed by Ohana Rodrigues
for ($i=1; $i<=10; ++$i)
{

    # +6cm shift (towards the south wall) as requested by Michelle Stancari
    $FromCryoCorner_x = (347.48+380.75)/2 - ($i-1)*$EastWallInt_x + 6;
    if($i%2==1){ $rotation = "rAuxDetEWWallUp"; 
		 $FromCryoCorner_y = ( 47    + 0.12*$inch) + $AuxDetScintillatorHeight/2; }
    else       { $rotation = "rAuxDetEWWallDown"; 
		 $FromCryoCorner_y = (112.45 - 0.12*$inch) - $AuxDetScintillatorHeight/2; }

    $posAuxDet_y = $posCryoInDetEnc_y - 375.92/2 + $FromCryoCorner_y;

print ENCL <<EOF;
      <physvol>
        <volumeref ref="volAuxDetTrap-East-$i"/>
        <position name="posAuxDet-East-$i" unit="cm" 
	  x=" $posCryoInDetEnc_x - 541.53/2 + $FromCryoCorner_x" 
	  y=" $posAuxDet_y"
          z=" $posCryoInDetEnc_z - $CryoWithPadding_z/2 - $AuxDetHousingThickness/2"/>
        <rotationref ref="$rotation"/>
      </physvol>
EOF

}


# place the 10 West Wall AuxDets
#$WestWallInt_x = ((134+147)/2 - (26.3+39.3)/2)*$inch/9;
$WestWallInt_x = ( (344.42+312.42)/2 - (66.37+99.82)/2 )/8;
 # these numbers from a spreadsheet emailed by Ohana Rodrigues
for ($i=1; $i<=10; ++$i)
{

    # +6cm shift (towards the south wall) as requested by Michelle Stancari
    $FromCryoCorner_x = (66.11+99.82)/2 + ($i-1)*$WestWallInt_x + 6;
    if($i%2==1){ $rotation = "rAuxDetEWWallUp";  
                 $FromCryoCorner_y = (264.44 + 0.12*$inch) + $AuxDetScintillatorHeight/2 }
    else       { $rotation = "rAuxDetEWWallDown";
                 $FromCryoCorner_y = (329.95 - 0.12*$inch) - $AuxDetScintillatorHeight/2 }

    $posAuxDet_y = $posCryoInDetEnc_y - 375.92/2 + $FromCryoCorner_y;

print ENCL <<EOF;
      <physvol>
        <volumeref ref="volAuxDetTrap-West-$i"/>
        <position name="posAuxDet-West-$i" unit="cm" 
	  x=" $posCryoInDetEnc_x - 541.53/2 + $FromCryoCorner_x" 
	  y=" $posAuxDet_y"
          z=" $posCryoInDetEnc_z + $CryoWithPadding_z/2 + $AuxDetHousingThickness/2"/>
        <rotationref ref="$rotation"/>
      </physvol>
EOF

}

# place the 13 Layer 1 telescope counters
for ($i=1; $i<=13; ++$i)
{
print ENCL <<EOF;
      <physvol>
        <volumeref ref="volAuxDetBoxBSU-L1-$i"/>
        <position name="posAuxDet-BSU-L1-$i" unit="cm" 
	  x="$BSULayer1_xpos[$i-1] - $OriginXSet" 
	  y="$BSULayer1_ypos - $OriginYSet"
          z="$BSULayer1_zpos - $OriginZSet"/>
      </physvol>
EOF
}

# place the 9 Layer 2 telescope counters
for ($i=1; $i<=10; ++$i)
{
print ENCL <<EOF;
      <physvol>
        <volumeref ref="volAuxDetBoxBSU-L2-$i"/>
        <position name="posAuxDet-BSU-L2-$i" unit="cm" 
	  x="$BSULayer2_xpos - $OriginXSet" 
	  y="$BSULayer2_ypos - $OriginYSet"
          z="$BSULayer2_zpos[$i-1] - $OriginZSet"/>
        <rotationref ref="rPlus90AboutY"/>
      </physvol>
EOF
}



print ENCL <<EOF;
    </volume>
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

$Dirt_y = $OriginYSet + $World_y/2 + $DetEnc_y/2;

# All the World solids.
print WORLD <<EOF;
<solids>
    <box name="World" lunit="cm" 
      x="$World_x" 
      y="$World_y" 
      z="$World_z"/>

    <box name="AuxDetBoxBSU" lunit="cm" 
      x="$AuxDetBSUScintWidth"
      y="$AuxDetBSUScintHeight"
      z="$AuxDetBSUScintLength"/>

    <box name="DirtBlock" lunit="cm" 
      x="$World_x" 
      y="$Dirt_y" 
      z="$World_z"/>
    <subtraction name="DirtWithHole">
      <first ref="DirtBlock"/>
      <second ref="DetEnclosure"/>
      <position name="posHoleInDirt" unit="cm" 
        x="$OriginXSet" 
        y="$Dirt_y/2 - $DetEnc_y/2" 
        z="$OriginZSet"/> 
    </subtraction>

    <tube name="Berm"
      rmax="$BermRadius"
      z="$DetEnc_x"
      deltaphi="180"
      aunit="deg"
      lunit="cm"/>


</solids>
EOF

# World structure
print WORLD <<EOF;
<structure>

    <volume name="volDirtWithHole" >
      <materialref ref="Dirt"/>
      <solidref ref="DirtWithHole"/>
    </volume>
    <volume name="volBerm" >
      <materialref ref="Dirt"/>
      <solidref ref="Berm"/>
    </volume>

    <volume name="volWorld" >
      <materialref ref="Air"/>
      <solidref ref="World"/>
      <physvol>
        <volumeref ref="volDetEnclosure"/>
        <position name="posDetEnclosure" unit="cm" 
          x="$OriginXSet" 
          y="$OriginYSet" 
          z="$OriginZSet"/>
      </physvol>
      <physvol>
        <volumeref ref="volDirtWithHole"/>
        <position name="posDirtWithHole" unit="cm" 
          x="0" 
          y=" -$World_y/2 + $Dirt_y/2 " 
          z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volBerm"/>
        <position name="posBerm" unit="cm" 
          x="0"
          y="$OriginYSet + $DetEnc_y/2 " 
          z="$OriginZSet + $DetEnc_z/2 + $BermRadius"/>
        <rotationref ref="rPlus90AboutY"/>
      </physvol>

EOF


# place the 16 Layer 3 telescope counters
for ($i=1; $i<=16; ++$i)
{
print WORLD <<EOF;
      <physvol>
        <volumeref ref="volAuxDetBoxBSU-L3-$i"/>
        <position name="posAuxDetBox-BSU-L3-$i" unit="cm" 
	  x="$BSULayer3_xpos[$i-1]" 
	  y="$BSULayer3_ypos"
          z="$BSULayer3_zpos"/>
      </physvol>
EOF
}

# place the 10 Layer 4 telescope counters
for ($i=1; $i<=10; ++$i)
{
print WORLD <<EOF;
      <physvol>
        <volumeref ref="volAuxDetBoxBSU-L4-$i"/>
        <position name="posAuxDetBox-BSU-L4-$i" unit="cm" 
	  x="$BSULayer4_xpos" 
	  y="$BSULayer4_ypos"
          z="$BSULayer4_zpos[$i-1]"/>
        <rotationref ref="rPlus90AboutY"/>
      </physvol>
EOF
}

# place the 16 Layer 5 telescope counters
#for ($i=1; $i<=16; ++$i)
#{
#print WORLD <<EOF;
#<!--
#      <physvol>
#        <volumeref ref="volAuxDetBoxBSU-L5-$i"/>
#        <position name="posAuxDet-BSU-L5-$i" unit="cm" 
#	  x="$BSULayer5_xpos[$i-1]" 
#	  y="$BSULayer5_ypos"
#          z="$BSULayer5_zpos"/>
#      </physvol>
#-->
#EOF
#}

print WORLD <<EOF;

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


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ place_bars +++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub place_bar()
{
print CRYO <<EOF;
	<!-- Bar$_[0]-->
EOF
#	for($j=1; $j<$numberofbars+1; $j++)
#	{

	$pos_ref_name="Bar" . $_[0] . "Pos";
print CRYO <<EOF;
	<physvol>
	<volumeref ref="volOpDetSensitive_Bar"/>
	<positionref ref="$pos_ref_name"/>
	</physvol>

EOF
#	}

}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ place_plank ++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub place_plank()
{
print CRYO <<EOF;
	<!-- Plank-->
EOF

$pos_ref_name="PlankPos";
print CRYO <<EOF;
	<physvol>
	<volumeref ref="volOpDetSensitive_Plank"/>
	<positionref ref="$pos_ref_name"/>
	</physvol>
EOF
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ place_fibers +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub place_fiber()
{
print CRYO <<EOF;
	<!-- 32 Fiber Module $_[0]-->

EOF

############################# Fiber ##############################
$pos_ref_name="Fiber" . $_[0] . "Pos";
print CRYO <<EOF;
	<physvol>
	<volumeref ref="volOpDetSensitive_Fiber"/>
	<positionref ref="$pos_ref_name"/>
	</physvol>

EOF
}
