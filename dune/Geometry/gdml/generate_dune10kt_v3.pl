#!/usr/bin/perl

# contact tylerdalion@gmail.com for any GDML/generate questions
# I would love to help!

# Each subroutine generates a fragment GDML file, and the last subroutine
# creates an XML file that make_gdml.pl will use to appropriately arrange
# the fragment GDML files to create the final desired DUNE GDML file, 
# to be named by make_gdml output command

# If you are playing with different geometries, you can use the
# suffix command to help organize your work.

#use warnings;
#use lib "gdmlMaterials";
use lib '.';
use gdmlMaterials;
#use gdmlMaterials::gen_Materials;
use Math::Trig;
use Getopt::Long;
use Math::BigFloat;
Math::BigFloat->precision(-16);

open(my $wout, '>', 'gdmlWireCenters.txt');

GetOptions( "help|h" => \$help,
          "suffix|s:s" => \$suffix,
          "output|o:s" => \$output,
          "wires|w:s" => \$wires,  
          "protoDune|p:s" => \$protoDune,  
            "workspace|k:s" => \$workspace); 

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

if ( ! defined $workspace )
{
    $workspace = 0;
    if ( ! defined $protoDune )
    {
      $protoDune = 0;
      print "\t\tCreating full geometry.\n";
    }
    elsif ( $protoDune == 1 )
    {
      print "\t\tCreating rough version of protoDUNE.\n";      
    }
}
elsif ( $workspace == 1)
{
    print "\t\tCreating smaller (1x2x2-APA) workspace geometry.\n";
}
elsif ( $workspace == 2)
{
    print "\t\tCreating 1x2x6-APA geometry.\n";
}

# set wires on to be the default, unless given an input by the user
$wires_on = 1; # 1=on, 0=off
if (defined $wires)
{
$wires_on = $wires
}

$tpc_on=1;

$killOuterActive = 0;

# options for different parameters
$Pitch3mmVersion = 0;
$UVAngle45Option = 0;



#++++++++++++++++++++++++ Begin defining variables +++++++++++++++++++++++++

## The GDML is build from the inside out to avoid overlaps, but that means
## that the larger dimensions (such as the cryostat dimensions) are calculated
## to be slightly different than the parameters calculated and reported by 
## engineers. All of this error is within reason and for the better, since
## we don't want this Geometry to be overly complex or long.


$inch = 2.54;

##################################################################
##################### Wire Plane Parameters ######################

# dune10kt
$UWirePitch             =   0.4667;
$VWirePitch             =   0.4667;
$XWirePitch             =   0.479;
$UAngle                 =   35.710;
$VAngle                 =   35.710;
$nZChannelsPerAPA       =   960;

# dune10kt 3mm version
if($Pitch3mmVersion==1){
    $UWirePitch         =   0.3;
    $VWirePitch         =   0.3;
    $XWirePitch         =   0.3;
    $nZChannelsPerAPA   =   2*(229.44/0.3 + 1);
}

# dune10kt ~45 deg UV wires version
if($UVAngle45Option==1){
    $UAngle             =   45.7;
    $VAngle             =   44.3;
}

$nAPAWide              =     3; 
$nAPAHigh              =     2;
$nAPALong              =     25;

if($protoDune==1){
    $nAPAWide              =     2; 
    $nAPAHigh              =     1;
    $nAPALong              =     3;    
}

$SinUAngle              =   sin( deg2rad($UAngle) );
$CosUAngle              =   cos( deg2rad($UAngle) );
$TanUAngle              =   tan( deg2rad($UAngle) );

$SinVAngle              =   sin( deg2rad($VAngle) );
$CosVAngle              =   cos( deg2rad($VAngle) );
$TanVAngle              =   tan( deg2rad($VAngle) );

$UWire_yint             =   $UWirePitch/$SinUAngle;
$UWire_zint             =   $UWirePitch/$CosUAngle;

$VWire_yint             =   $VWirePitch/$SinVAngle;
$VWire_zint             =   $VWirePitch/$CosVAngle;

$TPCWireThickness       =   0.015;

$TPCWirePlaneThickness  =   $TPCWireThickness;
#height and length defined lower



##################################################################
################### Cryostat/APA parameters ######################

if($workspace==1){
    $nAPAWide              =     1; 
    $nAPALong              =     2;
}
elsif($workspace==2){
    $nAPAWide              =     1; 
    $nAPALong              =     6;
}

$nAPAs                 =     $nAPAWide*$nAPAHigh*$nAPALong;


$G10thickness = $inch/8;
$WrapCover    = $inch/16;

$SpaceAPAToCryoWall    =     15;
$SpaceAPAToFloor       =     59.4  - 10.2; # 10.2cm just adjusts for how spaces are reported
$SpaceAPAToTopLAr      =     50.9  - 10.2;
$UpstreamLArPadding    =     311.4 - 10.2;
$DownstreamLArPadding  =     59.4  - 10.2;

if($protoDune==1){
    $SpaceAPAToCryoWall    =     20;
    $SpaceAPAToFloor       =     45;
    $SpaceAPAToTopLAr      =     45;
    $UpstreamLArPadding    =     80;
    $DownstreamLArPadding  =     80;
}

 #InnerDrift is the max distance form the edge of the CPA to the edge of the first wire plane
$InnerDrift              = 359.4;
$OuterDrift              = $SpaceAPAToCryoWall;
$APAFrame_x              = 5.0661; # ~2in -- this does not include the wire spacing

$TPCWireThickness        = 0.015;
$TPCWirePlaneThickness   = $TPCWireThickness;
#$APAWirePlaneSpacing     = 0.4730488 + $TPCWirePlaneThickness; # center to center spacing between all of the wire planes (g, u, v, and x)

if($Pitch3mmVersion==1){
    $APAWirePlaneSpacing     = 0.3;
} else {
    $APAWirePlaneSpacing     = 0.476;
}

# At creation of the plane volumes, the y and z boundaries will be increased
# by this much at each of the 4 edges. this is so the corners of the wire 
# tubes don't extrude. For all other purposes, the plane dimensions stay as originally defined
$UVPlaneBoundNudge = $TPCWireThickness;


# APA z-dimensions:
# The following are all widths about the z center of an APA
$Zactive_z    = ($nZChannelsPerAPA/2-1)*$XWirePitch + $TPCWireThickness;
$APAFrame_z   = 231.59 - 2*(2*$G10thickness+$WrapCover);
$Vactive_z    = $APAFrame_z;
$Uactive_z    = $APAFrame_z + 2*$G10thickness;
$APAphys_z    = 231.59; 


# APA y-dimensions:
$ReadoutBoardOverlap = 7.61; #board overlaps wires, chop this off of their active height
$APAFrame_y = 606.0;
$Zactive_y    = $APAFrame_y + 0*$G10thickness - $ReadoutBoardOverlap;
$Vactive_y    = $APAFrame_y + 1*$G10thickness - $ReadoutBoardOverlap; 
$Uactive_y    = $APAFrame_y + 2*$G10thickness - $ReadoutBoardOverlap;
  # the last G10 board for the grid, then a cover. This is not "covered" by the board
$APAphys_y    = $APAFrame_y + 4*$G10thickness + $WrapCover; 

$APAGap_y     =    0.4;  #separation between APAs (cover to cover) along the incident beam axis
$APAGap_z     =    0.8;  #separation between APAs (cover to cover) along the vertical axis


  # include APA spacing in y and z so volTPCs touch in y and z directions with correct APA
  # spacing - this makes for smoother event generation. 

$APA_UtoU_x   = $APAFrame_x + 6*$APAWirePlaneSpacing + $TPCWirePlaneThickness; # thickness of APA between its 
                                                                               # outer wire planes (center to center)
$TPCInner_x   = $InnerDrift   + 3*$APAWirePlaneSpacing + $TPCWirePlaneThickness;
$TPCOuter_x   = $OuterDrift + 3*$APAWirePlaneSpacing + $TPCWirePlaneThickness;

$TPC_z    =   $APAphys_z + $APAGap_z;
$TPC_y    =   $APAphys_y + $APAGap_y;

$CPATube_OD = 5.066;
#$CPATube_ID = 4.747;

$Cathode_x                 =    0.016;
$Cathode_y                 =    $APAphys_y - $CPATube_OD;
$Cathode_z                 =    $APAphys_z - $CPATube_OD;   


$APAToAPA =   $APAFrame_x 
            + 2*$TPCInner_x
            + $Cathode_x; # center to center

$CPAToAPA =   $APAFrame_x/2
            + $TPCInner_x
            + $Cathode_x/2; # center to center

$SteelThickness            =      0.5*$inch; #half inch
$HeightGaseousAr        =       50;

$Argon_x            =       ($nAPAWide-1)*$APAToAPA
                              + $APA_UtoU_x + $TPCWirePlaneThickness
                              + 2*$SpaceAPAToCryoWall;


#FIELD CAGE IN (FCI)                        
$FC_x                   =       ($nAPAWide-1)*$APAToAPA
                              + $APA_UtoU_x + $TPCWirePlaneThickness;

if($workspace==1||$workspace==2){ # this is arbitrary for the workspace, but size down a little
$Argon_x = 2*$CPAToAPA + $Cathode_x + 2*$SpaceAPAToCryoWall;
$FC_x    = 2*$CPAToAPA + $Cathode_x;
}


$LArgon_y              =      $nAPAHigh*$APAphys_y 
                                + ($nAPAHigh-1)*$APAGap_y
                                + $SpaceAPAToFloor + $SpaceAPAToTopLAr;
$Argon_y               =      $LArgon_y                        
                                + $HeightGaseousAr; 
               
$Argon_z               =      $nAPALong*$APAphys_z 
                                + ($nAPALong-1)*$APAGap_z
                                + $UpstreamLArPadding + $DownstreamLArPadding;

$FCI_y                 =      $nAPAHigh*$TPC_y+$FC_b; #The +FC_b is an arbitrary +2cm to make sure the FC doesn't contact the APA or CPAs.
$FCO_y                 =      $FCI_y+(0.1*$inch)*2; #The outer boundary of the Fied Cage
$FCI_z                 =      $nAPALong*$TPC_z
                                + $FC_b; #the plus FC_b is +2 cm arbitrarily. This makes sure the FC doesn't contact the APAs.
$FCO_z                 =      $FCI_z+(0.1*$inch)*2;                        

$FC_b                  =      2; #2 cm buffer for FC

$Cryostat_x            =      $Argon_x  + 2*$SteelThickness;
$Cryostat_y            =      $Argon_y + 2*$SteelThickness;
$Cryostat_z            =      $Argon_z + 2*$SteelThickness;

### Add dimensions for the source cylinder.
$Source_x_center              =        $FC_x_center + 220; # 220 cm past the center plane of the TPCs.
$Source_y_center              =        $FC_y_center + 300 -20 ; # 300 cm above the center of the active volume (-20 unexplained correction taken from viewer.
$Source_z_center              =        $FC_z_center - ($FCO_z/2) - 40+ 130; # 40 cm outside the field cage in Z. (+130 unexplained correction from viewer.
$SourceModerator_d              =      20; #20 cm diameter cylinder
$SourceModerator_z              =      20; # 20 cm tall cylinder  

##TO DO.
#Make Source VolumeReference
#Make Source position reference.
#Make SourceModerator
#Make Source child of Source Moderator

##################################################################
############## DetEnc and World relevant parameters  #############


$SteelSupport_x        =      100;
$SteelSupport_y        =      50;
$SteelSupport_z        =      100;
$FoamPadding           =        80;
$FracMassOfSteel       =        0.998; #The steel support is not a solid block, but a mixture of air and steel
$FracMassOfAir         =        1-$FracMassOfSteel;

$SpaceSteelSupportToWall    = 100;
$SpaceSteelSupportToCeiling = 100;

$DetEncWidth              =      $Cryostat_x
                              + 2*($SteelSupport_x + $FoamPadding) + $SpaceSteelSupportToCeiling;
$DetEncHeight              =      $Cryostat_y
                              + 2*($SteelSupport_y + $FoamPadding) + 2*$SpaceSteelSupportToWall;
$DetEncLength           =       $Cryostat_z
                              + 2*($SteelSupport_z + $FoamPadding) + 2*$SpaceSteelSupportToWall;

$posCryoInDetEnc_y = - $DetEncHeight/2 + $SteelSupport_y + $FoamPadding + $Cryostat_y/2;

$RockThickness              =       3000;
if ($workspace == 0){
    # Initially, the world dimensions and the OriginZSet
    #   left only ~222cm to the world boundary from the 
    #   +z boundary of volDetEnclosure. Bump that up to 
    #   at least 10m for good measure.
    # This is in an if statement so that it does not affect
    #   the workspace geometries, which are already being used
    $RockThickness += 800
}

  # We want the world origin to be at the very front of the fiducial volume.
  # move it to the front of the enclosure, then back it up through the concrete/foam, 
  # then through the Cryostat shell, then through the upstream dead LAr (including the
  # dead LAr on the edge of the TPC, but this is covered in $UpstreamLArPadding).
  # This is to be added to the z position of every volume in volWorld

$OriginZSet             =       $DetEncLength/2 
                              - $SpaceSteelSupportToWall
                              - $SteelSupport_z
                              - $FoamPadding
                              - $SteelThickness
                              - $UpstreamLArPadding
                              - ($APAphys_z - $Uactive_z)/2;

  # We want the world origin to be vertically centered between the stacked APAs.
  # This is to be added to the y position of every volume in volWorld

$OriginYSet             =       $DetEncHeight/2
                              - $SteelSupport_y
                              - $FoamPadding
                              - $SteelThickness
                              - $SpaceAPAToFloor
                              - $APAphys_y
                              - $APAGap_y/2;

if($protoDune==1){ # bring the origin to the bottom of the APAs for protoDUNE
    $OriginYSet = $OriginYSet + $APAphys_y + $APAGap_y/2;
}

$OriginXSet             =     0; # centered for now






####################################################################
########################  Paddle Dimensions ########################

$APAFrameZSide_y = 4*$inch;
$APAFrameYSide_z = 4*$inch;
$LightPaddle_x            =    0.476;
$LightPaddle_y            =    4*$inch;
$LightPaddle_z            =    $APAFrame_z - 2*$APAFrameYSide_z;
$nLightPaddlesPerAPA      =    10;   # 10, or 20 for double coverage (for now)
$PaddleYInterval          =    (2*$APAphys_y+$APAGap_y-$LightPaddle_y-2*$APAFrameZSide_y)
                              /(2*$nLightPaddlesPerAPA-1);
$FrameToPaddleSpace       =    ($PaddleYInterval-$APAGap_y)/2;

$SiPM_z                   =    0;

# $PaddleYInterval is defined so that the center-to-center distance in the 
# y direction between paddles is uniform between vertically stacked APAs.
# $FrameToPaddleSpace is from the BOTTOM of the APA frame (4" in y direction)
# to the CENTER of a paddle, including the 4" part of the frame. This variable's
# primary purpose is to position the lowest paddle in each APA.




#+++++++++++++++++++++++++ End defining variables ++++++++++++++++++++++++++


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
    $DEF = "dune10kt_v3_Def" . $suffix . ".gdml";
    push (@gdmlFiles, $DEF);
    $DEF = ">" . $DEF;
    open(DEF) or die("Could not open file $DEF for writing");


print DEF <<EOF;
<?xml version='1.0'?>
<gdml>
<define>

<!--



-->

   <position name="posCryoInDetEnc"     unit="cm" x="0" y="@{[$posCryoInDetEnc_y]}" z="0"/>
   <position name="posCenter"           unit="cm" x="0" y="0" z="0"/>
   <rotation name="rPlus90AboutX"       unit="deg" x="90" y="0" z="0"/>
   <rotation name="rMinus90AboutY"      unit="deg" x="0" y="270" z="0"/>
   <rotation name="rMinus90AboutYMinus90AboutX"       unit="deg" x="270" y="270" z="0"/>
   <rotation name="rPlusUAngleAboutX"      unit="deg" x="@{[90-$UAngle]}" y="0"   z="0"/>
   <rotation name="rPlusVAngleAboutX"      unit="deg" x="@{[90+$VAngle]}" y="0"   z="0"/>
   <rotation name="rPlus180AboutX"      unit="deg" x="180" y="0"   z="0"/>
   <rotation name="rPlus180AboutY"      unit="deg" x="0" y="180"   z="0"/>
   <rotation name="rPlus180AboutXPlus180AboutY"      unit="deg" x="180" y="180"   z="0"/>
   <rotation name="rIdentity"            unit="deg" x="0" y="0"   z="0"/>
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
    $MAT = "dune10kt_v3_Materials" . $suffix . ".gdml";
    push (@gdmlFiles, $MAT);
    $MAT = ">" . $MAT;

    open(MAT) or die("Could not open file $MAT for writing");

    # Add any materials special to this geometry by defining a mulitline string
    # and passing it to the gdmlMaterials::gen_Materials() function.

    $DensityAirSteel = 1/(0.001205/$FracMassOfAir + 7.9300/$FracMassOfSteel);

my $asmix = <<EOF;
  <!-- preliminary values -->
  <material name="AirSteelMixture" formula="AirSteelMixture">
   <D value="$DensityAirSteel" unit="g/cm3"/>
   <fraction n="$FracMassOfSteel" ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
   <fraction n="$FracMassOfAir"   ref="Air"/>
  </material>
EOF

    # add the general materials used anywere
    print MAT gdmlMaterials::gen_Materials( $asmix );

    close(MAT);
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++ gen_TPC ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


sub gen_TPC
{

# $_[0] = $TPC_x
# $_[1] = $TPC_y
# $_[2] = $TPC_z
# $_[3] = 'name'

    my $TPCActive_x   =  $_[0]-(3*$APAWirePlaneSpacing);
    my $TPCActive_y   =  $_[1] - $APAGap_y/2 - $ReadoutBoardOverlap ; #TODO: make the Active height more accurate
    my $TPCActive_z   =  $_[2];


#constructs everything inside volTPC, namely
# (moving from left to right, or from +x to -x)
#  -volCPActive
#  -volTPCPlaneU: with wires angled from vertical slightly different than in V
#  -volTPCPlaneV: with wires angled from vertical slightly differently than in U
#  -volTPCPlaneX: with vertical wires


# Create the TPC fragment file name,
# add file to list of output GDML fragments,
# and open it
    $TPC = "dune10kt_v3_TPC_${_[3]}" . $suffix . ".gdml";
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
      x="@{[$_[0]]}" 
      y="@{[$_[1]]}" 
      z="@{[$_[2]]}"/>
    <box name="${_[3]}UPlane" lunit="cm" 
      x="@{[$TPCWirePlaneThickness]}" 
      y="@{[$Uactive_y + $UVPlaneBoundNudge]}" 
      z="@{[$Uactive_z + $UVPlaneBoundNudge]}"/>
    <box name="${_[3]}VPlane" lunit="cm" 
      x="@{[$TPCWirePlaneThickness]}" 
      y="@{[$Vactive_y + $UVPlaneBoundNudge]}" 
      z="@{[$Vactive_z + $UVPlaneBoundNudge]}"/>
    <box name="${_[3]}ZPlane" lunit="cm" 
      x="@{[$TPCWirePlaneThickness]}" 
      y="@{[$Zactive_y]}" 
      z="@{[$Zactive_z]}"/>
    <box name="${_[3]}Active" lunit="cm"
      x="@{[$TPCActive_x]}"
      y="@{[$TPCActive_y]}"
      z="@{[$TPCActive_z]}"/>
EOF


#++++++++++++++++++++++++++++ Wire Solids ++++++++++++++++++++++++++++++

print TPC <<EOF;

    <tube name="${_[3]}WireVert"
      rmax="@{[0.5*$TPCWireThickness]}"
      z="@{[$Zactive_y]}"               
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
$NumberCornerUWires = int( $APAFrame_z/($UWirePitch/$CosUAngle) );

$NumberCornerVWires = int( $APAFrame_z/($VWirePitch/$CosVAngle) );


   # Total number of wires touching one vertical (longer) side
   # Note that the total number of wires per plane is this + another set of corner wires
$NumberSideUWires = int( $Uactive_y/$UWire_yint );
if($Pitch3mmVersion==1){ $NumberSideUWires = $NumberSideUWires-1; }

$NumberSideVWires = int( $Vactive_y/$VWire_yint );

   # Number of wires per side that aren't cut off by the corner
$NumberCommonUWires = $NumberSideUWires - $NumberCornerUWires;

$NumberCommonVWires = $NumberSideVWires - $NumberCornerVWires;

   # Number of wires on the vertical plane
   #   Since APA Active z is defined in docdb 3383 to be distance 
   #   between outer vertical wires, + 1 since the floor of this
   #   division will be one under, giving the amt of spaces, not wires
$NumberVerticalWires = int( $Zactive_z/$XWirePitch ) + 1;
#$NumberVerticalWires = 960; 


$nUchans = 2*$NumberCornerUWires;
$nVchans = 2*$NumberCornerVWires;

print $wout "$nUchans U channels\n";
print $wout "$nVchans V channels\n";
print $wout "$NumberVerticalWires Z channels per side\n";

print $wout "$NumberCornerUWires U corner wires\n";
print $wout "$NumberCommonUWires U common wires\n";
print $wout "$NumberCornerVWires V corner wires\n";
print $wout "$NumberCommonVWires V common wires\n";

}

# hard codeed number will be a factor determined from engineering spreadsheets on wire endpoints,
# but since that won't exist for a while, use this number to avoid overlaps
my $FirstUWireOffset = .55 + $G10thickness + 2*$G10thickness*$TanUAngle - $UWire_zint;
my $FirstVWireOffset = .5; # doesnt include a G10 board in width

if($Pitch3mmVersion==1){
    my $FirstUWireOffset = .15 + $G10thickness + 2*$G10thickness*$TanUAngle - $UWire_zint;
    my $FirstVWireOffset = .15; # doesnt include a G10 board in width
}

if($UVAngle45Option==1){$FirstVWireOffset = .7;}


my $FirstTopUWire_yspan =
    $Uactive_y/2
    - ( - $Uactive_y/2
        + $FirstUWireOffset/$TanUAngle      # walk us up to the first wire
        + $UWire_yint*($NumberSideUWires-1) # up to the top of the top common wire
        - $Uactive_z/$TanUAngle             # back to the bottom of the top common wire
      + $UWire_yint);                     # nudge up to bottom of the first top corner wire

my $FirstTopVWire_yspan =
    $Vactive_y/2
    - ( - $Vactive_y/2
        + $FirstVWireOffset/$TanVAngle      # walk us up to the first wire
        + $VWire_yint*($NumberSideVWires-1) # up to the top of the top common wire
        - $Vactive_z/$TanVAngle             # back to the bottom of the top common wire
      + $VWire_yint);                     # nudge up to bottom of the first top corner wire


# The corner wires for the U plane
if ($wires_on==1)
{
    for ($i = 0; $i < $NumberCornerUWires; $i++)
    {
      $CornerUWireLength[$i] = ($FirstUWireOffset + $i*$UWire_zint)/$SinUAngle;

   print TPC <<EOF;
    <tube name="${_[3]}WireU$i"
      rmax="@{[0.5*$TPCWireThickness]}"
      z="@{[$CornerUWireLength[$i]]}"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

    }

    $CommonUWireLength = $Uactive_z/$SinUAngle;

   print TPC <<EOF;
    <tube name="${_[3]}WireUCommon"
      rmax="@{[0.5*$TPCWireThickness]}"
      z="@{[$CommonUWireLength]}"
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
      rmax="@{[0.5*$TPCWireThickness]}"
      z="@{[$TopCornerUWireLength[$i]]}"
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
      $CornerVWireLength[$i] = ($FirstVWireOffset + $i*$VWire_zint)/$SinVAngle;

      print TPC <<EOF;

    <tube name="${_[3]}WireV$i"
      rmax="@{[0.5*$TPCWireThickness]}"
      z="@{[$CornerVWireLength[$i]]}"
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
      rmax="@{[0.5*$TPCWireThickness]}"
      z="@{[$CommonVWireLength]}"
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
      rmax="@{[0.5*$TPCWireThickness]}"
      z="@{[$TopCornerVWireLength[$i]]}"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

    }


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
print $wout " Uplane_y: $Uactive_y\n";
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

    $FirstU_ypos = - $Uactive_y/2 + $FirstUWireOffset/$TanUAngle/2;
    $FirstU_zpos = + $Uactive_z/2 - $FirstUWireOffset/2;

for ($i = 0; $i < $NumberCornerUWires; ++$i)
{

my $ypos = $FirstU_ypos + ($i)*0.5*$UWire_yint;
my $zpos = $FirstU_zpos - ($i)*0.5*$UWire_zint;

$pitch  =  ($ypos - $lastYpos) * $SinUAngle
         - ($zpos - $lastZpos) * $CosUAngle;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireU$i${_[3]}"/>
        <position name="pos${_[3]}WireU$i" unit="cm" x="0" y="@{[$ypos ]}" z="@{[$zpos]}"/>
        <rotation name="rUAngle$i"   unit="deg" x="@{[90-$UAngle]}" y="0"   z="0"/>
      </physvol>
EOF

$topY = $ypos + ($CosUAngle*$CornerUWireLength[$i]/2);
$bottomY = $ypos - ($CosUAngle*$CornerUWireLength[$i]/2);
$edgeZ_p = $zpos + ($SinUAngle*$CornerUWireLength[$i]/2);
$edgeZ_m = $zpos - ($SinUAngle*$CornerUWireLength[$i]/2);
print $wout "U$i: ( $ypos , $zpos ) (pitch = $pitch)\n";
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

for ($i = $NumberCornerUWires; $i < $NumberSideUWires; ++$i)
{

    $j = $i - $NumberCornerUWires;
    my $ypos = $StartCommonUWires_ypos + $UWire_yint*($j);

    $lastWnum = $i-1;
    if ( $ypos <  $lastYpos ){  print "WARNING: y position dropped from $lastYpos (wire U$lastWnum) to $ypos (wire U$i)\n"; }
    if ( $ypos == $lastYpos ){  print "WARNING: y position between wire U$lastWnum and U$lastWnum did not move: $ypos\n"; }


$pitch  =   ($ypos - $lastYpos) * $SinUAngle - ($zpos - $lastZpos) * $CosUAngle ;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireUCommon${_[3]}"/>
        <position name="pos${_[3]}WireU$i" unit="cm" x="0" y="@{[$ypos ]}" z="0"/>
      <rotation name="rUAngle$i"   unit="deg" x="@{[90-$UAngle]}" y="0"   z="0"/>
      </physvol>
EOF

$topY = $ypos + ($CosUAngle*$CommonUWireLength/2);
$bottomY = $ypos - ($CosUAngle*$CommonUWireLength/2);
$edgeZ_p = + ($SinUAngle*$CommonUWireLength/2);
$edgeZ_m = - ($SinUAngle*$CommonUWireLength/2);
print $wout "U$i: ( $ypos , 0 ) (pitch = $pitch)\n";
print $wout "  -- Y: $bottomY to $topY -- Z: $edgeZ_m to $edgeZ_p \n";

$lastYpos = $ypos;
$lastZpos = 0;

}



my $FirstTopUWire_zspan = $FirstTopUWire_yspan*$TanUAngle;
my $StartTopUWires_ypos =  + $Uactive_y/2 - $FirstTopUWire_yspan/2;
my $StartTopUWires_zpos =  - $Uactive_z/2 + $FirstTopUWire_zspan/2;
 
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
my $zpos = $StartTopUWires_zpos - ($i)*0.5*$UWire_zint;

    $lastWnum = $j-1;
    if ( $ypos <  $lastYpos ){  print "WARNING: y position dropped from $lastYpos (wire U$lastWnum) to $ypos (wire U$j)\n"; }
    if ( $ypos == $lastYpos ){  print "WARNING: y position between wire U$lastWnum and U$j did not move: $ypos\n"; }

$pitch  =   ($ypos - $lastYpos) * $SinUAngle - ($zpos - $lastZpos) * $CosUAngle ;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireU$j${_[3]}"/>
        <position name="pos${_[3]}WireU$j" unit="cm" x="0" y="@{[$ypos ]}" z="@{[$zpos]}"/>
        <rotation name="rUAngle$j"   unit="deg" x="@{[90-$UAngle]}" y="0"   z="0"/>
      </physvol>
EOF

$topY = $ypos + ($CosUAngle*$TopCornerUWireLength[$i]/2);
$bottomY = $ypos - ($CosUAngle*$TopCornerUWireLength[$i]/2);
$edgeZ_p = $zpos + ($SinUAngle*$TopCornerUWireLength[$i]/2);
$edgeZ_m = $zpos - ($SinUAngle*$TopCornerUWireLength[$i]/2);
print $wout "U$j: ( $ypos , $zpos ) (pitch = $pitch)\n";
print $wout "  -- Y: $bottomY to $topY -- Z: $edgeZ_m to $edgeZ_p \n";

$lastYpos = $ypos;
$lastZpos = $zpos;

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
print $wout " Vplane_y: $Vactive_y\n";
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

    $FirstV_ypos = - $Vactive_y/2 + $FirstVWireOffset/$TanVAngle/2;
    $FirstV_zpos = - $Vactive_z/2 + $FirstVWireOffset/2;

for ($i = 0; $i < $NumberCornerVWires; ++$i)
{

my $ypos = $FirstV_ypos + ($i)*0.5*$VWire_yint;
my $zpos = $FirstV_zpos + ($i)*0.5*$VWire_zint;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireV$i${_[3]}"/>
        <position name="pos${_[3]}WireV$i" unit="cm" x="0" y="@{[$ypos ]}" z="@{[$zpos]}"/>
        <rotation name="rVAngle$i"   unit="deg" x="@{[90+$VAngle]}" y="0"   z="0"/>
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
        <position name="pos${_[3]}WireV$i" unit="cm" x="0" y="@{[$ypos ]}" z="0"/>
        <rotation name="rVAngle$i"   unit="deg" x="@{[90+$VAngle]}" y="0"   z="0"/>
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
my $StartTopVWires_ypos =  + $Vactive_y/2 - $FirstTopVWire_yspan/2;
my $StartTopVWires_zpos =  + $Vactive_z/2 - $FirstTopVWire_zspan/2;

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
my $zpos = $StartTopVWires_zpos + ($i)*0.5*$VWire_zint;

print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireV$j${_[3]}"/>
        <position name="pos${_[3]}WireV$j" unit="cm" x="0" y="@{[$ypos ]}" z="@{[$zpos]}"/>
        <rotation name="rVAngle$j"   unit="deg" x="@{[90+$VAngle]}" y="0"   z="0"/>
      </physvol>
EOF

$topY = $ypos + ($CosVAngle*$TopCornerVWireLength[$i]/2);
$bottomY = $ypos - ($CosVAngle*$TopCornerVWireLength[$i]/2);
$edgeZ_p = $zpos + ($SinVAngle*$TopCornerVWireLength[$i]/2);
$edgeZ_m = $zpos - ($SinVAngle*$TopCornerVWireLength[$i]/2);
print $wout "V$j: ( $ypos , $zpos ) \n";
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
        <position name="pos${_[3]}WireZ$i" unit="cm" x="0" y="0 " z="@{[$zpos]}"/>
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

    my $BottomOfAPA = - $TPC_y/2 + $APAGap_y/2;


    $posZplane[0]   = -$_[0]/2 + $APAWirePlaneSpacing - $TPCWirePlaneThickness/2;
    $posZplane[1]   = $BottomOfAPA + $WrapCover + 4*$G10thickness + $Zactive_y/2;
    $posZplane[2]   = 0;

    $posVplane[0]   = $posZplane[0] + $APAWirePlaneSpacing;
    $posVplane[1]   = $BottomOfAPA + $WrapCover + 3*$G10thickness + $Vactive_y/2;
    $posVplane[2]   = $posZplane[2];

    $posUplane[0]   = $posVplane[0] + $APAWirePlaneSpacing;
    $posUplane[1]   = $BottomOfAPA + $WrapCover + 2*$G10thickness + $Uactive_y/2;
    $posUplane[2]   = $posZplane[2];

    $posTPCActive[0] = $posUplane[0] + ($TPCWirePlaneThickness/2 + $TPCActive_x/2);
    $posTPCActive[1] = -$_[1]/2 + $TPCActive_y/2;
    $posTPCActive[2] = 0;

#wrap up the TPC file
print TPC <<EOF;
    <volume name="volTPC${_[3]}">
      <materialref ref="LAr"/>
      <solidref ref="${_[3]}"/>
     <physvol>
       <volumeref ref="volTPCPlaneZ${_[3]}"/>
       <position name="pos${_[3]}PlaneZ" unit="cm" 
         x="@{[$posZplane[0]]}" y="@{[$posZplane[1]]}" z="@{[$posZplane[2]]}"/>
       <rotationref ref="rIdentity"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneV${_[3]}"/>
       <position name="pos${_[3]}PlaneV" unit="cm" 
         x="@{[$posVplane[0]]}" y="@{[$posVplane[1]]}" z="@{[$posVplane[2]]}"/>
       <rotationref ref="rIdentity"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneU${_[3]}"/>
       <position name="pos${_[3]}PlaneU" unit="cm" 
         x="@{[$posUplane[0]]}" y="@{[$posUplane[1]]}" z="@{[$posUplane[2]]}"/>
       <rotationref ref="rIdentity"/>
     </physvol>
EOF

# Option to omit the active volume on the outside for the prototype.
if( !($killOuterActive==1 && $_[3] eq 'Outer') ){
print TPC <<EOF;
     <physvol>
       <volumeref ref="volTPCActive${_[3]}"/>
       <position name="pos${_[3]}Active" unit="cm" 
         x="@{[$posTPCActive[0]]}" y="@{[$posTPCActive[1]]}" z="@{[$posTPCActive[2]]}"/>
       <rotationref ref="rIdentity"/>
     </physvol>
EOF
}

print TPC <<EOF;
    </volume>
</structure>
</gdml>
EOF

    close(TPC);

} #end of sub gen_TPC



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_Cryostat +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Cryostat()
{

# Create the cryostat fragment file name,
# add file to list of output GDML fragments,
# and open it
    $CRYO = "dune10kt_v3_Cryostat" . $suffix . ".gdml";
    push (@gdmlFiles, $CRYO);
    $CRYO = ">" . $CRYO;
    open(CRYO) or die("Could not open file $CRYO for writing");


# The standard XML prefix and starting the gdml
    print CRYO <<EOF;
<?xml version='1.0'?>
<gdml>
EOF

$APAFrameZSide_x = $APAFrame_x;
$APAFrameZSide_y = 4*$inch;
$APAFrameZSide_z = $APAFrame_z;

$APAFrameYSide_x = $APAFrame_x;
$APAFrameYSide_y = $APAFrame_y-2*$APAFrameZSide_y;
$APAFrameYSide_z = 4*$inch;

# Two outer Y supports will sandwich the light paddles
$APAFrameYOuterSupport_x = ($APAFrame_x-$LightPaddle_x)/2;
$APAFrameYOuterSupport_y = $APAFrame_y-2*$APAFrameZSide_y;
$APAFrameYOuterSupport_z = 4*$inch;

$EdgeFrameSteelThickness = 0.12*$inch;
$InnerFrameSteelThickness = 0.062*$inch;


$G10BoardYSide_x = $APAFrame_x;
$G10BoardYSide_y = $APAFrame_y;
$G10BoardYSide_z = $G10thickness;

$G10BoardZSide_x = $APAFrame_x;
$G10BoardZSide_y = $G10thickness;
$G10BoardZSide_z = $APAFrame_z;

##################################################################
################### Field Cage Y and Z centers ###################


$FC_x_center           = 0.0;
$FC_y_center           = 0.0;
$FC_z_center           = 0.0;
for($k=0 ; $k<$nAPALong ; $k++) #TODO!!! Factor out these loops so that the calculation of each APA center for placing APAs (happens later) will always agree with placing the center of the Field Cage.
{
	$tmp_APACenter_z  =  - $Argon_z/2 
	+ $UpstreamLArPadding 
	+ $APAphys_z/2 
	+ $k*($APAphys_z+$APAGap_z);
	$FC_z_center = $FC_z_center + $tmp_APACenter_z;
}
$FC_z_center = $FC_z_center / $nAPALong;

for($j=0 ; $j<$nAPAHigh ; $j++)  # nAPAHigh always going to be 2
{
	$tmp_APACenter_y  =  - $Argon_y/2 + $SpaceAPAToFloor
	+ $APAphys_y/2 
	+ $j*($APAphys_y+$APAGap_y);
	$FC_y_center = $FC_y_center + $tmp_APACenter_y;
}
$FC_y_center = $FC_y_center / $nAPAHigh;

# All the cryostat solids.
print CRYO <<EOF;
<solids>
    <tube name="SourceModerator"
      rmax="@{[0.5*$SourceModerator_d]}"
      z="@{[$SourceModerator_z]}"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
    <box name="Cryostat" lunit="cm" 
      x="@{[$Cryostat_x]}" 
      y="@{[$Cryostat_y]}" 
      z="@{[$Cryostat_z]}"/>
    <box name="ArgonInterior" lunit="cm" 
      x="@{[$Argon_x]}"
      y="@{[$Argon_y]}"
      z="@{[$Argon_z]}"/>
    <box name="FieldCageOut" lunit="cm"
      x="@{[$FC_x]}"
      y="@{[$FCO_y]}"
      z="@{[$FCO_z]}"/>
    <box name="FieldCageIn" lunit="cm"
      x="@{[$FC_x+1]}"
      y="@{[$FCI_y]}"
      z="@{[$FCI_z]}"/>
    <subtraction name="FieldCage">
      <first ref="FieldCageOut"/>
      <second ref="FieldCageIn"/>
    </subtraction>
    <box name="GaseousArgon" lunit="cm" 
      x="@{[$Argon_x]}"
      y="@{[$HeightGaseousAr]}"
      z="@{[$Argon_z]}"/>
    <subtraction name="SteelShell">
      <first ref="Cryostat"/>
      <second ref="ArgonInterior"/>
    </subtraction>
    <box name="Cathode" lunit="cm"
      x="@{[$Cathode_x]}"
      y="@{[$Cathode_y]}"
      z="@{[$Cathode_z]}"/>

    <box name="LightPaddle" lunit="cm"
      x="@{[$LightPaddle_x]}"
      y="@{[$LightPaddle_y]}"
      z="@{[$LightPaddle_z + $SiPM_z]}"/>

     <box name="APAFrameYSideHollow" lunit="cm"
      x="@{[$APAFrameYSide_x-2*$EdgeFrameSteelThickness]}"
      y="@{[$APAFrameYSide_y-2*$EdgeFrameSteelThickness]}"
      z="@{[$APAFrameYSide_z]}"/>
     <box name="APAFrameYSideShell" lunit="cm"
      x="@{[$APAFrameYSide_x]}"
      y="@{[$APAFrameYSide_y]}"
      z="@{[$APAFrameYSide_z]}"/>
     <subtraction name="APAFrameYSide">
      <first  ref="APAFrameYSideShell"/>
      <second ref="APAFrameYSideHollow"/>
      <positionref ref="posCenter"/>
      <rotationref ref="rIdentity"/>
      </subtraction>

     <box name="APAFrameZSideHollow" lunit="cm"
      x="@{[$APAFrameZSide_x-2*$EdgeFrameSteelThickness]}"
      y="@{[$APAFrameZSide_y-2*$EdgeFrameSteelThickness]}"
      z="@{[$APAFrameZSide_z]}"/>
     <box name="APAFrameZSideShell" lunit="cm"
      x="@{[$APAFrameZSide_x]}"
      y="@{[$APAFrameZSide_y]}"
      z="@{[$APAFrameZSide_z]}"/>
     <subtraction name="APAFrameZSide">
      <first  ref="APAFrameZSideShell"/>
      <second ref="APAFrameZSideHollow"/>
      <positionref ref="posCenter"/>
      <rotationref ref="rIdentity"/>
      </subtraction>

     <box name="APAFrameYOuterSupport" lunit="cm"
      x="@{[$EdgeFrameSteelThickness]}"
      y="@{[$APAFrameYOuterSupport_y]}"
      z="@{[$APAFrameYOuterSupport_z]}"/>


     <box name="G10BoardYSideCenterSeg" lunit="cm"
      x="@{[$G10BoardYSide_x]}"
      y="@{[$G10BoardYSide_y]}"
      z="@{[$G10BoardYSide_z]}"/>

     <box name="G10BoardZSideCenterSeg" lunit="cm"
      x="@{[$G10BoardZSide_x]}"
      y="@{[$G10BoardZSide_y]}"
      z="@{[$G10BoardZSide_z]}"/>

</solids>
EOF

# Cryostat structure
print CRYO <<EOF;
<structure>
    <volume name="volSteelShell">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="SteelShell" />
    </volume>
    <volume name="volSourceModerator">
      <materialref ref="Delrin"/>
      <solidref ref="SourceModerator"/>
    </volume>
    <volume name="volFieldCage">
      <materialref ref="ALUMINUM_Al"/>
      <solidref ref="FieldCage"/>
    </volume>
    <volume name="volGaseousArgon">
      <materialref ref="ArGas"/>
      <solidref ref="GaseousArgon"/>
    </volume>
    <volume name="volCathode">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="Cathode" />
    </volume>
EOF

for($i=0 ; $i<$nAPAs ; $i++){
for($p=0 ; $p<10 ; $p++){
print CRYO <<EOF;
   <volume name="volOpDetSensitive_$i\-$p">
     <materialref ref="Acrylic"/>
     <solidref ref="LightPaddle"/>
   </volume>
EOF
}
}


print CRYO <<EOF;

    <volume name="volAPAFrameYSide">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="APAFrameYSide"/>
    </volume>

    <volume name="volAPAFrameZSide">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="APAFrameZSide"/>
    </volume>

    <volume name="volAPAFrameYOuterSupport">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="APAFrameYOuterSupport"/>
    </volume>

    <volume name="volG10BoardYSideCenterSeg">
      <materialref ref="G10"/>
      <solidref ref="G10BoardYSideCenterSeg"/>
    </volume>

    <volume name="volG10BoardZSideCenterSeg">
      <materialref ref="G10"/>
      <solidref ref="G10BoardZSideCenterSeg"/>
    </volume>

    <volume name="volCryostat">
      <materialref ref="LAr" />
      <solidref ref="Cryostat" />
      <physvol>
        <volumeref ref="volGaseousArgon"/>
        <position name="posGaseousArgon" unit="cm" x="0" y="@{[$Argon_y/2-$HeightGaseousAr/2]}" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volSourceModerator"/>
        <position name="posSourceModerator" unit="cm"
          x="@{[$Source_x_center]}"
	        y="@{[$Source_y_center]}"
	        z="@{[$Source_z_center]}"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>    
      <physvol>
        <volumeref ref="volFieldCage"/>
        <position name="posFieldCage" unit="cm" 
          x="@{[$FC_x_center]}"
	        y="@{[$FC_y_center]}"
	        z="@{[$FC_z_center]}"/>
      </physvol>
      <physvol>
        <volumeref ref="volSteelShell"/>
        <position name="posSteelShell" unit="cm" x="0" y="0" z="0"/>
      </physvol>
EOF

# nested for loops to place the non-rotated AND rotated volTPC
   # x loop rotation: There are six drift volumes. Looking into the 
   # detector from incident direction, and counting from right (-x) to
   # left (+x), odd number volumes need to be rotated 180deg about Y in
   # order for the cathode to be on the right of the APA.

if ($tpc_on==1) {

    $apa_i=0;
    $cpa_i=0;

  for($k=0 ; $k<$nAPALong ; $k++)
  {
        $APACenter_z  =  - $Argon_z/2 
                       + $UpstreamLArPadding 
                   + $APAphys_z/2 
                   + $k*($APAphys_z+$APAGap_z);


    for($j=0 ; $j<$nAPAHigh ; $j++)  # nAPAHigh always going to be 2
    {        # $j=0 for bottom APAs
        # $j=1 for top APAs
        $APACenter_y  =  - $Argon_y/2 + $SpaceAPAToFloor
                       + $APAphys_y/2 
                       + $j*($APAphys_y+$APAGap_y);


      if( $workspace == 0 )
      {

          for($i=0 ; $i<$nAPAWide ; $i++)
          {


            $APACenter_x  =  - $Argon_x/2
                             + $SpaceAPAToCryoWall + $APA_UtoU_x/2 + $TPCWirePlaneThickness/2
                                 + $i*$APAToAPA;

            $CPA_x        =  $APACenter_x  +  $CPAToAPA;


            place_APA($APACenter_x, $APACenter_y, $APACenter_z, $apa_i, $j);
            place_OpDets($APACenter_x, $APACenter_y, $APACenter_z, $apa_i);
            $tpc_0 = 2*$apa_i+0;
            $tpc_1 = 2*$apa_i+1;
            $apa_i++;
            
            $SelectTPC_0 = "Inner";
            $SelectTPC_1 = "Inner";
            $TPC_0_x     = $TPCInner_x;
            $TPC_1_x     = $TPCInner_x;
            $rot_0       = "rPlus180AboutY";
            $rot_1       = "rIdentity";
            if($i == 0)           { $SelectTPC_0 = "Outer"; $TPC_0_x = $TPCOuter_x; }
            if($i == $nAPAWide-1) { $SelectTPC_1 = "Outer"; $TPC_1_x = $TPCOuter_x; }
            if($j == 0)           { $rot_0 = "rPlus180AboutXPlus180AboutY"; 
                              $rot_1 = "rPlus180AboutX"; } #put the readout end at the bottom for bottom APAs



print CRYO <<EOF;

      <physvol>
        <volumeref ref="volTPC$SelectTPC_0"/>
        <position name="posTPC\-$tpc_0" unit="cm" 
        x="@{[$APACenter_x - $APAFrame_x/2 - $TPC_0_x/2]}" 
        y="@{[$APACenter_y]}" 
        z="@{[$APACenter_z]}"/> 
        <rotationref ref="$rot_0"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPC$SelectTPC_1"/>
        <position name="posTPC\-$tpc_1" unit="cm" 
        x="@{[$APACenter_x + $APAFrame_x/2 + $TPC_1_x/2]}" 
        y="@{[$APACenter_y]}" 
        z="@{[$APACenter_z]}"/> 
        <rotationref ref="$rot_1"/>
      </physvol>

EOF


if( $i < $nAPAWide - 1 ){ # avoid placeing the last row of CPAs since the APAs are on the outside

print CRYO <<EOF;

      <physvol>
        <volumeref ref="volCathode"/>
        <position name="posCathode\-$cpa_i" unit="cm" 
      x="@{[$CPA_x]}" 
      y="@{[$APACenter_y]}" 
      z="@{[$APACenter_z]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>

EOF
          $cpa_i++;
}


          } #wide

      } # end if not the smaller workspace


        # Make the workspace have only one center APA with CPAs and the full drift on either side
        elsif( $workspace == 1 || $workspace==2 ){

            $APACenter_x  =  0;
            $CPA_0_x      =  $APACenter_x  -  $CPAToAPA;
            $CPA_1_x      =  $APACenter_x  +  $CPAToAPA;


            place_APA($APACenter_x, $APACenter_y, $APACenter_z, $apa_i, $j);
            place_OpDets($APACenter_x, $APACenter_y, $APACenter_z, $apa_i);

            $tpc_0 = 2*$apa_i+0;
            $tpc_1 = 2*$apa_i+1;
            $apa_i++;

            $rot_0       = "rPlus180AboutY";
            $rot_1       = "rIdentity";
            if($j == 0)           { $rot_0 = "rPlus180AboutXPlus180AboutY"; 
                              $rot_1 = "rPlus180AboutX"; } #put the readout end at the bottom for bottom APAs
            
print CRYO <<EOF;

      <physvol>
        <volumeref ref="volTPCInner"/>
        <position name="posTPC\-$tpc_0" unit="cm" 
        x="@{[$APACenter_x - $APAFrame_x/2 - $TPCInner_x/2]}" 
        y="@{[$APACenter_y]}" 
        z="@{[$APACenter_z]}"/> 
        <rotationref ref="$rot_0"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPCInner"/>
        <position name="posTPC\-$tpc_1" unit="cm"
        x="@{[$APACenter_x + $APAFrame_x/2 + $TPCInner_x/2]}" 
        y="@{[$APACenter_y]}"
        z="@{[$APACenter_z]}"/>
        <rotationref ref="$rot_1"/>
      </physvol>

      <physvol>
        <volumeref ref="volCathode"/>
        <position name="posCathode\-$apa_i-0" unit="cm" 
      x="@{[$CPA_0_x]}" 
      y="@{[$APACenter_y]}" 
      z="@{[$APACenter_z]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volCathode"/>
        <position name="posCathode\-$apa_i-1" unit="cm" 
      x="@{[$CPA_1_x]}" 
      y="@{[$APACenter_y]}" 
      z="@{[$APACenter_z]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>

EOF

        } # if the smaller workspace


    } #high
  } #long

}# end if tpc

print CRYO <<EOF;
    </volume>
</structure>
</gdml>
EOF

close(CRYO);
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++ place_OpDets +++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub place_OpDets()
{

    $APACenter_x = $_[0];
    $APACenter_y = $_[1];
    $APACenter_z = $_[2];
    $apa_i = $_[3];

for ($paddle = 0; $paddle<$nLightPaddlesPerAPA; $paddle++)
{

             # All Light Paddle centers will have the same
             #      X coordinate as the center of the current APA
             #      Z coordinate as the current TPC pair
             # The Y coordinate must be looped over:

             #the multiplication by j here is a temporary dirty trick to get around some strange behavior

             $Paddle_Y   =    $APACenter_y 
                            - $APAphys_y/2 
                            + $j*$FrameToPaddleSpace 
                            + (1-$j)*($LightPaddle_y/2 + $APAFrameZSide_y) 
                            + $PaddleYInterval*$paddle; 

             # Alternate the paddle orientations
             if ( $paddle % 2 == 0 ) { $rot = "rIdentity"; }
             else                    { $rot = "rPlus180AboutY"; }

             
             print CRYO <<EOF;
     <physvol>
       <volumeref ref="volOpDetSensitive_$apa_i\-$paddle"/>
       <position name="posPaddle\-$paddle\-TPC\-$i\-$j\-$k" unit="cm" 
         x="@{[$APACenter_x]}" 
       y="@{[$Paddle_Y]}" 
         z="@{[$APACenter_z + $SiPM_z/2]}"/>
       <rotationref ref="$rot"/>
     </physvol>
EOF

}#end Paddle for-loop

}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ place_APA ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Must be called only within gen_Cryostat(), 


# $_[0] = x APA physical center
# $_[1] = y APA physical center
# $_[2] = z APA physical center
# $_[3] = APA number
# $_[4] = 0 for APA with readout on bottom
#         1 for APA with readout on top

sub place_APA()
{

    $APA_i = $_[3];

####################################################################
################# APA Frame and Paddle Dimensions ##################

# The center passed to this function is the physical APA center,
# which is not quite the frame's center, since there are more boards
# at the bottom. Transform them:

    $APAFrameCenter_x =   $_[0];
    $APAFrameCenter_z =   $_[2];

# The bottom APAs are upside down relative to how the TPC is initially generated,
#   and the center needs to be slightly shifted since there is extra G10 on the
#   end that doesn't read out.

    if($_[4]==1)    # top APAs
    {
      $APAFrameCenter_y =   $_[1] - $APAphys_y/2 
                          + $WrapCover + 4*$G10thickness
                      + $APAFrame_y/2;
    }
    elsif($_[4]==0) # bottom APAs
    {
      $APAFrameCenter_y =   $_[1] + $APAphys_y/2 
                            - $WrapCover - 4*$G10thickness
                           - $APAFrame_y/2;
    }
    else{ print "APA not labeled as top or bottom"; }

$APAFrameZSide_x = $APAFrame_x;
$APAFrameZSide_y = 4*$inch;
$APAFrameZSide_z = $APAFrame_z;

$APAFrameYSide_x = $APAFrame_x;
$APAFrameYSide_y = $APAFrame_y-2*$APAFrameZSide_y;
$APAFrameYSide_z = 4*$inch;

# Two outer Y supports will sandwich the light paddles
$APAFrameYOuterSupport_x = ($APAFrame_x-$LightPaddle_x)/2;
$APAFrameYOuterSupport_y = $APAFrame_y-2*$APAFrameZSide_y;
$APAFrameYOuterSupport_z = 4*$inch;

# if there were an inner support to fill the hole
$APAFrameYInnerSupport_x = $LightPaddle_x;

$EdgeFrameSteelThickness = 0.12*$inch;
$InnerFrameSteelThickness = 0.062*$inch;


$G10BoardYSide_x = $APAFrame_x;
$G10BoardYSide_y = $APAFrame_y;
$G10BoardYSide_z = $G10thickness;

$G10BoardZSide_x = $APAFrame_x;
$G10BoardZSide_y = $G10thickness;
$G10BoardZSide_z = $APAFrame_z;

if($_[4]==1)    # top APAs
{
    $posG10ZSideZ_y    = $APAFrameCenter_y - $APAFrame_y/2 - (0+.5)*($G10BoardZSide_y);
    $posG10ZSideV_y    = $APAFrameCenter_y - $APAFrame_y/2 - (1+.5)*($G10BoardZSide_y);
    $posG10ZSideU_y    = $APAFrameCenter_y - $APAFrame_y/2 - (2+.5)*($G10BoardZSide_y);
    $posG10ZSideGrid_y = $APAFrameCenter_y - $APAFrame_y/2 - (3+.5)*($G10BoardZSide_y);
}
elsif($_[4]==0) # bottom APAs
{
    $posG10ZSideZ_y    = $APAFrameCenter_y + $APAFrame_y/2 + (0+.5)*($G10BoardZSide_y);
    $posG10ZSideV_y    = $APAFrameCenter_y + $APAFrame_y/2 + (1+.5)*($G10BoardZSide_y);
    $posG10ZSideU_y    = $APAFrameCenter_y + $APAFrame_y/2 + (2+.5)*($G10BoardZSide_y);
    $posG10ZSideGrid_y = $APAFrameCenter_y + $APAFrame_y/2 + (3+.5)*($G10BoardZSide_y);
}
else{ print "APA not labeled as top or bottom"; }

   # First put in the frame
#  print CRYO <<EOF;
#<!-- 
#      <physvol>
#        <volumeref ref="volAPAFrameYOuterSupport\-$APA_i"/>
#        <position name="posAPAFrameYOuterSupportNeg\-$APA_i" unit="cm" 
#      x="$APAFrameCenter_x - ($APAFrameYOuterSupport_x + $APAFrameYInnerSupport_x/2 - $EdgeFrameSteelThickness/2)" 
#      y="$APAFrameCenter_y" 
#      z="$APAFrameCenter_z"/> 
#        <rotationref ref="rIdentity"/>
#      </physvol>
#      <physvol>
#        <volumeref ref="volAPAFrameYOuterSupport\-$APA_i"/>
#        <position name="posAPAFrameYOuterSupportPos\-$APA_i" unit="cm" 
#      x="$APAFrameCenter_x + ($APAFrameYOuterSupport_x + $APAFrameYInnerSupport_x/2 - $EdgeFrameSteelThickness/2)" 
#      y="$APAFrameCenter_y" 
#      z="$APAFrameCenter_z"/> 
#        <rotationref ref="rIdentity"/>
#      </physvol>
#-->
#EOF
  print CRYO <<EOF;
      <physvol>
        <volumeref ref="volAPAFrameYSide"/>
        <position name="posAPAFrameYSideNeg\-$APA_i" unit="cm" 
      x="@{[$APAFrameCenter_x]}" 
      y="@{[$APAFrameCenter_y]}" 
      z="@{[$APAFrameCenter_z - $APAFrame_z/2 + $APAFrameYSide_z/2]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameYSide"/>
        <position name="posAPAFrameYSidePos\-$APA_i" unit="cm" 
      x="@{[$APAFrameCenter_x]}" 
      y="@{[$APAFrameCenter_y]}" 
      z="@{[$APAFrameCenter_z + $APAFrame_z/2 - $APAFrameYSide_z/2]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameZSide"/>
        <position name="posAPAFrameZSidePos\-$APA_i" unit="cm" 
      x="@{[$APAFrameCenter_x]}" 
      y="@{[$APAFrameCenter_y + $APAFrame_y/2 - $APAFrameZSide_y/2]}" 
      z="@{[$APAFrameCenter_z]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameZSide"/>
        <position name="posAPAFrameZSideNeg\-$APA_i" unit="cm" 
      x="@{[$APAFrameCenter_x]}" 
      y="@{[$APAFrameCenter_y  - $APAFrame_y/2 + $APAFrameZSide_y/2]}" 
      z="@{[$APAFrameCenter_z]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>


      <!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        - Add the *parts* of the G10 boards that exist directly in volCryostat.
      - There are two boards on each the up and downstream end, 
          one each to wrap the U and V views around the APA frame
      - There are 4 on the bottom which anchor the U V and Z wires and the grid plane
      - The rest of the parts of the G10 boards must be placed directly in volTPC   -->

      <physvol>
        <volumeref ref="volG10BoardYSideCenterSeg"/>
        <position name="posG10BoardYSideCenterSeg\-Vup\-$APA_i" unit="cm" 
      x="@{[$APAFrameCenter_x]}" 
      y="@{[$APAFrameCenter_y]}" 
      z="@{[$APAFrameCenter_z - $APAFrame_z/2 - (0+.5)*($G10BoardYSide_z)]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardYSideCenterSeg"/>
        <position name="posG10BoardYSideCenterSeg\-Uup\-$APA_i" unit="cm" 
      x="@{[$APAFrameCenter_x]}" 
      y="@{[$APAFrameCenter_y]}" 
      z="@{[$APAFrameCenter_z - $APAFrame_z/2 - (1+.5)*($G10BoardYSide_z)]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>

      <physvol>
        <volumeref ref="volG10BoardYSideCenterSeg"/>
        <position name="posG10BoardYSideCenterSeg\-Vdown\-$APA_i" unit="cm" 
      x="@{[$APAFrameCenter_x]}" 
      y="@{[$APAFrameCenter_y]}" 
      z="@{[$APAFrameCenter_z + $APAFrame_z/2 + (0+.5)*($G10BoardYSide_z)]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardYSideCenterSeg"/>
        <position name="posG10BoardYSideCenterSeg\-Udown\-$APA_i" unit="cm" 
      x="@{[$APAFrameCenter_x]}" 
      y="@{[$APAFrameCenter_y]}" 
      z="@{[$APAFrameCenter_z + $APAFrame_z/2 + (1+.5)*($G10BoardYSide_z)]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>

      <physvol>
        <volumeref ref="volG10BoardZSideCenterSeg"/>
        <position name="posG10BoardZSideCenterSeg\-Z\-$APA_i" unit="cm" 
      x="@{[$APAFrameCenter_x]}" 
      y="@{[$posG10ZSideZ_y]}" 
      z="@{[$APAFrameCenter_z]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardZSideCenterSeg"/>
        <position name="posG10BoardZSideCenterSeg\-V\-$APA_i" unit="cm" 
      x="@{[$APAFrameCenter_x]}" 
      y="@{[$posG10ZSideV_y]}" 
      z="@{[$APAFrameCenter_z]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardZSideCenterSeg"/>
        <position name="posG10BoardZSideCenterSeg\-U\-$APA_i" unit="cm" 
      x="@{[$APAFrameCenter_x]}" 
      y="@{[$posG10ZSideU_y]}" 
      z="@{[$APAFrameCenter_z]}"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardZSideCenterSeg"/>
        <position name="posG10BoardZSideCenterSeg\-Grid\-$APA_i" unit="cm" 
      x="@{[$APAFrameCenter_x]}" 
      y="@{[$posG10ZSideGrid_y]}"
      z="@{[$APAFrameCenter_z]}"/> 
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
    $ENCL = "dune10kt_v3_DetEnclosure" . $suffix . ".gdml";
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

    <box name="FoamPadBlock" lunit="cm"
      x="@{[$Cryostat_x + 2*$FoamPadding]}"
      y="@{[$Cryostat_y + 2*$FoamPadding]}"
      z="@{[$Cryostat_z + 2*$FoamPadding]}" />

    <box name="FoamPadInner" lunit="cm"
      x="@{[$Cryostat_x + 0.01]}"
      y="@{[$Cryostat_y + 0.01]}"
      z="@{[$Cryostat_z + 0.01]}" />

    <subtraction name="FoamPadding">
      <first ref="FoamPadBlock"/>
      <second ref="FoamPadInner"/>
      <positionref ref="posCenter"/>
    </subtraction>

    <box name="SteelSupportBlock" lunit="cm"
      x="@{[$Cryostat_x + 2*$FoamPadding + 2*$SteelSupport_x]}"
      y="@{[$Cryostat_y + 2*$FoamPadding + 2*$SteelSupport_y]}"
      z="@{[$Cryostat_z + 2*$FoamPadding + 2*$SteelSupport_z]}" />

    <box name="SteelSupportInner" lunit="cm"
      x="@{[$Cryostat_x + 2*$FoamPadding + 0.01]}"
      y="@{[$Cryostat_y + 2*$FoamPadding + 0.01]}"
      z="@{[$Cryostat_z + 2*$FoamPadding + 0.01]}" />

    <subtraction name="SteelSupport">
      <first ref="SteelSupportBlock"/>
      <second ref="SteelSupportInner"/>
      <positionref ref="posCenter"/>
    </subtraction>

    <box name="DetEnclosure" lunit="cm" 
      x="@{[$DetEncWidth]}"
      y="@{[$DetEncHeight]}"
      z="@{[$DetEncLength]}"/>

</solids>
EOF



# Detector enclosure structure
    print ENCL <<EOF;
<structure>
    <volume name="volFoamPadding">
      <materialref ref="FD_foam"/>
      <solidref ref="FoamPadding"/>
    </volume>

    <volume name="volSteelSupport">
      <materialref ref="AirSteelMixture"/>
      <solidref ref="SteelSupport"/>
    </volume>

    <volume name="volDetEnclosure">
      <materialref ref="Air"/>
      <solidref ref="DetEnclosure"/>

       <physvol>
           <volumeref ref="volFoamPadding"/>
           <positionref ref="posCryoInDetEnc"/>
       </physvol>
       <physvol>
           <volumeref ref="volSteelSupport"/>
           <positionref ref="posCryoInDetEnc"/>
       </physvol>
       <physvol>
           <volumeref ref="volCryostat"/>
           <positionref ref="posCryoInDetEnc"/>
       </physvol>
EOF



    print ENCL <<EOF;
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
    $WORLD = "dune10kt_v3_World" . $suffix . ".gdml";
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
      x="@{[$DetEncWidth+2*$RockThickness]}" 
      y="@{[$DetEncHeight+2*$RockThickness]}" 
      z="@{[$DetEncLength+2*$RockThickness]}"/>
</solids>
EOF

# World structure
print WORLD <<EOF;
<structure>
    <volume name="volWorld" >
      <materialref ref="DUSEL_Rock"/>
      <solidref ref="World"/>

      <physvol>
        <volumeref ref="volDetEnclosure"/>
      <position name="posDetEnclosure" unit="cm" x="@{[$OriginXSet]}" y="@{[$OriginYSet]}" z="@{[$OriginZSet]}"/>
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

   # This code is taken straigh from the similar MicroBooNE generate script, Thank you.

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



# run the sub routines that generate the fragments

gen_Define();        # generates definitions at beginning of GDML
gen_Materials(); # generates materials to be used


    gen_TPC( $TPCInner_x,  $TPC_y,  $TPC_z,  'Inner');
    gen_TPC( $TPCOuter_x,  $TPC_y,  $TPC_z,  'Outer');

close $wout;

gen_Cryostat();       # places (2*nAPAWide x nAPAHigh x nAPALong) volTPC,
             # half rotated 180 about Y
gen_Enclosure(); # places two cryostats and concrete volumes

gen_World();       # places the enclosure among DUSEL Rock


write_fragments(); # writes the XML input for make_gdml.pl
                  # which zips together the final GDML
exit;
