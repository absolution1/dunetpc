#!/usr/bin/perl

##################################################################################
#
#  This script is based on generate_dunedphase10kt_v2.pl (see description 
#  in there file for more info). Some parameters were changed to generate 
#  a first version of the protodune dual phase geometry (e.g. driftTPCActive
#  and number of CRM's). The beam window is not included in this version! 
#  For more info, please contact Christoph Alt: christoph.alt@cern.ch 
#
##################################################################################


#use warnings;
use gdmlMaterials;
use Math::Trig;
use Getopt::Long;
use Math::BigFloat;
Math::BigFloat->precision(-16);

GetOptions( "help|h" => \$help,
	    "suffix|s:s" => \$suffix,
	    "output|o:s" => \$output,
	    "wires|w:s" => \$wires,
	    "optical|opt:s" => \$optical);

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

#if ( ! defined $workspace ) # not done 
#{
#    $workspace = 0;
#}
#elsif ( $workspace == 1)
#{
#    print "\t\tCreating smaller workspace geometry.\n";
#}
$basename = "protodunedphase";


my $pmt_switch="on";
my $FieldCage_switch="on";
$GroundGrid_switch="off";
$Cathode_switch="off";
$ExtractionGrid_switch="off";
$LEMs_switch="off";

if (defined $optical)
{
    if ($optical==1)
    {
      $GroundGrid_switch="on";
      $Cathode_switch="on";
      $ExtractionGrid_switch="on";
      $basename = $basename."_optical";
      $LEMs_switch="on";	
    }
}


# set wires on to be the default, unless given an input by the user

$wires_on = 1; # 1=on, 0=off
if (defined $wires)
{
    $wires_on = $wires
}

$tpc_on = 1;


if ( $wires_on == 0 )
{
    $basename = $basename."_nowires";
}

#if ( $workspace == 1 )
#{
#    $basename = $basename."_workspace";
#}


##################################################################
############## Parameters for Charge Readout Plane ###############

# dune10kt dual-phase
$wirePitch           = 0.3125;  # channel pitch
$nChannelsViewPerCRM = 960;     # channels per collection view
$borderCRM           = 0.5;     # dead space at the border of each CRM

# dimensions of a single Charge Readout Module (CRM)
$widthCRM_active  = $wirePitch * $nChannelsViewPerCRM;
$lengthCRM_active = $wirePitch * $nChannelsViewPerCRM;

$widthCRM  = $widthCRM_active + 2 * $borderCRM;
$lengthCRM = $lengthCRM_active + 2 * $borderCRM;

# number of CRMs in y and z
$nCRM_y   = 2;
$nCRM_z   = 2;

# calculate tpc area based on number of CRMs and their dimensions
$widthTPCActive  = $nCRM_y * $widthCRM;  # around 600
$lengthTPCActive = $nCRM_z * $lengthCRM; # around 600

# active volume dimensions 
$driftTPCActive  = 597.0;

# model anode strips as wires
$padWidth  = 0.015;
$ReadoutPlane = 2 * $padWidth;

#$padHeight = 0.0035; 

##################################################################
############## Parameters for TPC and inner volume ###############

# inner volume dimensions of the cryostat
$Argon_x = 790.0;
$Argon_y = 854.8;
$Argon_z = 854.8;


# width of gas argon layer on top
#$HeightGaseousAr = 81.8;
$HeightGaseousAr = 51.5;

# size of liquid argon buffer
$xLArBuffer = $Argon_x - $driftTPCActive - $HeightGaseousAr - $ReadoutPlane;
$yLArBuffer = 0.5 * ($Argon_y - $widthTPCActive);
$zLArBuffer = 0.5 * ($Argon_z - $lengthTPCActive);

# cryostat 
#$SteelThickness = 0.12; # membrane
$SteelThickness = 0.2; # membrane

$Cryostat_x = $Argon_x + 2*$SteelThickness;
$Cryostat_y = $Argon_y + 2*$SteelThickness;
$Cryostat_z = $Argon_z + 2*$SteelThickness;

##################################################################
############## Field Cage Parameters ###############

$FieldShaperSizeX = 304.9;
$FieldShaperSizeY = 5.999;
$FieldShaperSizeZ = 4.603;

$FieldCageShaperProfileShift=-4.82;
$NFieldShapers = 100.0;
$FFSSeparation = 6.0;

$FieldCageSizeX = $FieldShaperSizeX+0.2;
$FieldCageSizeY = $FieldShaperSizeY+0.2;
$FieldCageSizeZ = $NFieldShapers*$FFSSeparation;

$FFSPositionX = 277.5;
$FFSPositionY = 153.521;
$FFSPositionZ = 309.593;


##################################################################
############## Parameters for PMTs ###############

#pos in cm inside the cryostat
 #$pmt_pos_x  =  -300-(5.5 * 2.54)+5.115; #just below the active volume
 #$pmt_pos_x  =  -300-(5.5 * 2.54)+5.115 + 22 -141.5; just on the top of the cryostat (new version)
 ##$pmt_pos_x  =  -300 + 5.115 - 1.27*2.54 + 22 -107.8;  #8.4cm from the bottom volTPCActive
 $HeightPMT = 37.0;

 #   $pmt_pos_x =  -$Argon_x/2 - 0.5*($HeightPMT);

#pmts not equally spaced:
 @pmt_pos = ( ' y="-238" z="-170" ', #pmt1
	' y="-238" z="-34" ', #pmt2
	' y="-238" z="34" ', #pmt3
	' y="-238" z="170" ', #pmt4
	' y="-170" z="-238" ', #pmt5
	' y="-170" z="-102" ', #pmt6
	' y="-170" z="102" ', #pmt7
	' y="-170" z="238" ', #pmt8
	' y="-102" z="-170" ', #pmt9
	' y="-102" z="-34" ', #pmt10
	' y="-102" z="34" ', #pmt11
	' y="-102" z="170" ', #pmt12
	' y="-34" z="-238" ', #pmt13
	' y="-34" z="-102" ', #pmt14
	' y="-34" z="-34" ', #pmt15
	' y="-34" z="34" ', #pmt16
	' y="-34" z="102" ', #pmt17
	' y="-34" z="238" ', #pmt18
	' y="34" z="-238" ', #pmt19
	' y="34" z="-102" ', #pmt20
	' y="34" z="-34" ', #pmt21
	' y="34" z="34" ', #pmt22
	' y="34" z="102" ', #pmt23
	' y="34" z="238" ', #pmt24
	' y="102" z="-170" ', #pmt25
	' y="102" z="-34" ', #pmt26
	' y="102" z="34" ', #pmt27
	' y="102" z="170" ', #pmt28
	' y="170" z="-238" ', #pmt29
	' y="170" z="-102" ', #pmt30
	' y="170" z="102" ', #pmt31
	' y="170" z="238" ', #pmt32
	' y="238" z="-170" ', #pmt33
	' y="238" z="-34" ', #pmt34
	' y="238" z="34" ', #pmt35
	' y="238" z="170" '); #PMT36

#pmts equally spaced:
# @pmt_pos = ( ' y="-162.5" z="-162.5"', #PMT1
#             ' y="-162.5" z="-97.5"', #PMT2
#             ' y="-162.5" z="-32.5"', #PMT3
#             ' y="-162.5" z="32.5"', #PMT4
#              ' y="-162.5" z="97.5"', #PMT5
#              ' y="-162.5" z="162.5"', #PMT6
#              ' y="-97.5" z="-162.5"', #PMT7
#              ' y="-97.5" z="-97.5"', #PMT8
#              ' y="-97.5" z="-32.5"', #PMT9
#              ' y="-97.5" z="32.5"', #PMT10#
#	      ' y="-97.5" z="97.5"', #PMT11
#              ' y="-97.5" z="162.5"', #PMT12
#              ' y="-32.5" z="-162.5"', #PMT13
#              ' y="-32.5" z="-97.5"', #PMT14
#              ' y="-32.5" z="-32.5"', #PMT15
#              ' y="-32.5" z="32.5"', #PMT16
#              ' y="-32.5" z="97.5"', #PMT17
#              ' y="-32.5" z="162.5"', #PMT18
#              ' y="32.5" z="-162.5"', #PMT19
#              ' y="32.5" z="-97.5"', #PMT20
#              ' y="32.5" z="-32.5"', #PMT21
#              ' y="32.5" z="32.5"', #PMT22
#              ' y="32.5" z="97.5"', #PMT23
#              ' y="32.5" z="162.5"', #PMT24
#              ' y="97.5" z="-162.5"', #PMT25
#              ' y="97.5" z="-97.5"', #PMT26
#              ' y="97.5" z="-32.5"', #PMT27
#              ' y="97.5" z="32.5"', #PMT27
#              ' y="97.5" z="97.5"', #PMT29
#              ' y="97.5" z="162.5"', #PMT30
#              ' y="162.5" z="-162.5"', #PMT31
#              ' y="162.5" z="-97.5"', #PMT32
#              ' y="162.5" z="-32.5"', #PMT33
#              ' y="162.5" z="32.5"', #PMT34
#              ' y="162.5" z="97.5"', #PMT35
#              ' y="162.5" z="162.5"'); #PMT36

##################################################################
############## DetEnc and World relevant parameters  #############

$SteelSupport_x  =  1;
$SteelSupport_y  =  1;
$SteelSupport_z  =  1;
$FoamPadding     =  80;

$FracMassOfSteel =  0.5; #The steel support is not a solid block, but a mixture of air and steel
$FracMassOfAir   =  1 - $FracMassOfSteel;


$SpaceSteelSupportToWall    = 200;
$SpaceSteelSupportToCeiling = 200;

$DetEncWidth   =    $Cryostat_x
                  + 2*($SteelSupport_x + $FoamPadding) + 2*$SpaceSteelSupportToWall;
$DetEncHeight  =    $Cryostat_y
                  + 2*($SteelSupport_y + $FoamPadding) + $SpaceSteelSupportToCeiling;
$DetEncLength  =    $Cryostat_z
                  + 2*($SteelSupport_z + $FoamPadding) + 2*$SpaceSteelSupportToWall;

$posCryoInDetEnc_y = 0;
$posCryoInDetEnc_y2 = - $DetEncHeight/2 + $SteelSupport_y + $FoamPadding + $Cryostat_y/2;

$posTopSteelStruct = $Argon_x/2+$FoamPadding+$SteelSupport_x;
$posBotSteelStruct = -$Argon_x/2-$FoamPadding-$SteelSupport_x;
$posZBackSteelStruct = $Argon_z/2+$FoamPadding+$SteelSupport_z;
$posZFrontSteelStruct = -$Argon_z/2-$FoamPadding-$SteelSupport_z;
$posLeftSteelStruct = -$Argon_y/2-$FoamPadding-$SteelSupport_y;
$posRightSteelStruct = $Argon_y/2+$FoamPadding+$SteelSupport_y;

$RockThickness = 3000;

  # We want the world origin to be at the very front of the fiducial volume.
  # move it to the front of the enclosure, then back it up through the concrete/foam, 
  # then through the Cryostat shell, then through the upstream dead LAr (including the
  # dead LAr on the edge of the TPC)
  # This is to be added to the z position of every volume in volWorld

$OriginZSet =   $DetEncLength/2.0 
              - $SpaceSteelSupportToWall
              - $SteelSupport_z
              - $FoamPadding
              - $SteelThickness
              - $zLArBuffer
              - $borderCRM;

  # We want the world origin to be vertically centered on active TPC
  # This is to be added to the y position of every volume in volWorld

$OriginYSet =   $DetEncHeight/2.0
              - $SteelSupport_y
              - $FoamPadding
              - $SteelThickness
              - $yLArBuffer
              - $widthTPCActive/2.0
              + $posCryoInDetEnc_y2;

$OriginXSet =  $DetEncWidth/2.0
              -$SpaceSteelSupportToWall
              -$SteelSupport_x
              -$FoamPadding
              -$SteelThickenss
              -$xLArBuffer
              -$driftTPCActive/2.0;


##################################################################
############### Parameters for det elements ######################

# cathode plane
#$Cathode_x = $widthTPCActive;
#$Cathode_y = 1.0;
#$Cathode_z = $lengthTPCActive;



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
    $DEF = $basename."_Def" . $suffix . ".gdml";
    push (@gdmlFiles, $DEF);
    $DEF = ">" . $DEF;
    open(DEF) or die("Could not open file $DEF for writing");


print DEF <<EOF;
<?xml version='1.0'?>
<gdml>
<define>

<!--



-->

   <position name="posCryoInDetEnc"     unit="cm" x="0" y="$posCryoInDetEnc_y" z="0"/>
   <position name="posCenter"           unit="cm" x="0" y="0" z="0"/>
   <rotation name="rPlus90AboutX"       unit="deg" x="90" y="0" z="0"/>
   <rotation name="rMinus90AboutX"       unit="deg" x="270" y="0" z="0"/>

   <rotation name="rMinus90AboutYMinus90AboutX"       unit="deg" x="270" y="270" z="0"/>
   <rotation name="rPlus180AboutXPlus180AboutY"	unit="deg" x="180" y="180"   z="0"/>
   <rotation name="rIdentity"		unit="deg" x="0" y="0"   z="0"/>

   <rotation name="rPlus180AboutX"	unit="deg" x="180" y="0"   z="0"/>
   <rotation name="rPlus180AboutY"	unit="deg" x="0" y="180"   z="0"/>
   <rotation name="rPlus180AboutZ"	unit="deg" x="0" y="0"   z="180"/>

   <rotation name="rPlus90AboutY"       unit="deg" x="0" y="90" z="0"/>
   <rotation name="rMinus90AboutY"      unit="deg" x="0" y="270" z="0"/>


   <rotation name="rot04"      unit="deg" x="0" y="270" z="90"/>
   <rotation name="rot07"      unit="deg" x="0" y="90" z="90"/>

   <rotation name="rot03"      unit="deg" x="0" y="90" z="270"/>
   <rotation name="rot08"      unit="deg" x="0" y="270" z="270"/>

   <rotation name="rot06"      unit="deg" x="180" y="270" z="0"/>
   <rotation name="rot05"      unit="deg" x="180" y="90" z="0"/>

   <rotation name="rPlus90AboutZ"      unit="deg" x="0" y="0" z="90"/>
   <rotation name="rMinus90AboutZ"      unit="deg" x="0" y="0" z="270"/>

    
    <rotation name="rBeamW3"             unit="deg" x="0" y="-16.1824077173347" z="43.7252002244654"/>



EOF

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
    $MAT = $basename."_Materials" . $suffix . ".gdml";
    push (@gdmlFiles, $MAT);
    $MAT = ">" . $MAT;

    open(MAT) or die("Could not open file $MAT for writing");

    # Add any materials special to this geometry by defining a mulitline string
    # and passing it to the gdmlMaterials::gen_Materials() function.
my $asmix = <<EOF;
  <!-- preliminary values -->
  <material name="AirSteelMixture" formula="AirSteelMixture">
   <D value=" @{[0.001205*(1-$FracMassOfSteel) + 7.9300*$FracMassOfSteel]} " unit="g/cm3"/>
   <fraction n="$FracMassOfSteel" ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
   <fraction n="$FracMassOfAir"   ref="Air"/>
  </material>
EOF

    # add the general materials used anywere
    print MAT gdmlMaterials::gen_Materials( $asmix );

    close(MAT);
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_ExtractionGrid +++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_ExtractionGrid {

    $ExtractionGrid = $basename."_ExtractionGrid" . $suffix . ".gdml";
    push (@gdmlFiles, $ExtractionGrid);
    $ExtractionGrid = ">" . $ExtractionGrid;
    open(ExtractionGrid) or die("Could not open file $ExtractionGrid for writing");

# The standard XML prefix and starting the gdml
print ExtractionGrid <<EOF;
<?xml version='1.0'?>
<gdml>
EOF

#GroundGrid SOLIDS
$ExtractionGridRadious = 0.05;
$ExtractionGridPitch = 0.3;

$ExtractionGridSizeX = 2*$ExtractionGridRadious;
$ExtractionGridSizeY = 600;
$ExtractionGridSizeZ = 600;


print ExtractionGrid <<EOF;

<solids>
      <tube name="solExtractionGridCable" rmin="0" rmax="$ExtractionGridRadious" z="$ExtractionGridSizeZ" deltaphi="360" startphi="0" aunit="deg" lunit="cm"/>
     <box name="solExtractionGrid" x="@{[$ExtractionGridSizeX]}" y="$ExtractionGridSizeY" z="@{[$ExtractionGridSizeZ]}" lunit="cm"/>
</solids>

EOF


print ExtractionGrid <<EOF;

<structure>

<volume name="volExtractionGridCable">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="solExtractionGridCable"/>
</volume>

<volume name="volExtractionGrid">
  <materialref ref="LAr"/>
  <solidref ref="solExtractionGrid"/>
EOF

for($ii=0;$ii<$$ExtractionGridSizeY;$ii=$ii+$ExtractionGridPitch)
{
	print ExtractionGrid <<EOF;
  <physvol>
   <volumeref ref="volExtractionGridCable"/>
   <position name="posExtractionGridCable$ii" unit="cm" x="0" y="@{[$ii-0.5*$ExtractionGridSizeY]}" z="0"/>
   <rotation name="GGrot$aux2" unit="deg" x="0" y="90" z="0" /> 
   </physvol>
EOF
 
}

for($jj=0;$jj<$$ExtractionGridSizeZ;$jj=$jj+$ExtractionGridPitch)
{
	print ExtractionGrid <<EOF;
  <physvol>
   <volumeref ref="volExtractionGridCable"/>
   <position name="posExtractionGridCableLat$jj" unit="cm" x="0" y="0" z="@{[$jj-0.5*$ExtractionGridSizeZ]}"/>
   <rotation name="GGrot$aux2" unit="deg" x="90" y="0" z="0" /> 
   </physvol>
EOF
 
}

	print ExtractionGrid <<EOF;
  
  </volume>
</structure>
</gdml>
EOF
close(ExtractionGrid);
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++ gen_TPC ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub gen_TPC()
{
    # CRM active volume
    my $TPCActive_x = $driftTPCActive;
    my $TPCActive_y = $widthCRM_active;
    my $TPCActive_z = $lengthCRM_active;

    # CRM total volume
    my $TPC_x = $TPCActive_x + $ReadoutPlane;
    my $TPC_y = $widthCRM;
    my $TPC_z = $lengthCRM;

	
    $TPC = $basename."_TPC" . $suffix . ".gdml";
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
   <box name="CRM" lunit="cm" 
      x="$TPC_x" 
      y="$TPC_y" 
      z="$TPC_z"/>
   <box name="CRMVPlane" lunit="cm" 
      x="$padWidth" 
      y="$TPCActive_y" 
      z="$TPCActive_z"/>
   <box name="CRMZPlane" lunit="cm" 
      x="$padWidth"
      y="$TPCActive_y"
      z="$TPCActive_z"/>
   <box name="CRMActive" lunit="cm"
      x="$TPCActive_x"
      y="$TPCActive_y"
      z="$TPCActive_z"/>
EOF


#++++++++++++++++++++++++++++ Wire Solids ++++++++++++++++++++++++++++++
# in principle we only need only one wire solid, since CRM is a square 
# but to be more general ...

print TPC <<EOF;

   <tube name="CRMWireV"
     rmax="@{[0.5*$padWidth]}"
     z="$TPCActive_z"               
     deltaphi="360"
     aunit="deg"
     lunit="cm"/>
   <tube name="CRMWireZ"
     rmax="@{[0.5*$padWidth]}"
     z="$TPCActive_y"               
     deltaphi="360"
     aunit="deg"
     lunit="cm"/>
</solids>

EOF


# Begin structure and create wire logical volumes
print TPC <<EOF;
<structure>

    <volume name="volTPCActive">
      <materialref ref="LAr"/>
      <solidref ref="CRMActive"/>
    </volume>
EOF

if ($wires_on==1) 
{ 
  print TPC <<EOF;
    <volume name="volTPCWireV">
      <materialref ref="Copper_Beryllium_alloy25"/>
      <solidref ref="CRMWireV"/>
    </volume>

    <volume name="volTPCWireZ">
      <materialref ref="Copper_Beryllium_alloy25"/>
      <solidref ref="CRMWireZ"/>
    </volume>
EOF
}

print TPC <<EOF;

    <volume name="volTPCPlaneV">
      <materialref ref="LAr"/>
      <solidref ref="CRMVPlane"/>
EOF

if ($wires_on==1) # add wires to Z plane
{
for($i=0;$i<$nChannelsViewPerCRM;++$i)
{
my $ypos = -0.5 * $TPCActive_y + ($i+0.5)*$wirePitch + 0.5*$padWidth;

print TPC <<EOF;
    <physvol>
      <volumeref ref="volTPCWireV"/> 
      <position name="posWireV$i" unit="cm" x="0" y="$ypos" z="0"/>
      <rotationref ref="rIdentity"/> 
    </physvol>
EOF
}
}

print TPC <<EOF;
   </volume>
   
   <volume name="volTPCPlaneZ">
     <materialref ref="LAr"/>
     <solidref ref="CRMZPlane"/>
EOF


if ($wires_on==1) # add wires to X plane
{
for($i=0;$i<$nChannelsViewPerCRM;++$i)
{

my $zpos = -0.5 * $TPCActive_z + ($i+0.5)*$wirePitch + 0.5*$padWidth;
print TPC <<EOF;
    <physvol>
     <volumeref ref="volTPCWireZ"/>
     <position name="posWireZ$i" unit="cm" x="0" y="0" z="$zpos"/>
     <rotationref ref="rPlus90AboutX"/>
    </physvol>
EOF
}
}


print TPC <<EOF;
   </volume>
EOF


$posVplane[0] = 0.5*$TPC_x - 1.5*$padWidth;
$posVplane[1] = 0;
$posVplane[2] = 0;

$posZplane[0] = 0.5*$TPC_x - 0.5*$padWidth;
$posZplane[1] = 0; 
$posZplane[2] = 0;

$posTPCActive[0] = -$ReadoutPlane;
$posTPCActive[1] = 0;
$posTPCActive[2] = 0;

$ExtractionGridX = 0.5*$Argon_x-$HeightGaseousAr-0.5-0.5*$ExtractionGridSizeX;
$ExtractionGridY = 0;
$ExtractionGridZ = 0;

$ExtractionGridX = 0.5*$TPC_x-0.5-0.5*$ExtractionGridSizeX;
$ExtractionGridY = 0;
$ExtractionGridZ = 0;

#wrap up the TPC file
print TPC <<EOF;

   <volume name="volTPC">
     <materialref ref="LAr"/>
     <solidref ref="CRM"/>
     <physvol>
       <volumeref ref="volTPCPlaneV"/>
       <position name="posPlaneV" unit="cm" 
         x="$posVplane[0]" y="$posVplane[1]" z="$posVplane[2]"/>
       <rotationref ref="rIdentity"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneZ"/>
       <position name="posPlaneZ" unit="cm" 
         x="$posZplane[0]" y="$posZplane[1]" z="$posZplane[2]"/>
       <rotationref ref="rIdentity"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCActive"/>
       <position name="posActive" unit="cm" 
         x="$posTPCActive[0]" y="$posTPCActive[1]" z="$posTPCActive[2]"/>
       <rotationref ref="rIdentity"/>
     </physvol>
EOF

  if ( $ExtractionGrid_switch eq "on" )
  {

      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volExtractionGrid"/>
   <position name="posExtractionGrid" unit="cm" x="$ExtractionGridX" y="$ExtractionGridY" z="$ExtractionGridZ"/>
  </physvol>
EOF

  }

print TPC <<EOF;
   </volume>
EOF


print TPC <<EOF;
</structure>
</gdml>
EOF

close(TPC);
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_FieldCage +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_FieldCage {


#FieldCage


    $FieldCage = $basename."_FieldCage" . $suffix . ".gdml";
    push (@gdmlFiles, $FieldCage);
    $FieldCage = ">" . $FieldCage;
    open(FieldCage) or die("Could not open file $FieldCage for writing");

# Create the <define> fragment file name, 
# add file to list of fragments,
# and open it




# The standard XML prefix and starting the gdml
print FieldCage <<EOF;
   <?xml version='1.0'?>
   <gdml>
EOF
# The printing solids used in the Field Cage
print FieldCage <<EOF;

    <solids>
        
	<xtru name="FIELD_SHAPER-SOL" lunit="mm">
	<twoDimVertex x="3.174889" y="4.182078"/>
	<twoDimVertex x="3.287192" y="4.163508"/>
	<twoDimVertex x="3.392385" y="4.120019"/>
	<twoDimVertex x="3.485016" y="4.053866"/>
	<twoDimVertex x="3.560284" y="3.968476"/>
	<twoDimVertex x="3.614289" y="3.868273"/>
	<twoDimVertex x="3.640993" y="3.777188"/>
	<twoDimVertex x="3.650003" y="3.682696"/>
	<twoDimVertex x="3.650003" y="0.765841000000002"/>
	<twoDimVertex x="3.635299" y="0.645475999999999"/>
	<twoDimVertex x="3.59205" y="0.532188999999999"/>
	<twoDimVertex x="3.522803" y="0.432644"/>
	<twoDimVertex x="3.431628" y="0.352698"/>
	<twoDimVertex x="3.323888" y="0.297053000000002"/>
	<twoDimVertex x="3.226003" y="0.271650999999999"/>
	<twoDimVertex x="3.125008" y="0.266466999999999"/>
	<twoDimVertex x="3.044094" y="0.263959"/>
	<twoDimVertex x="2.964646" y="0.248422999999999"/>
	<twoDimVertex x="2.852753" y="0.201684"/>
	<twoDimVertex x="2.755402" y="0.129382"/>
	<twoDimVertex x="2.67832" y="0.0357690000000019"/>
	<twoDimVertex x="2.626041" y="-0.0736460000000001"/>
	<twoDimVertex x="2.60164" y="-0.192428"/>
	<twoDimVertex x="2.606552" y="-0.313592"/>
	<twoDimVertex x="2.640488" y="-0.430011"/>
	<twoDimVertex x="2.701452" y="-0.534835000000001"/>
	<twoDimVertex x="2.785858" y="-0.6219"/>
	<twoDimVertex x="2.888742" y="-0.686086"/>
	<twoDimVertex x="3.004052" y="-0.723614999999999"/>
	<twoDimVertex x="3.125005" y="-0.732282000000001"/>
	<twoDimVertex x="3.245956" y="-0.740949000000001"/>
	<twoDimVertex x="3.361264" y="-0.778479999999998"/>
	<twoDimVertex x="3.464146" y="-0.842666000000001"/>
	<twoDimVertex x="3.548551" y="-0.929731"/>
	<twoDimVertex x="3.609514" y="-1.034553"/>
	<twoDimVertex x="3.639776" y="-1.131048"/>
	<twoDimVertex x="3.650003" y="-1.231656"/>
	<twoDimVertex x="3.650003" y="-4.051359"/>
	<twoDimVertex x="3.636188" y="-4.168081"/>
	<twoDimVertex x="3.595506" y="-4.278353"/>
	<twoDimVertex x="3.530206" y="-4.376082"/>
	<twoDimVertex x="3.443895" y="-4.455866"/>
	<twoDimVertex x="3.341345" y="-4.513298"/>
	<twoDimVertex x="3.247548" y="-4.54175"/>
	<twoDimVertex x="3.150003" y="-4.551358"/>
	<twoDimVertex x="2.599999" y="-4.551358"/>
	<twoDimVertex x="2.521782" y="-4.557514"/>
	<twoDimVertex x="2.44549" y="-4.57583"/>
	<twoDimVertex x="2.338749" y="-4.625038"/>
	<twoDimVertex x="2.246446" y="-4.697805"/>
	<twoDimVertex x="2.173679" y="-4.790109"/>
	<twoDimVertex x="2.124471" y="-4.89685"/>
	<twoDimVertex x="2.10154" y="-5.012128"/>
	<twoDimVertex x="2.106155" y="-5.129575"/>
	<twoDimVertex x="2.138059" y="-5.2427"/>
	<twoDimVertex x="2.19549" y="-5.345251"/>
	<twoDimVertex x="2.275275" y="-5.431561"/>
	<twoDimVertex x="2.373003" y="-5.496861"/>
	<twoDimVertex x="2.483276" y="-5.537543"/>
	<twoDimVertex x="2.599999" y="-5.551358"/>
	<twoDimVertex x="4.449998" y="-5.551358"/>
	<twoDimVertex x="4.547543" y="-5.54175"/>
	<twoDimVertex x="4.64134" y="-5.513298"/>
	<twoDimVertex x="4.743891" y="-5.455866"/>
	<twoDimVertex x="4.830201" y="-5.376082"/>
	<twoDimVertex x="4.895501" y="-5.278353"/>
	<twoDimVertex x="4.936183" y="-5.168081"/>
	<twoDimVertex x="4.949998" y="-5.051358"/>
	<twoDimVertex x="4.949998" y="3.827153"/>
	<twoDimVertex x="4.965356" y="3.950123"/>
	<twoDimVertex x="5.010485" y="4.065541"/>
	<twoDimVertex x="5.082615" y="4.166314"/>
	<twoDimVertex x="5.177314" y="4.246251"/>
	<twoDimVertex x="5.288764" y="4.300442"/>
	<twoDimVertex x="5.389517" y="4.323481"/>
	<twoDimVertex x="5.492853" y="4.325318"/>
	<twoDimVertex x="9.961787" y="3.782179"/>
	<twoDimVertex x="13.85859" y="2.875176"/>
	<twoDimVertex x="17.583122" y="1.413806"/>
	<twoDimVertex x="21.056972" y="-0.571162999999999"/>
	<twoDimVertex x="21.410387" y="-0.867622000000001"/>
	<twoDimVertex x="21.693333" y="-1.231946"/>
	<twoDimVertex x="21.893082" y="-1.647748"/>
	<twoDimVertex x="22.000649" y="-2.096324"/>
	<twoDimVertex x="22.011195" y="-2.557496"/>
	<twoDimVertex x="21.924246" y="-3.01052"/>
	<twoDimVertex x="21.743713" y="-3.435018"/>
	<twoDimVertex x="21.477717" y="-3.811897"/>
	<twoDimVertex x="21.138223" y="-4.124202"/>
	<twoDimVertex x="20.740501" y="-4.357886"/>
	<twoDimVertex x="20.302442" y="-4.502438"/>
	<twoDimVertex x="19.84375" y="-4.551358"/>
	<twoDimVertex x="14.900003" y="-4.551358"/>
	<twoDimVertex x="14.821786" y="-4.557514"/>
	<twoDimVertex x="14.745494" y="-4.57583"/>
	<twoDimVertex x="14.638754" y="-4.625038"/>
	<twoDimVertex x="14.54645" y="-4.697805"/>
	<twoDimVertex x="14.473683" y="-4.790109"/>
	<twoDimVertex x="14.424475" y="-4.89685"/>
	<twoDimVertex x="14.401544" y="-5.012128"/>
	<twoDimVertex x="14.406159" y="-5.129575"/>
	<twoDimVertex x="14.438063" y="-5.2427"/>
	<twoDimVertex x="14.495494" y="-5.345251"/>
	<twoDimVertex x="14.575279" y="-5.431561"/>
	<twoDimVertex x="14.673007" y="-5.496862"/>
	<twoDimVertex x="14.78328" y="-5.537543"/>
	<twoDimVertex x="14.900003" y="-5.551358"/>
	<twoDimVertex x="19.843727" y="-5.551355"/>
	<twoDimVertex x="20.513314" y="-5.479951"/>
	<twoDimVertex x="21.152782" y="-5.268942"/>
	<twoDimVertex x="21.733367" y="-4.927821"/>
	<twoDimVertex x="22.228955" y="-4.471929"/>
	<twoDimVertex x="22.617251" y="-3.921775"/>
	<twoDimVertex x="22.880791" y="-3.302104"/>
	<twoDimVertex x="23.007719" y="-2.640792"/>
	<twoDimVertex x="22.992327" y="-1.967586"/>
	<twoDimVertex x="22.835306" y="-1.312766"/>
	<twoDimVertex x="22.543721" y="-0.705788999999999"/>
	<twoDimVertex x="22.130685" y="-0.173957000000001"/>
	<twoDimVertex x="21.614779" y="0.258807000000001"/>
	<twoDimVertex x="18.014949" y="2.315762"/>
	<twoDimVertex x="14.155346" y="3.830128"/>
	<twoDimVertex x="10.117223" y="4.770024"/>
	<twoDimVertex x="4.51123" y="5.404618"/>
	<twoDimVertex x="-1.128658" y="5.551358"/>
	<twoDimVertex x="-10.117225" y="4.770024"/>
	<twoDimVertex x="-14.155346" y="3.830128"/>
	<twoDimVertex x="-18.014946" y="2.315762"/>
	<twoDimVertex x="-21.614774" y="0.258808999999999"/>
	<twoDimVertex x="-22.130682" y="-0.173954999999999"/>
	<twoDimVertex x="-22.543718" y="-0.705787000000001"/>
	<twoDimVertex x="-22.835302" y="-1.312765"/>
	<twoDimVertex x="-22.992326" y="-1.967584"/>
	<twoDimVertex x="-23.007719" y="-2.64079"/>
	<twoDimVertex x="-22.880788" y="-3.302103"/>
	<twoDimVertex x="-22.61725" y="-3.921773"/>
	<twoDimVertex x="-22.228952" y="-4.471928"/>
	<twoDimVertex x="-21.733366" y="-4.92782"/>
	<twoDimVertex x="-21.15278" y="-5.268942"/>
	<twoDimVertex x="-20.513313" y="-5.479951"/>
	<twoDimVertex x="-19.843727" y="-5.551358"/>
	<twoDimVertex x="-14.900001" y="-5.551358"/>
	<twoDimVertex x="-14.821785" y="-5.545203"/>
	<twoDimVertex x="-14.745493" y="-5.526886"/>
	<twoDimVertex x="-14.638752" y="-5.477678"/>
	<twoDimVertex x="-14.546447" y="-5.404911"/>
	<twoDimVertex x="-14.47368" y="-5.312607"/>
	<twoDimVertex x="-14.424472" y="-5.205866"/>
	<twoDimVertex x="-14.401542" y="-5.090588"/>
	<twoDimVertex x="-14.406156" y="-4.973141"/>
	<twoDimVertex x="-14.438062" y="-4.860016"/>
	<twoDimVertex x="-14.495493" y="-4.757466"/>
	<twoDimVertex x="-14.575277" y="-4.671155"/>
	<twoDimVertex x="-14.673006" y="-4.605854"/>
	<twoDimVertex x="-14.783279" y="-4.565174"/>
	<twoDimVertex x="-14.900001" y="-4.551358"/>
	<twoDimVertex x="-19.84375" y="-4.551355"/>
	<twoDimVertex x="-20.30244" y="-4.502438"/>
	<twoDimVertex x="-20.740499" y="-4.357886"/>
	<twoDimVertex x="-21.13822" y="-4.124202"/>
	<twoDimVertex x="-21.477714" y="-3.811899"/>
	<twoDimVertex x="-21.743711" y="-3.43502"/>
	<twoDimVertex x="-21.924244" y="-3.010522"/>
	<twoDimVertex x="-22.011194" y="-2.557498"/>
	<twoDimVertex x="-22.000646" y="-2.096328"/>
	<twoDimVertex x="-21.893081" y="-1.647752"/>
	<twoDimVertex x="-21.693334" y="-1.23195"/>
	<twoDimVertex x="-21.410389" y="-0.867626000000001"/>
	<twoDimVertex x="-21.056972" y="-0.571165999999998"/>
	<twoDimVertex x="-17.583122" y="1.413806"/>
	<twoDimVertex x="-13.858589" y="2.875176"/>
	<twoDimVertex x="-9.961782" y="3.782178"/>
	<twoDimVertex x="-5.492853" y="4.325312"/>
	<twoDimVertex x="-5.369017" y="4.320551"/>
	<twoDimVertex x="-5.250157" y="4.285479"/>
	<twoDimVertex x="-5.143572" y="4.222252"/>
	<twoDimVertex x="-5.055813" y="4.134755"/>
	<twoDimVertex x="-4.992267" y="4.028359"/>
	<twoDimVertex x="-4.960678" y="3.929951"/>
	<twoDimVertex x="-4.949997" y="3.827153"/>
	<twoDimVertex x="-4.949997" y="-5.051358"/>
	<twoDimVertex x="-4.940388" y="-5.148903"/>
	<twoDimVertex x="-4.911935" y="-5.2427"/>
	<twoDimVertex x="-4.854504" y="-5.345251"/>
	<twoDimVertex x="-4.77472" y="-5.431561"/>
	<twoDimVertex x="-4.676991" y="-5.496861"/>
	<twoDimVertex x="-4.566719" y="-5.537543"/>
	<twoDimVertex x="-4.449997" y="-5.551358"/>
	<twoDimVertex x="-2.599997" y="-5.551358"/>
	<twoDimVertex x="-2.52178" y="-5.545203"/>
	<twoDimVertex x="-2.445489" y="-5.526886"/>
	<twoDimVertex x="-2.338748" y="-5.477678"/>
	<twoDimVertex x="-2.246443" y="-5.404911"/>
	<twoDimVertex x="-2.173677" y="-5.312607"/>
	<twoDimVertex x="-2.124468" y="-5.205866"/>
	<twoDimVertex x="-2.101538" y="-5.090588"/>
	<twoDimVertex x="-2.106153" y="-4.973141"/>
	<twoDimVertex x="-2.138057" y="-4.860016"/>
	<twoDimVertex x="-2.195489" y="-4.757466"/>
	<twoDimVertex x="-2.275273" y="-4.671155"/>
	<twoDimVertex x="-2.373002" y="-4.605854"/>
	<twoDimVertex x="-2.483274" y="-4.565174"/>
	<twoDimVertex x="-2.599997" y="-4.551358"/>
	<twoDimVertex x="-3.150001" y="-4.551358"/>
	<twoDimVertex x="-3.266724" y="-4.537543"/>
	<twoDimVertex x="-3.376997" y="-4.496862"/>
	<twoDimVertex x="-3.474725" y="-4.431561"/>
	<twoDimVertex x="-3.55451" y="-4.345251"/>
	<twoDimVertex x="-3.611941" y="-4.2427"/>
	<twoDimVertex x="-3.640394" y="-4.148904"/>
	<twoDimVertex x="-3.650001" y="-4.051359"/>
	<twoDimVertex x="-3.650001" y="-1.231656"/>
	<twoDimVertex x="-3.635296" y="-1.111289"/>
	<twoDimVertex x="-3.592049" y="-0.998002"/>
	<twoDimVertex x="-3.522801" y="-0.89846"/>
	<twoDimVertex x="-3.431626" y="-0.818512999999999"/>
	<twoDimVertex x="-3.323887" y="-0.762865999999999"/>
	<twoDimVertex x="-3.226001" y="-0.737466000000001"/>
	<twoDimVertex x="-3.125007" y="-0.732282000000001"/>
	<twoDimVertex x="-3.049645" y="-0.722729999999999"/>
	<twoDimVertex x="-2.976599" y="-0.701872999999999"/>
	<twoDimVertex x="-2.875039" y="-0.650507999999999"/>
	<twoDimVertex x="-2.787724" y="-0.577507999999998"/>
	<twoDimVertex x="-2.719181" y="-0.486654000000002"/>
	<twoDimVertex x="-2.672959" y="-0.382652"/>
	<twoDimVertex x="-2.651454" y="-0.270893000000001"/>
	<twoDimVertex x="-2.65578" y="-0.157163000000001"/>
	<twoDimVertex x="-2.685712" y="-0.0473590000000002"/>
	<twoDimVertex x="-2.739699" y="0.0528300000000002"/>
	<twoDimVertex x="-2.814946" y="0.138217000000001"/>
	<twoDimVertex x="-2.907552" y="0.204374999999999"/>
	<twoDimVertex x="-3.012722" y="0.247876000000002"/>
	<twoDimVertex x="-3.125003" y="0.266466999999999"/>
	<twoDimVertex x="-3.245954" y="0.275136"/>
	<twoDimVertex x="-3.361263" y="0.312666"/>
	<twoDimVertex x="-3.464144" y="0.376850999999998"/>
	<twoDimVertex x="-3.548549" y="0.463915"/>
	<twoDimVertex x="-3.609512" y="0.568739000000001"/>
	<twoDimVertex x="-3.639774" y="0.665232"/>
	<twoDimVertex x="-3.650001" y="0.765841000000002"/>
	<twoDimVertex x="-3.650001" y="3.682696"/>
	<twoDimVertex x="-3.637045" y="3.795785"/>
	<twoDimVertex x="-3.598845" y="3.903012"/>
	<twoDimVertex x="-3.537384" y="3.998822"/>
	<twoDimVertex x="-3.455846" y="4.078246"/>
	<twoDimVertex x="-3.358455" y="4.137171"/>
	<twoDimVertex x="-3.268812" y="4.168376"/>
	<twoDimVertex x="-3.174888" y="4.182084"/>
	
	
	<section zOrder="0" zPosition="-1523.5077" xOffset="0" yOffset="0" scalingFactor="1"/>
	<section zOrder="1" zPosition="1431.2323" xOffset="0" yOffset="0" scalingFactor="1"/>
	<section zOrder="2" zPosition="1523.5077" xOffset="0" yOffset="@{[10*$FieldCageShaperProfileShift]}" scalingFactor="1"/>


	</xtru>
	
   <box name="FieldCageModule" lunit="cm"
      x="$FieldCageSizeX"
      y="$FieldCageSizeY"
      z="$FieldCageSizeZ"/>

    </solids>

EOF

print FieldCage <<EOF;

<structure>



  <volume name="volFIELDSHAPER">
  <materialref ref="ALUMINUM_Al"/>
  <solidref ref="FIELD_SHAPER-SOL"/>
  </volume>


  <volume name="volFieldCageModule">
  <materialref ref="LAr"/>
  <solidref ref="FieldCageModule"/>

EOF
    for ( $i=0; $i<$NFieldShapers; $i=$i+1 ) { # pmts with coating
   
      $Zposition = +0.5*$FieldCageSizeZ-0.5*$FFSSeparation-$i*$FFSSeparation;
      print FieldCage <<EOF;

   <physvol>
   <volumeref ref="volFIELDSHAPER"/>
   <position name="posFieldShaper$i" unit="cm" x="0" y="@{[-0.5*$FieldCageShaperProfileShift]}" z="$Zposition"/>
   <rotation name="FShaperRot$i" unit="deg" x="0" y="270" z="0"/>
  </physvol>

EOF
    }

print FieldCage <<EOF;

</volume>
</structure>
</gdml>
EOF
close(FieldCage);
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_GroundGrid +++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_GroundGrid {

    $GroundGrid = $basename."_GroundGrid" . $suffix . ".gdml";
    push (@gdmlFiles, $GroundGrid);
    $GroundGrid = ">" . $GroundGrid;
    open(GroundGrid) or die("Could not open file $GroundGrid for writing");

# The standard XML prefix and starting the gdml
print GroundGrid <<EOF;
<?xml version='1.0'?>
<gdml>
EOF

#GroundGrid SOLIDS
$GroundGridTubeLength = 297.5;
$GroundGridInnerRadious = 1.85;
$GroundGridOuterRadious = 2.0;
$GroundGridTorRad = 16;

$GroundGridInnerStructureLength = 313.5;
$GroundGridInnerStructureWidth = 2;
$GroundGridInnerStructureHeight = 4;
$GroundGridInnerStructureSeparation = 65.0;
$GroundGridInnerStructureNumberOfBars = 4;

$GroundGridInnerStructureNumberOfCablesPerInnerSquare = 5.0;
$GroundGridInnerStructureCableRadious = 0.1;
$GroundGridInnerStructureCableSeparation = $GroundGridInnerStructureSeparation/($GroundGridInnerStructureNumberOfCablesPerInnerSquare+1);

print GroundGrid <<EOF;

<solids>
     <torus name="GroundGridCorner" rmin="$GroundGridInnerRadious" rmax="$GroundGridOuterRadious" rtor="$GroundGridTorRad" deltaphi="90" startphi="0" aunit="deg" lunit="cm"/>
     <tube name="GroundGridtube" rmin="$GroundGridInnerRadious" rmax="$GroundGridOuterRadious" z="$GroundGridTubeLength" deltaphi="360" startphi="0" aunit="deg" lunit="cm"/>
     <box name="GroundGridInnerBoxBorder" x="@{[0.5*$GroundGridInnerStructureWidth]}" y="$GroundGridInnerStructureHeight" z="@{[$GroundGridInnerStructureLength]}" lunit="cm"/>
     <box name="GroundGridInnerBox" x="@{[$GroundGridInnerStructureWidth]}" y="$GroundGridInnerStructureHeight" z="@{[$GroundGridInnerStructureLength]}" lunit="cm"/>
    <box name="GroundGridModule" x="@{[$GroundGridInnerStructureLength+2+$GroundGridOuterRadious]}" y="$GroundGridInnerStructureHeight"    z="@{[$GroundGridInnerStructureLength+2+$GroundGridOuterRadious]}" lunit="cm"/>
     <tube name="GroundGridCable" rmin="0" rmax="$GroundGridInnerStructureCableRadious" z="@{[$GroundGridInnerStructureLength]}" deltaphi="360" startphi="0"  aunit="deg" lunit="cm"/>



    <union name="GGunion1">
      <first ref="GroundGridtube"/>
      <second ref="GroundGridCorner"/>
   		<position name="GGcorner1" unit="cm" x="@{[-$GroundGridTorRad]}" y="0" z="@{[0.5*$GroundGridTubeLength]}"/>
		<rotation name="GGrot1" unit="deg" x="90" y="0" z="0" />
    </union>

    <union name="GGunion2">
      <first ref="GGunion1"/>
      <second ref="GroundGridtube"/>
   		<position name="GGcorner2" unit="cm" x="@{[-0.5*$GroundGridTubeLength-$GroundGridTorRad]}" y="0" z="@{[+0.5*$GroundGridTubeLength+$GroundGridTorRad]}"/>
   		<rotation name="GGrot2" unit="deg" x="0" y="90" z="0" />
    </union>

    <union name="GGunion3">
      <first ref="GGunion2"/>
      <second ref="GroundGridInnerBoxBorder"/>
   		<position name="GGcorner3" unit="cm" x="@{[-$GroundGridTubeLength-$GroundGridTorRad+0.25*$GroundGridInnerStructureWidth]}" y="0" z="@{[$GroundGridTorRad-0.5*($GroundGridInnerStructureLength-$GroundGridTubeLength)]}"/>
    </union>


    <union name="GGunion4">
      <first ref="GGunion3"/>
      <second ref="GroundGridInnerBoxBorder"/>
   		<position name="GGcorner4" unit="cm" x="@{[-0.5*$GroundGridTubeLength-$GroundGridTorRad+0.5*($GroundGridInnerStructureLength-$GroundGridTubeLength)]}" y="0" z="@{[-0.5*$GroundGridTubeLength+0.25*$GroundGridInnerStructureWidth]}"/>
		<rotation name="GGrot4" unit="deg" x="0" y="90" z="0" />
    </union>

</solids>

EOF

$xGGorigin=0.5*($GroundGridInnerStructureLength+2+$GroundGridOuterRadious)-$GroundGridOuterRadious-2;
$zGGoriging=-0.5*($GroundGridInnerStructureLength+2+$GroundGridOuterRadious)+0.5*$GroundGridTubeLength;


print GroundGrid <<EOF;

<structure>

<volume name="volGroundGridCable">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="GroundGridCable"/>
</volume>

<volume name="volGroundGridInnerBox">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="GroundGridInnerBox"/>
</volume>

<volume name="volGGunion">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="GGunion4"/>
</volume>

 <volume name="volGroundGrid">
  <materialref ref="LAr"/>
  <solidref ref="GroundGridModule"/>

  <physvol>
   <volumeref ref="volGGunion"/>
   <position name="posGG18" unit="cm" x="@{[$xGGorigin]}" y="0" z="@{[$zGGoriging]}"/>
  </physvol>

EOF

$aux=4;
$aux2=$aux+1;
for($ii=0;$ii<$GroundGridInnerStructureNumberOfBars;$ii++)
{
	$aux2=$aux+1;
	print GroundGrid <<EOF;
  <physvol>
   <volumeref ref="volGroundGridInnerBox"/>
   <position name="posGGInnerBox$ii" unit="cm" x="@{[$xGGorigin-$GroundGridTubeLength-$GroundGridTorRad+$GroundGridInnerStructureWidth+($ii+1)*$GroundGridInnerStructureSeparation]}" y="0" z="@{[zGGoriging+$GroundGridTorRad-0.5*($GroundGridInnerStructureLength-$GroundGridTubeLength)-0.5*$GroundGridTorRad- 0.5*$GroundGridOuterRadious]}"/>
  </physvol>
EOF
   $aux++; 
}

for($ii=0;$ii<=$GroundGridInnerStructureNumberOfBars;$ii++)
{
 for($jj=0;$jj<$GroundGridInnerStructureNumberOfCablesPerInnerSquare;$jj++)
 {
	print GroundGrid <<EOF;
  <physvol>
   <volumeref ref="volGroundGridCable"/>
   <position name="posGGCable$ii$jj" unit="cm" x="@{[$xGGorigin-$GroundGridTubeLength-$GroundGridTorRad+$GroundGridInnerStructureWidth+($ii)*$GroundGridInnerStructureSeparation + ($jj+1)*$GroundGridInnerStructureCableSeparation]}" y="0" z="@{[zGGoriging+$GroundGridTorRad-0.5*($GroundGridInnerStructureLength-$GroundGridTubeLength)-0.5*$GroundGridTorRad- 0.5*$GroundGridOuterRadious]}"/>
  </physvol>
EOF
   if($ii ==$GroundGridInnerStructureNumberOfBars) { if($jj == 3){ $jj=$GroundGridInnerStructureNumberOfCablesPerInnerSquare;}}
 }
   
}

for($ii=0;$ii<$GroundGridInnerStructureNumberOfBars;$ii++)
{
	$aux2=$aux+1;
	print GroundGrid <<EOF;
  <physvol>
   <volumeref ref="volGroundGridInnerBox"/>
   <position name="posGGInnerBoxLat$ii" unit="cm" x="@{[$xGGorigin-0.5*$GroundGridTubeLength-$GroundGridTorRad+0.5*($GroundGridInnerStructureLength-$GroundGridTubeLength)]}" y="0" z="@{[zGGoriging-0.5*$GroundGridTubeLength+$GroundGridInnerStructureWidth +($ii+1)*$GroundGridInnerStructureSeparation - 0.5*$GroundGridTorRad- 0.5*$GroundGridOuterRadious]}"/>
   <rotation name="GGrot$aux2" unit="deg" x="0" y="90" z="0" /> 
   </physvol>
EOF
    $aux++;   
}
	print GroundGrid <<EOF;
  
  </volume>
</structure>
</gdml>
EOF
close(GroundGrid);
}




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_Cathode +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Cathode {

    $Cathode = $basename."_Cathode" . $suffix . ".gdml";
    push (@gdmlFiles, $Cathode);
    $Cathode = ">" . $Cathode;
    open(Cathode) or die("Could not open file $Cathode for writing");

# The standard XML prefix and starting the gdml
print Cathode <<EOF;
<?xml version='1.0'?>
<gdml>
EOF

#Cathode SOLIDS
$CathodeTubeLength = 297.5;
$CathodeInnerRadious = 1.85;
$CathodeOuterRadious = 2.0;
$CathodeTorRad = 16;

$CathodeInnerStructureLength = 313.5;
$CathodeInnerStructureWidth = 2;
$CathodeInnerStructureHeight = 4;
$CathodeInnerStructureSeparation = 65.0;
$CathodeInnerStructureNumberOfBars = 4;

$CathodeInnerStructureNumberOfCablesPerInnerSquare = 5.0;
$CathodeInnerStructureCableRadious = 0.5;
$CathodeInnerStructureCableSeparation = $CathodeInnerStructureSeparation/($CathodeInnerStructureNumberOfCablesPerInnerSquare+1);
$CathodeInnerCableLength = 128.0;
$CathodeInnerCableLengthShort = 47.5;#48.5;

$CathodeSupportTubeLength = 252.0;
$CathodeSupportInnerRadious = 1.85;
$CathodeSupportOuterRadious = 2.0;
$CathodeSupportTorRad = 12;

$CathodeSupportInnerStructureLength = 264.0;
$CathodeSupportInnerStructureWidth = 2;
$CathodeSupportInnerStructureHeight = 4;
$CathodeSupportInnerStructureSeparation = 65.0;
$CathodeSupportInnerStructureNumberOfBars = 3;

$CathodeSupportDistance =20;
$CathodeModuleX =$CathodeInnerStructureLength+2+$CathodeOuterRadious;
$CathodeModuleY =$CathodeInnerStructureHeight + $CathodeSupportDistance;
$CathodeModuleZ =$CathodeInnerStructureLength+2+$CathodeOuterRadious;

$CathodeSupportSeparatorsHeight=$CathodeSupportDistance-$CathodeInnerStructureHeight;

$CathodeSupportSeparatorsBorderLength=33;#53.39;
print Cathode <<EOF;

<solids>
     <torus name="CathodeCorner" rmin="$CathodeInnerRadious" rmax="$CathodeOuterRadious" rtor="$CathodeTorRad" deltaphi="90" startphi="0" aunit="deg" lunit="cm"/>
     <tube name="Cathodetube" rmin="$CathodeInnerRadious" rmax="$CathodeOuterRadious" z="$CathodeTubeLength" deltaphi="360" startphi="0" aunit="deg" lunit="cm"/>
     <box name="CathodeInnerBoxBorder" x="@{[0.5*$CathodeInnerStructureWidth]}" y="$CathodeInnerStructureHeight" z="@{[$CathodeInnerStructureLength]}" lunit="cm"/>
     <box name="CathodeInnerBox" x="@{[$CathodeInnerStructureWidth]}" y="$CathodeInnerStructureHeight" z="@{[$CathodeInnerStructureLength]}" lunit="cm"/>
    <box name="CathodeModule" x="@{[$CathodeModuleX]}" y="$CathodeModuleY"    z="@{[$CathodeModuleZ]}" lunit="cm"/>
    <box name="CathodeSupportSeparators" x="@{[$CathodeSupportInnerStructureWidth]}" y="@{[$CathodeSupportSeparatorsHeight]}"    z="@{[$CathodeSupportInnerStructureHeight]}" lunit="cm"/>
    <box name="CathodeSupportSeparatorsBorder" x="@{[$CathodeSupportInnerStructureWidth]}" y="$CathodeSupportSeparatorsBorderLength"    z="@{[$CathodeSupportInnerStructureHeight]}" lunit="cm"/>
     <tube name="CathodeCable" rmin="0" rmax="$CathodeInnerStructureCableRadious" z="@{[$CathodeInnerCableLength]}" deltaphi="360" startphi="0"  aunit="deg" lunit="cm"/>
     <tube name="CathodeCableShort" rmin="0" rmax="$CathodeInnerStructureCableRadious" z="@{[$CathodeInnerCableLengthShort]}" deltaphi="360" startphi="0"  aunit="deg" lunit="cm"/>
     <tube name="CathodeCableShortShort" rmin="0" rmax="$CathodeInnerStructureCableRadious" z="@{[$CathodeInnerCableLengthShort-3]}" deltaphi="360" startphi="0"  aunit="deg" lunit="cm"/>



    <union name="Cathunion1">
      <first ref="Cathodetube"/>
      <second ref="CathodeCorner"/>
   		<position name="Cathcorner1" unit="cm" x="@{[-$CathodeTorRad]}" y="0" z="@{[0.5*$CathodeTubeLength]}"/>
		<rotation name="Cathrot1" unit="deg" x="90" y="0" z="0" />
    </union>

    <union name="Cathunion2">
      <first ref="Cathunion1"/>
      <second ref="Cathodetube"/>
   		<position name="Cathcorner2" unit="cm" x="@{[-0.5*$CathodeTubeLength-$CathodeTorRad]}" y="0" z="@{[+0.5*$CathodeTubeLength+$CathodeTorRad]}"/>
   		<rotation name="Cathrot2" unit="deg" x="0" y="90" z="0" />
    </union>

    <union name="Cathunion3">
      <first ref="Cathunion2"/>
      <second ref="CathodeInnerBoxBorder"/>
   		<position name="Cathcorner3" unit="cm" x="@{[-$CathodeTubeLength-$CathodeTorRad+0.25*$CathodeInnerStructureWidth]}" y="0" z="@{[$CathodeTorRad-0.5*($CathodeInnerStructureLength-$CathodeTubeLength)]}"/>
    </union>


    <union name="Cathunion4">
      <first ref="Cathunion3"/>
      <second ref="CathodeInnerBoxBorder"/>
   		<position name="Cathcorner4" unit="cm" x="@{[-0.5*$CathodeTubeLength-$CathodeTorRad+0.5*($CathodeInnerStructureLength-$CathodeTubeLength)]}" y="0" z="@{[-0.5*$CathodeTubeLength+0.25*$CathodeInnerStructureWidth]}"/>
		<rotation name="Cathrot4" unit="deg" x="0" y="90" z="0" />
    </union>


     <torus name="CathodeSupportCorner" rmin="$CathodeSupportInnerRadious" rmax="$CathodeSupportOuterRadious" rtor="$CathodeSupportTorRad" deltaphi="90" startphi="0" aunit="deg" lunit="cm"/>
     <tube name="CathodeSupporttube" rmin="$CathodeSupportInnerRadious" rmax="$CathodeSupportOuterRadious" z="$CathodeSupportTubeLength" deltaphi="360" startphi="0" aunit="deg" lunit="cm"/>
     <box name="CathodeSupportInnerBoxBorder" x="@{[0.5*$CathodeSupportInnerStructureWidth]}" y="$CathodeSupportInnerStructureHeight" z="@{[$CathodeSupportInnerStructureLength]}" lunit="cm"/>
     <box name="CathodeSupportInnerBox" x="@{[$CathodeSupportInnerStructureWidth]}" y="$CathodeSupportInnerStructureHeight" z="@{[$CathodeSupportInnerStructureLength]}" lunit="cm"/>
    <box name="CathodeSupportModule" x="@{[$CathodeSupportInnerStructureLength+2+$CathodeSupportOuterRadious]}" y="$CathodeSupportInnerStructureHeight"    z="@{[$CathodeSupportInnerStructureLength+2+$CathodeSupportOuterRadious]}" lunit="cm"/>



    <union name="CathSupunion1">
      <first ref="CathodeSupporttube"/>
      <second ref="CathodeSupportCorner"/>
   		<position name="CathSupcorner1" unit="cm" x="@{[-$CathodeSupportTorRad]}" y="0" z="@{[0.5*$CathodeSupportTubeLength]}"/>
		<rotation name="CathSuprot1" unit="deg" x="90" y="0" z="0" />
    </union>

    <union name="CathSupunion2">
      <first ref="CathSupunion1"/>
      <second ref="CathodeSupporttube"/>
   		<position name="CathSupcorner2" unit="cm" x="@{[-0.5*$CathodeSupportTubeLength-$CathodeSupportTorRad]}" y="0" z="@{[+0.5*$CathodeSupportTubeLength+$CathodeSupportTorRad]}"/>
   		<rotation name="CathSuprot2" unit="deg" x="0" y="90" z="0" />
    </union>

    <union name="CathSupunion3">
      <first ref="CathSupunion2"/>
      <second ref="CathodeSupportInnerBoxBorder"/>
   		<position name="CathSupcorner3" unit="cm" x="@{[-$CathodeSupportTubeLength-$CathodeSupportTorRad+0.25*$CathodeSupportInnerStructureWidth]}" y="0" z="@{[$CathodeSupportTorRad-0.5*($CathodeSupportInnerStructureLength-$CathodeSupportTubeLength)]}"/>
    </union>


    <union name="CathSupunion4">
      <first ref="CathSupunion3"/>
      <second ref="CathodeSupportInnerBoxBorder"/>
   		<position name="CathSupcorner4" unit="cm" x="@{[-0.5*$CathodeSupportTubeLength-$CathodeSupportTorRad+0.5*($CathodeSupportInnerStructureLength-$CathodeSupportTubeLength)]}" y="0" z="@{[-0.5*$CathodeSupportTubeLength+0.25*$CathodeSupportInnerStructureWidth]}"/>
		<rotation name="CathSuprot4" unit="deg" x="0" y="90" z="0" />
    </union>

</solids>

EOF


$xCathorigin=0.5*($CathodeModuleX)-$CathodeOuterRadious-2;
$yCathorigin=0.5*$CathodeModuleY - 0.5*$CathodeInnerStructureHeight;
$zCathoriging=-0.5*($CathodeModuleX)+0.5*$CathodeTubeLength;
$xCathSuporigin=$xCathorigin+($CathodeSupportInnerStructureLength-$CathodeInnerStructureLength);
$zCathSuporigin=-$zCathoriging-($CathodeInnerStructureLength+$CathodeOuterRadious-$CathodeSupportInnerStructureLength-$CathodeSupportOuterRadious)+7;

$xCathSuporigin=$xCathorigin+($CathodeSupportInnerStructureLength-$CathodeInnerStructureLength);
$zCathSuporigin=-0.5*($CathodeModuleZ)+0.5*$CathodeSupportTubeLength;
$yCathSuporigin=$yCathorigin-$CathodeSupportDistance;


print Cathode <<EOF;

<structure>



<volume name="volCathodeSupportSeparators">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="CathodeSupportSeparators"/>
</volume>

<volume name="volCathodeSupportSeparatorsBorder">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="CathodeSupportSeparatorsBorder"/>
</volume>




<volume name="volCathodeCable">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="CathodeCable"/>
</volume>

<volume name="volCathodeCableShort">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="CathodeCableShort"/>
</volume>

<volume name="volCathodeCableShortShort">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="CathodeCableShortShort"/>
</volume>

<volume name="volCathodeInnerBox">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="CathodeInnerBox"/>
</volume>

<volume name="volCathodeSupportInnerBox">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="CathodeSupportInnerBox"/>
</volume>


<volume name="volCathunion">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="Cathunion4"/>
</volume>

<volume name="volCathSupunion">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="CathSupunion4"/>
</volume>

 <volume name="volCathode">
  <materialref ref="LAr"/>
  <solidref ref="CathodeModule"/>

  <physvol>
   <volumeref ref="volCathunion"/>
   <position name="posCath18" unit="cm" x="@{[$xCathorigin]}" y="$yCathorigin" z="@{[$zCathoriging]}"/>
  </physvol>

  <physvol>
   <volumeref ref="volCathSupunion"/>
   <position name="posCathSup" unit="cm" x="@{[$xCathSuporigin]}" y="$yCathSuporigin" z="@{[$zCathSuporigin]}"/>
  </physvol>

EOF

$aux=4;
$aux2=$aux+1;
for($ii=0;$ii<$CathodeInnerStructureNumberOfBars;$ii++)
{
	$aux2=$aux+1;
	print Cathode <<EOF;
  <physvol>
   <volumeref ref="volCathodeInnerBox"/>
   <position name="posCathInnerBox$ii" unit="cm" x="@{[$xCathorigin-$CathodeTubeLength-$CathodeTorRad+$CathodeInnerStructureWidth+($ii+1)*$CathodeInnerStructureSeparation]}" y="$yCathorigin" z="@{[zCathoriging+$CathodeTorRad-0.5*($CathodeInnerStructureLength-$CathodeTubeLength)-0.5*$CathodeTorRad- 0.5*$CathodeOuterRadious]}"/>
  </physvol>
EOF
   $aux++; 
}



for($ii=0;$ii<$CathodeInnerStructureNumberOfBars;$ii++)
{
	$aux2=$aux+1;
	print Cathode <<EOF;
  <physvol>
   <volumeref ref="volCathodeInnerBox"/>
   <position name="posCathInnerBoxLat$ii" unit="cm" x="@{[$xCathorigin-0.5*$CathodeTubeLength-$CathodeTorRad+0.5*($CathodeInnerStructureLength-$CathodeTubeLength)]}" y="$yCathorigin" z="@{[zCathoriging-0.5*$CathodeTubeLength+$CathodeInnerStructureWidth +($ii+1)*$CathodeInnerStructureSeparation - 0.5*$CathodeTorRad- 0.5*$CathodeOuterRadious]}"/>
   <rotation name="Cathrot$aux2" unit="deg" x="0" y="90" z="0" /> 
   </physvol>
EOF
    $aux++;   
}


for($kk=0;$kk<$CathodeInnerStructureNumberOfBars;$kk=$kk+2)
{
for($ii=0;$ii<=$CathodeInnerStructureNumberOfBars;$ii=$ii+2)
{
if ($kk==2) {$ii=$ii+2;};

 for($jj=0;$jj<2*$CathodeInnerStructureNumberOfCablesPerInnerSquare+1;$jj++)
 {
	print Cathode <<EOF;
  <physvol>
   <volumeref ref="volCathodeCable"/>
   <position name="posCathCable$ii$jj$kk" unit="cm" x="@{[($xCathorigin-$CathodeTubeLength-$CathodeTorRad+$CathodeInnerStructureWidth+($ii)*$CathodeInnerStructureSeparation+($jj+1)*$CathodeInnerStructureCableSeparation)]}" y="$yCathorigin" z="@{[zCathoriging-0.5*$CathodeTubeLength+$CathodeInnerStructureWidth +($kk+1)*$CathodeInnerStructureSeparation - 0.5*$CathodeTorRad- 0.5*$CathodeOuterRadious]}"/>
 </physvol>
EOF
   if($ii ==$CathodeInnerStructureNumberOfCablesPerInnerSquare-1 && $jj ==$CathodeInnerStructureNumberOfCablesPerInnerSquare-2) {$jj=2*$CathodeInnerStructureNumberOfCablesPerInnerSquare+1;}
   if($jj ==$CathodeInnerStructureNumberOfCablesPerInnerSquare-1) {$jj++;}

 }
if ($kk==2) {$ii=$ii+2;};
if ($kk==0) {$ii=$ii+2;};

}
}


$kk=4;
for($ii=0;$ii<=$CathodeInnerStructureNumberOfBars;$ii=$ii+4)
{
 for($jj=0;$jj<2*$CathodeInnerStructureNumberOfCablesPerInnerSquare+1;$jj++)
 {
	print Cathode <<EOF;
  <physvol>
   <volumeref ref="volCathodeCableShort"/>
   <position name="posCathCableShort$ii$jj$kk" unit="cm" x="@{[($xCathorigin-$CathodeTubeLength-$CathodeTorRad+$CathodeInnerStructureWidth+($ii)*$CathodeInnerStructureSeparation+($jj+1)*$CathodeInnerStructureCableSeparation)]}" y="$yCathorigin" z="@{[zCathoriging-0.5*$CathodeTubeLength+$CathodeInnerStructureWidth +($kk)*$CathodeInnerStructureSeparation - 0.5*$CathodeTorRad- 0.5*$CathodeOuterRadious+0.5*$CathodeInnerCableLengthShort+1]}"/>
 </physvol>
EOF
   if($ii ==$CathodeInnerStructureNumberOfCablesPerInnerSquare-1 && $jj ==$CathodeInnerStructureNumberOfCablesPerInnerSquare-3) {$jj=2*$CathodeInnerStructureNumberOfCablesPerInnerSquare+1;}
   if($jj ==$CathodeInnerStructureNumberOfCablesPerInnerSquare-1) {$jj++;}
   #Falta el $ii 4 $jj 2 that needs to be shorter to not overlap.

#print "$ii $jj \n";

 }
}

$ii=4;
$jj=3;

	print Cathode <<EOF;
  <physvol>
   <volumeref ref="volCathodeCableShortShort"/>
   <position name="posCathCableShort$ii$jj$kk" unit="cm" x="@{[($xCathorigin-$CathodeTubeLength-$CathodeTorRad+$CathodeInnerStructureWidth+($ii)*$CathodeInnerStructureSeparation+($jj+1)*$CathodeInnerStructureCableSeparation)]}" y="$yCathorigin" z="@{[zCathoriging-0.5*$CathodeTubeLength+$CathodeInnerStructureWidth +($kk)*$CathodeInnerStructureSeparation - 0.5*$CathodeTorRad- 0.5*$CathodeOuterRadious+0.5*$CathodeInnerCableLengthShort+1-1]}"/>
 </physvol>
EOF

for($ii=0;$ii<$CathodeInnerStructureNumberOfBars;$ii=$ii+2)
{
for($kk=0;$kk<=$CathodeInnerStructureNumberOfBars;$kk=$kk+2)
{
if ($ii==0) {$kk=$kk+2;};

 for($jj=0;$jj<2*$CathodeInnerStructureNumberOfCablesPerInnerSquare+1;$jj++)
 {
	print Cathode <<EOF;
  <physvol>
   <volumeref ref="volCathodeCable"/>
   <position name="posCathCable$ii$jj$kk" unit="cm" x="@{[$xCathorigin-$CathodeTubeLength-$CathodeTorRad+$CathodeInnerStructureWidth+($ii+1)*$CathodeInnerStructureSeparation]}" y="$yCathorigin" z="@{[zCathoriging-0.5*$CathodeTubeLength+$CathodeInnerStructureWidth +($kk)*$CathodeInnerStructureSeparation - 0.5*$CathodeTorRad- 0.5*$CathodeOuterRadious +($jj+1)*$CathodeInnerStructureCableSeparation]}"/>
   <rotation name="CathCablerot$ii$jj$kk" unit="deg" x="0" y="90" z="0" /> 
 </physvol>
EOF
   if($kk ==$CathodeInnerStructureNumberOfCablesPerInnerSquare-1 && $jj ==$CathodeInnerStructureNumberOfCablesPerInnerSquare-2) {$jj=2*$CathodeInnerStructureNumberOfCablesPerInnerSquare+1;}
   if($jj ==$CathodeInnerStructureNumberOfCablesPerInnerSquare-1) {$jj++;}

 }
if ($ii==2) {$kk=$kk+2;};
if ($ii==0) {$kk=$kk+2;};

}
}

$ii=4;
$kk=2;
 for($jj=0;$jj<2*$CathodeInnerStructureNumberOfCablesPerInnerSquare+1;$jj++)
 {
	print Cathode <<EOF;
  <physvol>
   <volumeref ref="volCathodeCableShort"/>
   <position name="posCathCableShort$ii$jj$kk" unit="cm" x="@{[$xCathorigin-$CathodeTubeLength-$CathodeTorRad+$CathodeInnerStructureWidth+($ii)*$CathodeInnerStructureSeparation+0.5*$CathodeInnerCableLengthShort+1]}" y="$yCathorigin" z="@{[zCathoriging-0.5*$CathodeTubeLength+$CathodeInnerStructureWidth +($kk)*$CathodeInnerStructureSeparation - 0.5*$CathodeTorRad- 0.5*$CathodeOuterRadious +($jj+1)*$CathodeInnerStructureCableSeparation]}"/>
   <rotation name="CathCablerot$ii$jj$kk" unit="deg" x="0" y="90" z="0" /> 
 </physvol>
EOF

   if($jj ==$CathodeInnerStructureNumberOfCablesPerInnerSquare-1) {$jj++;}
 }


for($ii=0;$ii<$CathodeSupportInnerStructureNumberOfBars;$ii++)
{
	print Cathode <<EOF;
  <physvol>
   <volumeref ref="volCathodeSupportInnerBox"/>
   <position name="posCathSupInnerBoxLat$ii" unit="cm" x="@{[$xCathSuporigin-0.5*$CathodeSupportTubeLength-$CathodeSupportTorRad+0.5*($CathodeSupportInnerStructureLength-$CathodeSupportTubeLength)]}" y="$yCathSuporigin" z="@{[zCathoriging-0.5*$CathodeTubeLength+$CathodeInnerStructureWidth +($ii+1)*$CathodeInnerStructureSeparation - 0.5*$CathodeTorRad- 0.5*$CathodeOuterRadious]}"/>
   <rotation name="CathSupInnerBoxLatrot$ii" unit="deg" x="0" y="90" z="0" /> 
   </physvol>
EOF
}

for($ii=0;$ii<$CathodeSupportInnerStructureNumberOfBars;$ii++)
{
	$aux2=$aux+1;
	print Cathode <<EOF;
  <physvol>
   <volumeref ref="volCathodeSupportInnerBox"/>
   <position name="posCathSupInnerBox$ii" unit="cm" x="@{[$xCathSuporigin-$CathodeSupportTubeLength-$CathodeSupportTorRad+$CathodeSupportInnerStructureWidth+($ii+1)*$CathodeSupportInnerStructureSeparation]}" y="$yCathSuporigin" z="@{[$zCathSuporigin+$CathodeSupportTorRad-$CathodeSupportOuterRadious-4]}"/>
  </physvol>
EOF
   $aux++; 
}

for($ii=0;$ii<=$CathodeSupportInnerStructureNumberOfBars+1;$ii++)
{
   for($jj=0;$jj<=$CathodeSupportInnerStructureNumberOfBars+1;$jj++)
   {	print Cathode <<EOF;
    <physvol>
   <volumeref ref="volCathodeSupportSeparators"/>
   <position name="posCathodeSupportSeparators$ii$jj" unit="cm" x="@{[$xCathSuporigin-$CathodeSupportTubeLength-$CathodeSupportTorRad+$CathodeSupportInnerStructureWidth+($ii)*$CathodeSupportInnerStructureSeparation]}" y="@{[0.5*$CathodeSupportSeparatorsHeight+$yCathSuporigin+2]}" z="@{[zCathoriging-0.5*$CathodeTubeLength+$CathodeInnerStructureWidth +($jj)*$CathodeInnerStructureSeparation - 0.5*$CathodeTorRad- 0.5*$CathodeOuterRadious]}"/>
  </physvol>
EOF
	if(($jj==$CathodeSupportInnerStructureNumberOfBars+1 )&& $jj!=$ii)
	{
		
	print Cathode <<EOF;
    <physvol>
   <volumeref ref="volCathodeSupportSeparatorsBorder"/>
   <position name="posCathodeSupportSeparatorsBorder$ii$jj" unit="cm" x="@{[$xCathSuporigin-$CathodeSupportTubeLength-$CathodeSupportTorRad+$CathodeSupportInnerStructureWidth+($ii)*$CathodeSupportInnerStructureSeparation]}" y="@{[0.5*$CathodeSupportDistance+$yCathSuporigin]}" z="@{[zCathoriging-0.5*$CathodeTubeLength+$CathodeInnerStructureWidth +($jj)*$CathodeInnerStructureSeparation - 0.5*$CathodeTorRad- 0.5*$CathodeOuterRadious+0.5*49.5]}"/>
   <rotation name="CathSupSepRot$ii$jj" unit="deg" x="292" y="0" z="0" /> 
  </physvol>
EOF
		
	}

	if(($ii==$CathodeSupportInnerStructureNumberOfBars+1)&& $jj!=$ii)
	{
		
	print Cathode <<EOF;
    <physvol>
   <volumeref ref="volCathodeSupportSeparatorsBorder"/>
   <position name="posCathodeSupportSeparatorsBorder$ii$jj" unit="cm" x="@{[$xCathSuporigin-$CathodeSupportTubeLength-$CathodeSupportTorRad+$CathodeSupportInnerStructureWidth+($ii)*$CathodeSupportInnerStructureSeparation+0.5*49.5]}" y="@{[0.5*$CathodeSupportDistance+$yCathSuporigin]}" z="@{[zCathoriging-0.5*$CathodeTubeLength+$CathodeInnerStructureWidth +($jj)*$CathodeInnerStructureSeparation - 0.5*$CathodeTorRad- 0.5*$CathodeOuterRadious]}"/>
   <rotation name="CathSupSepRot$ii$jj" unit="deg" x="0" y="0" z="68" /> 

  </physvol>
EOF
		
	}


	if($ii==$CathodeSupportInnerStructureNumberOfBars+1){if($jj==$CathodeSupportInnerStructureNumberOfBars) {$j++;}}
   }
}

	print Cathode <<EOF;
  
  </volume>
</structure>
</gdml>
EOF
close(Cathode);
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_LEMs +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_LEMs {


    $LEMs = $basename."_LEMs" . $suffix . ".gdml";
    push (@gdmlFiles, $LEMs);
    $LEMs = ">" . $LEMs;
    open(LEMs) or die("Could not open file $LEMs for writing");

# The standard XML prefix and starting the gdml
print LEMs <<EOF;
<?xml version='1.0'?>
<gdml>
EOF

$LEMsSizeX=0.1;
$LEMsSizeY=600.64;
$LEMsSizeZ=600.64;

print LEMs <<EOF;

<solids>
     <box name="solLEMs" x="@{[$LEMsSizeX]}" y="$LEMsSizeY" z="@{[$LEMsSizeZ]}" lunit="cm"/>

</solids>

EOF
print LEMs <<EOF;


<structure>
 <volume name="volLEMs">
  <materialref ref="Copper_Beryllium_alloy25"/>
  <solidref ref="solLEMs"/>
 </volume>
</structure>
</gdml>
EOF
close(LEMs);
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_pmt +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_pmt {

#PMTs
#$PMT_COATING_THICKNESS=0.2;#
#$PMT_PLATE_THICKNESS=0.4;
#$PMT_GLASS_THICKNESS=0.2;

    $PMT = $basename."_PMT" . $suffix . ".gdml";
    push (@gdmlFiles, $PMT);
    $PMT = ">" . $PMT;
    open(PMT) or die("Could not open file $PMT for writing");

# The standard XML prefix and starting the gdml
print PMT <<EOF;
<?xml version='1.0'?>
<gdml>
EOF

#PMT SOLIDS
print PMT <<EOF;

<solids>
 <tube name="PMTVolume"
  rmax="(6.5*2.54)"
  z="(11.1*2.54)"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>

    <tube aunit="deg" deltaphi="360" lunit="mm" name="pmtMiddleCylinder" rmax="102.351822048586" rmin="100.351822048586" startphi="0" z="54"/>
    <sphere aunit="deg" deltaphi="360" deltatheta="50" lunit="mm" name="sphPartTop" rmax="133" rmin="131" startphi="0" starttheta="0"/>
    <union name="pmt0x7fb8f489dfe0">
      <first ref="pmtMiddleCylinder"/>
      <second ref="sphPartTop"/>
      <position name="pmt0x7fb8f489dfe0_pos" unit="mm" x="0" y="0" z="-57.2051768689367"/>
    </union>
    <sphere aunit="deg" deltaphi="360" deltatheta="31.477975238527" lunit="mm" name="sphPartBtm" rmax="133" rmin="131" startphi="0" starttheta="130"/>
    <union name="pmt0x7fb8f48a0d50">
      <first ref="pmt0x7fb8f489dfe0"/>
      <second ref="sphPartBtm"/>
      <position name="pmt0x7fb8f48a0d50_pos" unit="mm" x="0" y="0" z="57.2051768689367"/>
    </union>
    <tube aunit="deg" deltaphi="360" lunit="mm" name="pmtBtmTube" rmax="44.25" rmin="42.25" startphi="0" z="72"/>
    <union name="solidpmt">
      <first ref="pmt0x7fb8f48a0d50"/>
      <second ref="pmtBtmTube"/>
      <position name="solidpmt_pos" unit="mm" x="0" y="0" z="-104.905637496842"/>
    </union>
    <sphere aunit="deg" deltaphi="360" deltatheta="50" lunit="mm" name="pmt0x7fb8f48a1eb0" rmax="133.2" rmin="133" startphi="0" starttheta="0"/>
    <sphere aunit="deg" deltaphi="360" deltatheta="46.5" lunit="mm" name="pmt0x7fb8f48a4860" rmax="131" rmin="130.999" startphi="0" starttheta="0"/>


</solids>

EOF
print PMT <<EOF;


<structure>
 <volume name="pmtCoatVol">
  <materialref ref="LAr"/>
  <solidref ref="pmt0x7fb8f48a1eb0"/>
  </volume>

 <volume name="allpmt">
  <materialref ref="Glass"/>
  <solidref ref="solidpmt"/>
  </volume>

<volume name="volPMT_coated">
  <materialref ref="LAr"/>
  <solidref ref="PMTVolume"/>

  <physvol>
   <volumeref ref="allpmt"/>
   <position name="posallpmt" unit="cm" x="0" y="0" z="1.27*2.54"/>
  </physvol>

 <physvol name="volOpDetSensitive">
  <volumeref ref="pmtCoatVol"/>
  <position name="posOpDetSensitive" unit="cm" x="0" y="0" z="1.27*2.54- (2.23*2.54)"/>
  </physvol>

 </volume>
</structure>
</gdml>
EOF
close(PMT);
}

  #<physvol>
  #<volumeref ref="pmtCathVol"/>
  # <position name="pmtCath_pos" unit="cm" x="0" y="0" z="0"/>#8*2.54 + 0.9"/>
  #</physvol>
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_Cryostat +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Cryostat()
{

# Create the cryostat fragment file name,
# add file to list of output GDML fragments,
# and open it
    $CRYO = $basename."_Cryostat" . $suffix . ".gdml";
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

    <box name="GaseousArgon" lunit="cm" 
      x="$HeightGaseousAr"
      y="$Argon_y"
      z="$Argon_z"/>

    <subtraction name="SteelShell">
      <first ref="Cryostat"/>
      <second ref="ArgonInterior"/>
    </subtraction>

</solids>
EOF

# Cryostat structure
print CRYO <<EOF;
<structure>
    <volume name="volSteelShell">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="SteelShell" />
    </volume>
    <volume name="volGaseousArgon">
      <materialref ref="ArGas"/>
      <solidref ref="GaseousArgon"/>

EOF
  if ( $LEMs_switch eq "on" )
  {

$posLEMsX = -0.5*$HeightGaseousAr+0.5;
$posLEMsY = 0;
$posLEMsZ = 0;

      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volLEMs"/>
   <position name="posLEMs" unit="cm" x="$posLEMsX" y="$posLEMsY" z="$posLEMsZ"/>

  </physvol>
EOF

  }
      print CRYO <<EOF;
    </volume>

    <volume name="volCryostat">
      <materialref ref="LAr" />
      <solidref ref="Cryostat" />
      <physvol>
        <volumeref ref="volGaseousArgon"/>
        <position name="posGaseousArgon" unit="cm" x="$Argon_x/2-$HeightGaseousAr/2" y="0" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volSteelShell"/>
        <position name="posSteelShell" unit="cm" x="0" y="0" z="0"/>
      </physvol>
EOF


if ($tpc_on==1) # place TPC inside croysotat
{
$posX =  $Argon_x/2 - $HeightGaseousAr - 0.5*($driftTPCActive + $ReadoutPlane); 

for($ii=0;$ii<$nCRM_z;$ii++)
{
    $posZ = -0.5*$Argon_z + $zLArBuffer + ($ii+0.5)*$lengthCRM;

    for($jj=0;$jj<$nCRM_y;$jj++)
    {
	$posY = -0.5*$Argon_y + $yLArBuffer + ($jj+0.5)*$widthCRM;
	print CRYO <<EOF;
      <physvol>
        <volumeref ref="volTPC"/>
	<position name="posTPC\-$ii\-$jj" unit="cm"
           x="$posX" y="$posY" z="$posZ"/>
      </physvol>
EOF
    }
}

}


   
  if ( $pmt_switch eq "on" ) {
    $pmt_pos_x =  -$Argon_x/2 + 0.5*($HeightPMT);
    for ( $i=0; $i<36; $i=$i+1 ) { # pmts with coating
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volPMT_coated"/>
   <position name="posPMT$i" unit="cm" x="$pmt_pos_x" @pmt_pos[$i]/>
   <rotationref ref="rMinus90AboutY"/>
  </physvol>
EOF
    }

  }

$GroundGridPosX=-$Argon_x/2 + 47.5;
$GroundGridPosY=0.5*($GroundGridInnerStructureLength+2+$GroundGridOuterRadious);

  if ( $GroundGrid_switch eq "on" )
  {
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volGroundGrid"/>
   <position name="posGroundGrid01" unit="cm" x="$GroundGridPosX" y="@{[-$GroundGridPosY]}" z="@{[$GroundGridPosY]}"/>
   <rotation name="rotGG01" unit="deg" x="0" y="0" z="90" />
  </physvol>

  <physvol>
   <volumeref ref="volGroundGrid"/>
   <position name="posGroundGrid02" unit="cm" x="$GroundGridPosX" y="@{[-$GroundGridPosY]}" z="@{[-$GroundGridPosY]}"/>
   <rotation name="rotGG02" unit="deg" x="0" y="180" z="90"  />
  </physvol>

  <physvol>
   <volumeref ref="volGroundGrid"/>
   <position name="posGroundGrid03" unit="cm" x="$GroundGridPosX" y="@{[$GroundGridPosY]}" z="@{[-$GroundGridPosY]}"/>
   <rotation name="rotGG03" unit="deg" x="180" y="0" z="90" />
  </physvol>

  <physvol>
   <volumeref ref="volGroundGrid"/>
   <position name="posGroundGrid04" unit="cm" x="$GroundGridPosX" y="@{[$GroundGridPosY]}" z="@{[$GroundGridPosY]}"/>
   <rotation name="rotGG04" unit="deg" x="180" y="180" z="90"/>
  </physvol>
EOF

  }

$CathodePosX=-$Argon_x/2 + 47.5 + 94.0 - 0.5*$CathodeSupportDistance-5.1;
$CathodePosY=0.5*($CathodeInnerStructureLength+2+$CathodeOuterRadious);

  if ( $Cathode_switch eq "on" )
  {
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volCathode"/>
   <position name="posCathode01" unit="cm" x="$CathodePosX" y="@{[-$CathodePosY]}" z="@{[$CathodePosY]}"/>
   <rotation name="rotCath01" unit="deg" x="0" y="0" z="90" />
  </physvol>

  <physvol>
   <volumeref ref="volCathode"/>
   <position name="posCathode02" unit="cm" x="$CathodePosX" y="@{[-$CathodePosY]}" z="@{[-$CathodePosY]}"/>
   <rotation name="rotCath02" unit="deg" x="270" y="0" z="90"  />
  </physvol>

  <physvol>
   <volumeref ref="volCathode"/>
   <position name="posCathode03" unit="cm" x="$CathodePosX" y="@{[$CathodePosY]}" z="@{[-$CathodePosY]}"/>
   <rotation name="rotCath03" unit="deg" x="180" y="0" z="90" />
  </physvol>

  <physvol>
   <volumeref ref="volCathode"/>
   <position name="posCathode04" unit="cm" x="$CathodePosX" y="@{[$CathodePosY]}" z="@{[$CathodePosY]}"/>
   <rotation name="rotCath04" unit="deg" x="90" y="0" z="90"/>
  </physvol>
EOF

  }

  if ( $FieldCage_switch eq "on" )
  {

$FFSPositionX = -$OriginXSet-1.3;
$FFSPositionY = 309.593;
$FFSPositionZ = 153.521;
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volFieldCageModule"/>
   <position name="posFieldCageModule01" unit="cm" x="$FFSPositionX" y="$FFSPositionY" z="$FFSPositionZ"/>
   <rotationref ref="rPlus90AboutY"/>
  </physvol>
EOF
$FFSPositionX = -$OriginXSet-1.3;
$FFSPositionY = 309.593;
$FFSPositionZ = -153.521;
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volFieldCageModule"/>
   <position name="posFieldCageModule02" unit="cm" x="$FFSPositionX" y="$FFSPositionY" z="$FFSPositionZ"/>
   <rotationref ref="rMinus90AboutY"/>
  </physvol>
EOF
$FFSPositionX = -$OriginXSet-1.3;
$FFSPositionY = -309.593;
$FFSPositionZ = -153.521;
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volFieldCageModule"/>
   <position name="posFieldCageModule05" unit="cm" x="$FFSPositionX" y="$FFSPositionY" z="$FFSPositionZ"/>
   <rotationref ref="rot05"/>

  </physvol>
EOF
$FFSPositionX = -$OriginXSet-1.3;
$FFSPositionY = -309.593;
$FFSPositionZ = 153.521;
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volFieldCageModule"/>
   <position name="posFieldCageModule06" unit="cm" x="$FFSPositionX" y="$FFSPositionY" z="$FFSPositionZ"/>
   <rotationref ref="rot06"/>

  </physvol>
EOF

$FFSPositionX = -$OriginXSet-1.3;
$FFSPositionY = 153.521;
$FFSPositionZ = -309.593;
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volFieldCageModule"/>
   <position name="posFieldCageModule03" unit="cm" x="$FFSPositionX" y="$FFSPositionY" z="$FFSPositionZ"/>
   <rotationref ref="rot03"/>

  </physvol>
EOF
$FFSPositionX = -$OriginXSet-1.3;
$FFSPositionY = -153.521;
$FFSPositionZ = -309.593;
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volFieldCageModule"/>
   <position name="posFieldCageModule04" unit="cm" x="$FFSPositionX" y="$FFSPositionY" z="$FFSPositionZ"/>
   <rotationref ref="rot04"/>

  </physvol>
EOF
$FFSPositionX = -$OriginXSet-1.3;
$FFSPositionY = -153.521;
$FFSPositionZ = 309.593;
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volFieldCageModule"/>
   <position name="posFieldCageModule07" unit="cm" x="$FFSPositionX" y="$FFSPositionY" z="$FFSPositionZ"/>
   <rotationref ref="rot07"/>

  </physvol>
EOF
$FFSPositionX = -$OriginXSet-1.3;
$FFSPositionY = 153.521;
$FFSPositionZ = 309.593;
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volFieldCageModule"/>
   <position name="posFieldCageModule08" unit="cm" x="$FFSPositionX" y="$FFSPositionY" z="$FFSPositionZ"/>
   <rotationref ref="rot08"/>

  </physvol>
EOF

  }


print CRYO <<EOF;
    </volume>
</structure>
</gdml>
EOF

close(CRYO);
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++ gen_Enclosure +++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Enclosure()
{

# Create the detector enclosure fragment file name,
# add file to list of output GDML fragments,
# and open it
    $ENCL = $basename."_DetEnclosure" . $suffix . ".gdml";
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

    <subtraction name="FoamPadding">
      <first ref="FoamPadBlock"/>
      <second ref="Cryostat"/>
      <positionref ref="posCenter"/>
    </subtraction>

    
    <box name="SteelSupportBlock" lunit="cm"
        x="@{[$Cryostat_x + 2*$FoamPadding + 2*$SteelSupport_x]}"
        y="@{[$Cryostat_y + 2*$FoamPadding + 2*$SteelSupport_y]}"
        z="@{[$Cryostat_z + 2*$FoamPadding + 2*$SteelSupport_z]}" />

    <subtraction name="SteelSupport">
        <first ref="SteelSupportBlock"/>
        <second ref="FoamPadBlock"/>
        <positionref ref="posCenter"/>
    </subtraction>
    
    
    <box name="DetEnclosure" lunit="cm" 
      x="$DetEncWidth"
      y="$DetEncHeight"
      z="$DetEncLength"/>

    <box name="boxCryoTop" x="1016.8" y="1016.8" z="61.8" lunit="cm"/>
    <box name="boxCryoWallLg" x="1140.4" y="1075.6" z="61.8" lunit="cm"/>
    <box name="boxCryoWallSm" x="1016.8" y="1075.6" z="61.8" lunit="cm"/>
    <box name="box1" x="160" y="160" z="61.8" lunit="cm"/>
    <box name="box2" x="158.2" y="158.2" z="56.2" lunit="cm"/>
    <box name="box3" x="137.2" y="137.2" z="61.801" lunit="cm"/>
    <box name="box4" x="158.2" y="13.6" z="27.4" lunit="cm"/>
    <box name="box5" x="158.2" y="6.425" z="24.96" lunit="cm"/>
    <box name="box11" x="108.4" y="160" z="61.8" lunit="cm"/>
    <box name="box12" x="107.5" y="158.2" z="56.2" lunit="cm"/>
    <box name="box13" x="97" y="137.2" z="61.81" lunit="cm"/>
    <box name="box16" x="107.5" y="13.6" z="27.4" lunit="cm"/>
    <box name="box17" x="107.5" y="6.425" z="24.96" lunit="cm"/>
    <box name="box21" x="137.8" y="160" z="61.8" lunit="cm"/>
    <box name="box22" x="136.9" y="158.2" z="56.2" lunit="cm"/>
    <box name="box23" x="126.4" y="137.2" z="61.801" lunit="cm"/>
    <box name="box26" x="102.5" y="13.6" z="27.4" lunit="cm"/>
    <box name="box27" x="102.5" y="6.425" z="24.96" lunit="cm"/>
    <box name="box31" x="170.2" y="160" z="61.8" lunit="cm"/>
    <box name="box32" x="169.3" y="158.2" z="56.2" lunit="cm"/>
    <box name="box33" x="158.8" y="137.2" z="61.801" lunit="cm"/>
    <box name="box36" x="135.8" y="13.6" z="27.4" lunit="cm"/>
    <box name="box37" x="135.8" y="6.425" z="24.96" lunit="cm"/>
    
    
    <subtraction name="boxHoll">
    <first ref="box1"/>
    <second ref="box2"/>
    </subtraction>
    <subtraction name="boxLargeCent">
    <first ref="boxHoll"/>
    <second ref="box3"/>
    </subtraction>
    <subtraction name="boxBarI">
    <first ref="box4"/>
    <second ref="box5"/>
    <position name="posBoxBarI" x="0" y="3.5876" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="boxBarCent">
    <first ref="boxBarI"/>
    <second ref="box5"/>
    <position name="posBoxBarCent" x="0" y="-3.5876" z="0" unit="cm"/>
    </subtraction>
    <union name="boxUniCent">
    <first ref="boxLargeCent"/>
    <second ref="boxBarCent"/>
    <position name="posBoxUniCent" x="0" y="0" z="-17.2" unit="cm"/>
    </union>
    <union name="UnitCent">
    <first ref="boxUniCent"/>
    <second ref="boxBarCent"/>
    <position name="posUnitCent" x="0" y="0" z="-17.2" unit="cm"/>
    <rotation name="rotUnitCent" x="0" y="0" z="90" unit="deg"/>
    </union>
    
    <subtraction name="boxTopHoll">
    <first ref="box11"/>
    <second ref="box12"/>
    <position name="posboxTopHoll" x="0.451" y="0" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="boxLargeTop">
    <first ref="boxTopHoll"/>
    <second ref="box13"/>
    <position name="posboxLargeTop" x="5.701" y="0" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="boxBarTopI">
    <first ref="box16"/>
    <second ref="box17"/>
    <position name="posboxBarTopI" x="0" y="3.5876" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="boxBarTop">
    <first ref="boxBarTopI"/>
    <second ref="box17"/>
    <position name="posboxBarTop" x="0" y="-3.5876" z="0" unit="cm"/>
    </subtraction>
    <union name="boxUniTop">
    <first ref="boxLargeTop"/>
    <second ref="boxBarCent"/>
    <position name="posboxUni1" x="5.6" y="0" z="-17.2" unit="cm"/>
    <rotation name="rotUni1" x="0" y="0" z="90" unit="deg"/>
    </union>
    <union name="UnitTop">
    <first ref="boxUniTop"/>
    <second ref="boxBarTop"/>
    <position name="posUniTop" x="0.45" y="0" z="-17.2" unit="cm"/>
    </union>
    
    <!--
    <subtraction name="boxCryoWallSmUS">
        <first ref="boxCryoWallSm"/>
        <second ref="BeamWindowStSu"/>
        <position name="posBWCryoWallUS" x="-34.5198845542345" y="131.897988017573" z="0." unit="cm"/>
        < rotationref ref="rBeamW3"/>
    </subtraction>
    -->
    
    <subtraction name="boxWallUHoll">
    <first ref="box21"/>
    <second ref="box22"/>
    <position name="posboxWallUHoll" x="0.451" y="0" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="boxLargeWallU">
    <first ref="boxWallUHoll"/>
    <second ref="box23"/>
    <position name="posboxLargeWallU" x="5.701" y="0" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="boxBarWallUI">
    <first ref="box26"/>
    <second ref="box27"/>
    <position name="posboxBarWallUI" x="0" y="3.5876" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="boxBarWallU">
    <first ref="boxBarWallUI"/>
    <second ref="box27"/>
    <position name="posboxBarWallU" x="0" y="-3.5876" z="0" unit="cm"/>
    </subtraction>
    <union name="boxUniWallU">
    <first ref="boxLargeWallU"/>
    <second ref="boxBarCent"/>
    <position name="posboxUni2" x="-9.1" y="0" z="-17.2" unit="cm"/>
    <rotation name="rotUni2" x="0" y="0" z="90" unit="deg"/>
    </union>
    <union name="UnitWallU">
    <first ref="boxUniWallU"/>
    <second ref="boxBarWallU"/>
    <position name="posUniWallU" x="-16.75" y="0" z="-17.2" unit="cm"/>
    </union>
    
    <subtraction name="boxWallLHoll">
    <first ref="box31"/>
    <second ref="box32"/>
    <position name="posboxWallLHoll" x="0.451" y="0" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="boxLargeWallL">
    <first ref="boxWallLHoll"/>
    <second ref="box33"/>
    <position name="posboxLargeWallL" x="5.701" y="0" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="boxBarWallLI">
    <first ref="box36"/>
    <second ref="box37"/>
    <position name="posboxBarWallLI" x="0" y="3.5876" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="boxBarWallL">
    <first ref="boxBarWallLI"/>
    <second ref="box37"/>
    <position name="posboxBarWallL" x="0" y="-3.5876" z="0" unit="cm"/>
    </subtraction>
    <union name="boxUniWallL">
    <first ref="boxLargeWallL"/>
    <second ref="boxBarCent"/>
    <position name="posboxUni3" x="-25.3" y="0" z="-17.2" unit="cm"/>
    <rotation name="rotUni3" x="0" y="0" z="90" unit="deg"/>
    </union>
    <union name="UnitWallL">
    <first ref="boxUniWallL"/>
    <second ref="boxBarWallL"/>
    <position name="posUniWallL" x="-16.3" y="0" z="-17.2" unit="cm"/>
    </union>

    
</solids>
EOF


# Detector enclosure structure
    print ENCL <<EOF;
<structure>
    <volume name="volFoamPadding">
      <materialref ref="foam_protoDUNEdp"/>
      <solidref ref="FoamPadding"/>
    </volume>
    
    <volume name="volSteelSupport">
        <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
        <solidref ref="SteelSupport"/>
    </volume>

    <volume name="volUnitCent">
    <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
    <solidref ref="UnitCent"/>
    </volume>
    <volume name="volUnitTop">
    <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
    <solidref ref="UnitTop"/>
    </volume>
    <volume name="volUnitWallS">
    <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
    <solidref ref="UnitWallU"/>
    </volume>
    <volume name="volUnitWallL">
    <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
    <solidref ref="UnitWallL"/>
    </volume>
    
    <volume name="volSteelSupport_TB">
        <materialref ref="Air"/>
        <solidref ref="boxCryoTop"/>
        <physvol name="volUnitTBCent_0-0">
            <volumeref ref="volUnitCent"/>
            <position name="posUnitTBCent_0-0" x="-320" y="-320" z="0" unit="cm"/>
        </physvol>
        <physvol name="volUnitTBCent_0-1">
            <volumeref ref="volUnitCent"/>
            <position name="posUnitTBCent_0-1" x="-320" y="-160" z="0" unit="cm"/>
        </physvol>
    <physvol name="volUnitTBCent_0-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_0-2" x="-320" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_0-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_0-3" x="-320" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_0-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_0-4" x="-320" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBE_0">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBE_0" x="454.2" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBS_0">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBS_0" x="-320" y="454.2" z="0" unit="cm"/>
    <rotation name="rotUnitTBS_0" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBW_0">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBW_0" x="-454.2" y="-320" z="0" unit="cm"/>
    <rotation name="rotUnitTBW_0" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBN_0">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBN_0" x="-320" y="-454.2" z="0" unit="cm"/>
    <rotation name="rotUnitTBN_0" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBCent_1-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_1-0" x="-160" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_1-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_1-1" x="-160" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_1-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_1-2" x="-160" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_1-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_1-3" x="-160" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_1-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_1-4" x="-160" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBE_1">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBE_1" x="454.2" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBS_1">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBS_1" x="-160" y="454.2" z="0" unit="cm"/>
    <rotation name="rotUnitTBS_1" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBW_1">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBW_1" x="-454.2" y="-160" z="0" unit="cm"/>
    <rotation name="rotUnitTBW_1" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBN_1">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBN_1" x="-160" y="-454.2" z="0" unit="cm"/>
    <rotation name="rotUnitTBN_1" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBCent_2-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_2-0" x="0" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_2-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_2-1" x="0" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_2-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_2-2" x="0" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_2-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_2-3" x="0" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_2-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_2-4" x="0" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBE_2">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBE_2" x="454.2" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBS_2">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBS_2" x="0" y="454.2" z="0" unit="cm"/>
    <rotation name="rotUnitTBS_2" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBW_2">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBW_2" x="-454.2" y="0" z="0" unit="cm"/>
    <rotation name="rotUnitTBW_2" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBN_2">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBN_2" x="0" y="-454.2" z="0" unit="cm"/>
    <rotation name="rotUnitTBN_2" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBCent_3-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_3-0" x="160" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_3-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_3-1" x="160" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_3-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_3-2" x="160" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_3-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_3-3" x="160" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_3-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_3-4" x="160" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBE_3">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBE_3" x="454.2" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBS_3">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBS_3" x="160" y="454.2" z="0" unit="cm"/>
    <rotation name="rotUnitTBS_3" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBW_3">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBW_3" x="-454.2" y="160" z="0" unit="cm"/>
    <rotation name="rotUnitTBW_3" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBN_3">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBN_3" x="160" y="-454.2" z="0" unit="cm"/>
    <rotation name="rotUnitTBN_3" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBCent_4-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_4-0" x="320" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_4-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_4-1" x="320" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_4-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_4-2" x="320" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_4-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_4-3" x="320" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBCent_4-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitTBCent_4-4" x="320" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBE_4">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBE_4" x="454.2" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitTBS_4">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBS_4" x="320" y="454.2" z="0" unit="cm"/>
    <rotation name="rotUnitTBS_4" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBW_4">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBW_4" x="-454.2" y="320" z="0" unit="cm"/>
    <rotation name="rotUnitTBW_4" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitTBN_4">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitTBN_4" x="320" y="-454.2" z="0" unit="cm"/>
    <rotation name="rotUnitTBN_4" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    </volume>
    
    
    <volume name="volSteelSupport_US">
    <materialref ref="Air"/>
    <solidref ref="boxCryoWallSm"/>
    
    <physvol name="volUnitUSCent_0-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_0-0" x="-320" y="-320" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_0-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_0-1" x="-320" y="-160" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_0-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_0-2" x="-320" y="0" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_0-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_0-3" x="-320" y="160" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_0-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_0-4" x="-320" y="320" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    <physvol name="volUnitUSE_0">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitUSE_0" x="454.2" y="-320" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutX"/>
    </physvol>
    <physvol name="volUnitUSS_0">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitUSS_0" x="-320" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitUSS_0" x="0" y="180" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitUSW_0">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitUSW_0" x="-454.2" y="-320" z="0" unit="cm"/>
    <rotation  name="rotUnitUSW_0" x="180" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitUSN_0">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitUSN_0" x="-320" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitUSN_0" x="0" y="180" z="-270" unit="deg"/>
    </physvol>
    
    
    <physvol name="volUnitUSCent_1-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_1-0" x="-160" y="-320" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_1-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_1-1" x="-160" y="-160" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_1-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_1-2" x="-160" y="0" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_1-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_1-3" x="-160" y="160" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_1-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_1-4" x="-160" y="320" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    <physvol name="volUnitUSE_1">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitUSE_1" x="454.2" y="-160" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutX"/>
    </physvol>
    <physvol name="volUnitUSS_1">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitUSS_1" x="-160" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitUSS_1" x="0" y="180" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitUSW_1">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitUSW_1" x="-454.2" y="-160" z="0" unit="cm"/>
    <rotation  name="rotUnitUSW_1" x="180" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitUSN_1">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitUSN_1" x="-160" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitUSN_1" x="0" y="180" z="-270" unit="deg"/>
    </physvol>
    
    
    <physvol name="volUnitUSCent_2-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_2-0" x="0" y="-320" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_2-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_2-1" x="0" y="-160" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_2-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_2-2" x="0" y="0" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_2-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_2-3" x="0" y="160" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_2-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_2-4" x="0" y="320" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    <physvol name="volUnitUSE_2">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitUSE_2" x="454.2" y="0" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutX"/>
    </physvol>
    <physvol name="volUnitUSS_2">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitUSS_2" x="0" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitUSS_2" x="0" y="180" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitUSW_2">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitUSW_2" x="-454.2" y="0" z="0" unit="cm"/>
    <rotation  name="rotUnitUSW_2" x="180" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitUSN_2">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitUSN_2" x="0" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitUSN_2" x="0" y="180" z="-270" unit="deg"/>
    </physvol>
    
    
    <physvol name="volUnitUSCent_3-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_3-0" x="160" y="-320" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_3-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_3-1" x="160" y="-160" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_3-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_3-2" x="160" y="0" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_3-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_3-3" x="160" y="160" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_3-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_3-4" x="160" y="320" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    <physvol name="volUnitUSE_3">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitUSE_3" x="454.2" y="160" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutX"/>
    </physvol>
    <physvol name="volUnitUSS_3">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitUSS_3" x="160" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitUSS_3" x="0" y="180" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitUSW_3">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitUSW_3" x="-454.2" y="160" z="0" unit="cm"/>
    <rotation  name="rotUnitUSW_3" x="180" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitUSN_3">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitUSN_3" x="160" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitUSN_3" x="0" y="180" z="-270" unit="deg"/>
    </physvol>
    
    
    <physvol name="volUnitUSCent_4-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_4-0" x="320" y="-320" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_4-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_4-1" x="320" y="-160" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_4-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_4-2" x="320" y="0" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_4-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_4-3" x="320" y="160" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    
    <physvol name="volUnitUSCent_4-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitUSCent_4-4" x="320" y="320" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutY"/>
    </physvol>
    <physvol name="volUnitUSE_4">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitUSE_4" x="454.2" y="320" z="0" unit="cm"/>
    <rotationref ref="rPlus180AboutX"/>
    </physvol>
    <physvol name="volUnitUSS_4">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitUSS_4" x="320" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitUSS_4" x="0" y="180" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitUSW_4">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitUSW_4" x="-454.2" y="320" z="0" unit="cm"/>
    <rotation  name="rotUnitUSW_4" x="180" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitUSN_4">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitUSN_4" x="320" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitUSN_4" x="0" y="180" z="-270" unit="deg"/>
    </physvol>
    </volume>
    
    
    
    <volume name="volSteelSupport_WS">
    <materialref ref="Air"/>
    <solidref ref="boxCryoWallSm"/>
    <physvol name="volUnitWSCent_0-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_0-0" x="-320" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_0-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_0-1" x="-320" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_0-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_0-2" x="-320" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_0-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_0-3" x="-320" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_0-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_0-4" x="-320" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSE_0">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitWSE_0" x="454.2" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSS_0">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitWSS_0" x="-320" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitWSS_0" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSW_0">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitWSW_0" x="-454.2" y="-320" z="0" unit="cm"/>
    <rotation  name="rotUnitWSW_0" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSN_0">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitWSN_0" x="-320" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitWSN_0" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSCent_1-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_1-0" x="-160" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_1-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_1-1" x="-160" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_1-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_1-2" x="-160" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_1-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_1-3" x="-160" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_1-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_1-4" x="-160" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSE_1">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitWSE_1" x="454.2" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSS_1">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitWSS_1" x="-160" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitWSS_1" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSW_1">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitWSW_1" x="-454.2" y="-160" z="0" unit="cm"/>
    <rotation  name="rotUnitWSW_1" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSN_1">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitWSN_1" x="-160" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitWSN_1" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSCent_2-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_2-0" x="0" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_2-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_2-1" x="0" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_2-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_2-2" x="0" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_2-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_2-3" x="0" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_2-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_2-4" x="0" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSE_2">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitWSE_2" x="454.2" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSS_2">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitWSS_2" x="0" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitWSS_2" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSW_2">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitWSW_2" x="-454.2" y="0" z="0" unit="cm"/>
    <rotation  name="rotUnitWSW_2" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSN_2">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitWSN_2" x="0" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitWSN_2" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSCent_3-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_3-0" x="160" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_3-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_3-1" x="160" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_3-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_3-2" x="160" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_3-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_3-3" x="160" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_3-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_3-4" x="160" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSE_3">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitWSE_3" x="454.2" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSS_3">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitWSS_3" x="160" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitWSS_3" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSW_3">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitWSW_3" x="-454.2" y="160" z="0" unit="cm"/>
    <rotation  name="rotUnitWSW_3" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSN_3">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitWSN_3" x="160" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitWSN_3" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSCent_4-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_4-0" x="320" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_4-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_4-1" x="320" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_4-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_4-2" x="320" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_4-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_4-3" x="320" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSCent_4-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitWSCent_4-4" x="320" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSE_4">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitWSE_4" x="454.2" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitWSS_4">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitWSS_4" x="320" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitWSS_4" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSW_4">
    <volumeref ref="volUnitTop"/>
    <position name="posUnitWSW_4" x="-454.2" y="320" z="0" unit="cm"/>
    <rotation  name="rotUnitWSW_4" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitWSN_4">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitWSN_4" x="320" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitWSN_4" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    </volume>
    

    <volume name="volSteelSupport_LR">
    <materialref ref="Air"/>
    <solidref ref="boxCryoWallLg"/>
    <physvol name="volUnitLRCent_0-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_0-0" x="-320" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_0-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_0-1" x="-320" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_0-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_0-2" x="-320" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_0-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_0-3" x="-320" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_0-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_0-4" x="-320" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRE_0">
    <volumeref ref="volUnitWallL"/>
    <position name="posUnitLRE_0" x="485.1" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRS_0">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitLRS_0" x="-320" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitLRS_0" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRW_0">
    <volumeref ref="volUnitWallL"/>
    <position name="posUnitLRW_0" x="-485.1" y="-320" z="0" unit="cm"/>
    <rotation  name="rotUnitLRW_0" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRN_0">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitLRN_0" x="-320" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitLRN_0" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRCent_1-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_1-0" x="-160" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_1-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_1-1" x="-160" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_1-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_1-2" x="-160" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_1-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_1-3" x="-160" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_1-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_1-4" x="-160" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRE_1">
    <volumeref ref="volUnitWallL"/>
    <position name="posUnitLRE_1" x="485.1" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRS_1">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitLRS_1" x="-160" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitLRS_1" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRW_1">
    <volumeref ref="volUnitWallL"/>
    <position name="posUnitLRW_1" x="-485.1" y="-160" z="0" unit="cm"/>
    <rotation  name="rotUnitLRW_1" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRN_1">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitLRN_1" x="-160" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitLRN_1" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRCent_2-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_2-0" x="0" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_2-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_2-1" x="0" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_2-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_2-2" x="0" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_2-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_2-3" x="0" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_2-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_2-4" x="0" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRE_2">
    <volumeref ref="volUnitWallL"/>
    <position name="posUnitLRE_2" x="485.1" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRS_2">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitLRS_2" x="0" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitLRS_2" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRW_2">
    <volumeref ref="volUnitWallL"/>
    <position name="posUnitLRW_2" x="-485.1" y="0" z="0" unit="cm"/>
    <rotation  name="rotUnitLRW_2" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRN_2">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitLRN_2" x="0" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitLRN_2" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRCent_3-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_3-0" x="160" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_3-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_3-1" x="160" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_3-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_3-2" x="160" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_3-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_3-3" x="160" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_3-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_3-4" x="160" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRE_3">
    <volumeref ref="volUnitWallL"/>
    <position name="posUnitLRE_3" x="485.1" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRS_3">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitLRS_3" x="160" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitLRS_3" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRW_3">
    <volumeref ref="volUnitWallL"/>
    <position name="posUnitLRW_3" x="-485.1" y="160" z="0" unit="cm"/>
    <rotation  name="rotUnitLRW_3" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRN_3">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitLRN_3" x="160" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitLRN_3" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRCent_4-0">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_4-0" x="320" y="-320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_4-1">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_4-1" x="320" y="-160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_4-2">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_4-2" x="320" y="0" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_4-3">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_4-3" x="320" y="160" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRCent_4-4">
    <volumeref ref="volUnitCent"/>
    <position name="posUnitLRCent_4-4" x="320" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRE_4">
    <volumeref ref="volUnitWallL"/>
    <position name="posUnitLRE_4" x="485.1" y="320" z="0" unit="cm"/>
    </physvol>
    <physvol name="volUnitLRS_4">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitLRS_4" x="320" y="468.9" z="0" unit="cm"/>
    <rotation name="rotUnitLRS_4" x="0" y="0" z="-90" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRW_4">
    <volumeref ref="volUnitWallL"/>
    <position name="posUnitLRW_4" x="-485.1" y="320" z="0" unit="cm"/>
    <rotation  name="rotUnitLRW_4" x="0" y="0" z="-180" unit="deg"/>
    </physvol>
    <physvol name="volUnitLRN_4">
    <volumeref ref="volUnitWallS"/>
    <position name="posUnitLRN_4" x="320" y="-468.9" z="0" unit="cm"/>
    <rotation  name="rotUnitLRN_4" x="0" y="0" z="-270" unit="deg"/>
    </physvol>
    </volume>
    
    <volume name="volDetEnclosure">
      <materialref ref="Air"/>
      <solidref ref="DetEnclosure"/>

       <physvol>
           <volumeref ref="volFoamPadding"/>
           <positionref ref="posCryoInDetEnc"/>
       </physvol>
    
    
    <physvol name="volSteelSupport_Top">
        <volumeref ref="volSteelSupport_TB"/>
        <position name="posSteelSupport_Top" x="$posTopSteelStruct+31.1" y="$posCryoInDetEnc_y" z="0" unit="cm"/>
        <rotation name="rotSteelSupport_Top" x="0" y="-90" z="0" unit="deg"/>
    
    </physvol>
    
    <physvol name="volSteelSupport_Bottom">
        <volumeref ref="volSteelSupport_TB"/>
        <position name="posSteelSupport_Bottom" x="$posBotSteelStruct-31.1" y="$posCryoInDetEnc_y" z="0" unit="cm"/>
        <rotation name="rotSteelSupport_Bottom" x="0" y="90" z="0" unit="deg"/>
    </physvol>
    
    <physvol>
        <volumeref ref="volSteelSupport_US"/>
        <position name="posSteelSupport_US" x="0" y="$posCryoInDetEnc_y" z="$posZFrontSteelStruct-31.1" unit="cm"/>
        <rotation name="rotSteelSupport_Front" x="0" y="0" z="90" unit="deg"/>
    </physvol>
    
    <physvol name="volSteelSupport_DS">
        <volumeref ref="volSteelSupport_WS"/>
        <position name="posSteelSupport_DS" x="0" y="$posCryoInDetEnc_y" z="$posZBackSteelStruct+31.1" unit="cm"/>
        <rotation name="rotSteelSupport_Back" x="0" y="0" z="90" unit="deg"/>
    </physvol>
    
    
    <physvol name="volSteelSupport_LS">
        <volumeref ref="volSteelSupport_LR"/>
        <position name="posSteelSupport_LS" x="0" y="$posLeftSteelStruct+$posCryoInDetEnc_y-31.1" z="0" unit="cm"/>
        <rotation name="rotSteelSupport_LS" x="-90" y="0" z="90" unit="deg"/>
    </physvol>
    
    <physvol name="volSteelSupport_RS">
        <volumeref ref="volSteelSupport_LR"/>
        <position name="posSteelSupport_RS" x="0" y="$posRightSteelStruct+$posCryoInDetEnc_y+31.1" z="0" unit="cm"/>
        <rotation name="rotSteelSupport_RS" x="90" y="0" z="90" unit="deg"/>
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
    $WORLD = $basename."_World" . $suffix . ".gdml";
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


print "Some key parameters for dual-phase LAr TPC (unit cm unless noted otherwise)\n";
print "CRM active area    : $widthCRM_active x $lengthCRM_active\n";
print "CRM total area     : $widthCRM x $lengthCRM\n";
print "TPC active volume  : $driftTPCActive x $widthTPCActive x $lengthTPCActive\n";
print "Argon buffer       : ($xLArBuffer, $yLArBuffer, $zLArBuffer) \n"; 
print "Detector enclosure : $DetEncWidth x $DetEncHeight x $DetEncLength\n";
print "TPC Origin         : ($OriginXSet, $OriginYSet, $OriginZSet) \n";
print "Field Cage         : $FieldCage_switch \n";
print "Cathode            : $Cathode_switch \n";;
print "GroundGrid         : $GroundGrid_switch \n";
print "ExtractionGrid     : $ExtractionGrid_switch \n";
print "LEMs               : $LEMs_switch \n";
print "PMTs               : $pmt_switch \n";
# run the sub routines that generate the fragments




gen_Define(); 	 # generates definitions at beginning of GDML
gen_Materials(); # generates materials to be used

if ( $pmt_switch eq "on" ) {  gen_pmt();	}
if ( $FieldCage_switch eq "on" ) {  gen_FieldCage();	}
if ( $GroundGrid_switch eq "on" ) {  gen_GroundGrid();	}
if ( $Cathode_switch eq "on" ) {  gen_Cathode();	}
if ( $ExtractionGrid_switch eq "on" ) {  gen_ExtractionGrid();	}
if ( $LEMs_switch eq "on" ) {  gen_LEMs();	}

gen_TPC();       # generate TPC for a given unit CRM
gen_Cryostat();  # 
gen_Enclosure(); # 
gen_World();	 # places the enclosure among DUSEL Rockk
write_fragments(); # writes the XML input for make_gdml.pl
		   # which zips together the final GDML
exit;
