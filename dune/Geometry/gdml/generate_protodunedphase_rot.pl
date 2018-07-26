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
	    "wires|w:s" => \$wires);

my $pmt_switch="on";

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

# set wires on to be the default, unless given an input by the user
$wires_on = 1; # 1=on, 0=off
if (defined $wires)
{
    $wires_on = $wires
}

$tpc_on = 1;

$basename = "protodunedphase_rot";
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
$widthCRM_active  = $wirePitch * $nChannelsViewPerCRM; # 300
$lengthCRM_active = $wirePitch * $nChannelsViewPerCRM; # 300

$widthCRM  = $widthCRM_active + 2 * $borderCRM; # 301
$lengthCRM = $lengthCRM_active + 2 * $borderCRM; # 301

# number of CRMs in y and z
$nCRM_x   = 2;
$nCRM_z   = 2;

# calculate tpc area based on number of CRMs and their dimensions
$widthTPCActive  = $nCRM_x * $widthCRM;  # 602
$lengthTPCActive = $nCRM_z * $lengthCRM; # 602

# active volume dimensions 
$driftTPCActive  = 597.0;

# model anode strips as wires
$padWidth  = 0.015;
$ReadoutPlane = 2 * $padWidth; # 0.03

#$padHeight = 0.0035; 

##################################################################
############## Parameters for TPC and inner volume ###############

# inner volume dimensions of the cryostat
$Argon_x = 854.8;
$Argon_y = 790.0;
$Argon_z = 854.8;

# width of gas argon layer on top
$HeightGaseousAr = 51.5;

# size of liquid argon buffer
$xLArBuffer = 0.5 * ($Argon_x - $widthTPCActive); # 126.2
$yLArBuffer = $Argon_y - $driftTPCActive - $HeightGaseousAr - $ReadoutPlane; #107.77
$zLArBuffer = 0.5 * ($Argon_z - $lengthTPCActive); # 126.2

# cryostat 
$SteelThickness = 0.2; # membrane

$Cryostat_x = $Argon_x + 2*$SteelThickness; # 854.64
$Cryostat_y = $Argon_y + 2*$SteelThickness; # 789.84
$Cryostat_z = $Argon_z + 2*$SteelThickness; # 854.64


##################################################################
############## Parameters for PMTs ###############

#pos in cm inside the cryostat
 #$pmt_pos_y  =  -300-(5.5 * 2.54)+5.115; just below the active volume
 #$pmt_pos_y  =  -300-(5.5 * 2.54)+5.115 + 22 -141.5; just on the top of the cryostat (new version)
# $pmt_pos_y  =  -300+5.115 - 1.27*2.54 + 22 -107.8;  #8.8cm from the bottom volTPCActive
 $HeightPMT = 37.0;

#pmts not equally spaced:
 @pmt_pos = ( ' z="-238" x="-170" ', #pmt1
		' z="-238" x="-34" ', #pmt2
		' z="-238" x="34" ', #pmt3
		' z="-238" x="170" ', #pmt4
		' z="-170" x="-238" ', #pmt5
		' z="-170" x="-102" ', #pmt6
		' z="-170" x="102" ', #pmt7
		' z="-170" x="238" ', #pmt8
		' z="-102" x="-170" ', #pmt9
		' z="-102" x="-34" ', #pmt10
		' z="-102" x="34" ', #pmt11
		' z="-102" x="170" ', #pmt12
		' z="-34" x="-238" ', #pmt13
		' z="-34" x="-102" ', #pmt14
		' z="-34" x="-34" ', #pmt15
		' z="-34" x="34" ', #pmt16
		' z="-34" x="102" ', #pmt17
		' z="-34" x="238" ', #pmt18
		' z="34" x="-238" ', #pmt19
		' z="34" x="-102" ', #pmt20
		' z="34" x="-34" ', #pmt21
		' z="34" x="34" ', #pmt22
		' z="34" x="102" ', #pmt23
		' z="34" x="238" ', #pmt24
		' z="102" x="-170" ', #pmt25
		' z="102" x="-34" ', #pmt26
		' z="102" x="34" ', #pmt27
		' z="102" x="170" ', #pmt28
		' z="170" x="-238" ', #pmt29
		' z="170" x="-102" ', #pmt30
		' z="170" x="102" ', #pmt31
		' z="170" x="238" ', #pmt32
		' z="238" x="-170" ', #pmt33
		' z="238" x="-34" ', #pmt34
		' z="238" x="34" ', #pmt35
		' z="238" x="170" '); #PMT36

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

$posTopSteelStruct = $Argon_y/2+$FoamPadding+$SteelSupport_y;
$posBotSteelStruct = -$Argon_y/2-$FoamPadding-$SteelSupport_y;
$posZBackSteelStruct = $Argon_z/2+$FoamPadding+$SteelSupport_z;
$posZFrontSteelStruct = -$Argon_z/2-$FoamPadding-$SteelSupport_z;
$posLeftSteelStruct = $Argon_x/2+$FoamPadding+$SteelSupport_x;
$posRightSteelStruct = -$Argon_x/2-$FoamPadding-$SteelSupport_x;

# 2*AirThickness is added to the world volume in x, y and z
$AirThickness = 3000;

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
	      - $SpaceSteelSupportToCeiling/2.0
              - $SteelSupport_y
              - $FoamPadding
              - $SteelThickness
              - $yLArBuffer
              - $driftTPCActive/2.0;
            

$OriginXSet =  $DetEncWidth/2.0
              -$SpaceSteelSupportToWall
              -$SteelSupport_x
              -$FoamPadding
              -$SteelThickness
              -$xLArBuffer
              -$widthTPCActive/2.0;


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

   <position name="posCryoInDetEnc"     unit="cm" x="$posCryoInDetEnc_x" y="0" z="0"/>
   <position name="posCenter"           unit="cm" x="0" y="0" z="0"/>
   <rotation name="rPlus90AboutX"       unit="deg" x="90" y="0" z="0"/>
   <rotation name="rPlus90AboutY"       unit="deg" x="90" y="90" z="0"/>
   <rotation name="rMinus90AboutY"      unit="deg" x="0" y="270" z="0"/>
   <rotation name="rPlus90AboutX"      unit="deg" x="90" y="0" z="0"/>
   <rotation name="rMinus90AboutYMinus90AboutX"       unit="deg" x="270" y="270" z="0"/>
   <rotation name="rPlus180AboutX"	unit="deg" x="180" y="0"   z="0"/>
   <rotation name="rPlus180AboutY"	unit="deg" x="0" y="180"   z="0"/>
   <rotation name="rPlus180AboutXPlus180AboutY"	unit="deg" x="180" y="180"   z="0"/>
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
#++++++++++++++++++++++++++++++++++++++++ gen_TPC ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub gen_TPC()
{
    # CRM active volume
    my $TPCActive_y = $driftTPCActive;
    my $TPCActive_x = $widthCRM_active;
    my $TPCActive_z = $lengthCRM_active;

    # CRM total volume
    my $TPC_y = $TPCActive_y + $ReadoutPlane;
    my $TPC_x = $widthCRM;
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
     x="$TPCActive_x" 
     y="$padWidth"  
     z="$TPCActive_z"/>
   <box name="CRMZPlane" lunit="cm" 
      x="$TPCActive_x"
      y="$padWidth"
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
     z="$TPCActive_x"               
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
my $xpos = -0.5 * $TPCActive_x + ($i+0.5)*$wirePitch + 0.5*$padWidth;

print TPC <<EOF;
    <physvol>
      <volumeref ref="volTPCWireV"/> 
      <position name="posWireV$i" unit="cm" x="$xpos" y="0" z="0"/>
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
     <rotationref ref="rPlus90AboutY"/>
    </physvol>
EOF
}
}


print TPC <<EOF;
   </volume>
EOF


$posVplane[0] = 0;
$posVplane[1] = 0.5*$TPC_y - 1.5*$padWidth;
$posVplane[2] = 0;

$posZplane[0] = 0; 
$posZplane[1] = 0.5*$TPC_y - 0.5*$padWidth;
$posZplane[2] = 0;

$posTPCActive[0] = 0;
$posTPCActive[1] = -$ReadoutPlane;
$posTPCActive[2] = 0;


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
   </volume>
EOF


print TPC <<EOF;
</structure>
</gdml>
EOF

close(TPC);
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


	print PMT <<EOF;


<solids>
 <tube name="PMTVolume"
  rmax="@{[(6.5*2.54)]}"
  z="@{[(11.1*2.54)]}"
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

<structure>
 <volume name="pmtCoatVol">
  <materialref ref="TPB"/>
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
   <position name="posallpmt" unit="cm" x="0" y="0" z="@{[1.27*2.54]}"/>
  </physvol>

 <physvol name="volOpDetSensitive">
  <volumeref ref="pmtCoatVol"/>
  <position name="posOpDetSensitive" unit="cm" x="0" y="0" z="@{[1.27*2.54 - (2.23*2.54)]}"/>
  </physvol>

 </volume>

</structure>

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
      x="$Argon_x"
      y="$HeightGaseousAr"
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
    </volume>

    <volume name="volCryostat">
      <materialref ref="LAr" />
      <solidref ref="Cryostat" />
      <physvol>
        <volumeref ref="volGaseousArgon"/>
        <position name="posGaseousArgon" unit="cm" x="0" y="$Argon_y/2-$HeightGaseousAr/2"  z="0"/>
         </physvol>
      <physvol>
        <volumeref ref="volSteelShell"/>
        <position name="posSteelShell" unit="cm" x="0" y="0" z="0"/>
      </physvol>
EOF


if ($tpc_on==1) # place TPC inside croysotat
{

$posY =  $Argon_y/2 - $HeightGaseousAr - 0.5*($driftTPCActive + $ReadoutPlane); 
for($ii=0;$ii<$nCRM_z;$ii++)
{
    $posZ = -0.5*$Argon_z + $zLArBuffer + ($ii+0.5)*$lengthCRM;

    for($jj=0;$jj<$nCRM_x;$jj++)
    {
	$posX = -0.5*$Argon_x + $xLArBuffer + ($jj+0.5)*$widthCRM;
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
    $pmt_pos_y =  -$Argon_y/2 + 0.5*($HeightPMT);
    for ( $i=0; $i<36; $i=$i+1 ) { # pmts with coating
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volPMT_coated"/>
   
   <position name="posPMT$i" unit="cm" y="$pmt_pos_y" @pmt_pos[$i]/>
   <rotationref ref="rPlus90AboutX"/>
  </physvol>
EOF
    }

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
      <second ref="FoamPadding"/>
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
    <position name="posSteelSupport_Top" x="0" y="$posTopSteelStruct+30.9" z="0" unit="cm"/>
    <rotation name="rotSteelSupport_Top" x="90" y="0" z="0" unit="deg"/>
    
    </physvol>
    
    
    <physvol name="volSteelSupport_Bottom">
    <volumeref ref="volSteelSupport_TB"/>
    <position name="posSteelSupport_Bottom" x="0" y="$posBotSteelStruct-30.9" z="0" unit="cm"/>
    <rotation name="rotSteelSupport_Bottom" x="-90" y="0" z="0" unit="deg"/>
    </physvol>
    
    
    <physvol>
    <volumeref ref="volSteelSupport_US"/>
    <position name="posSteelSupport_US" x="0" y="0" z="$posZFrontSteelStruct-30.9" unit="cm"/>
    <rotation name="rotSteelSupport_Front" x="0" y="0" z="0" unit="deg"/>
    </physvol>
    
    
    <physvol name="volSteelSupport_DS">
    <volumeref ref="volSteelSupport_WS"/>
    <position name="posSteelSupport_DS" x="0" y="0" z="$posZBackSteelStruct+30.9" unit="cm"/>
    <rotation name="rotSteelSupport_Back" x="0" y="0" z="" unit="deg"/>
    </physvol>
    
    
    <physvol name="volSteelSupport_LS">
    <volumeref ref="volSteelSupport_LR"/>
    <position name="posSteelSupport_LS" x="$posLeftSteelStruct+30.9" y="0" z="0" unit="cm"/>
    <rotation name="rotSteelSupport_LS" x="0" y="-90" z="0" unit="deg"/>
    </physvol>
    
    
    <physvol name="volSteelSupport_RS">
    <volumeref ref="volSteelSupport_LR"/>
    <position name="posSteelSupport_RS" x="$posRightSteelStruct-30.9" y="0" z="0" unit="cm"/>
    <rotation name="rotSteelSupport_RS" x="0" y="90" z="0" unit="deg"/>
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
      x="@{[$DetEncWidth+2*$AirThickness]}" 
      y="@{[$DetEncHeight+2*$AirThickness]}" 
      z="@{[$DetEncLength+2*$AirThickness]}"/>
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
print "PMTs               : $pmt_switch \n";

# run the sub routines that generate the fragments
if ( $pmt_switch eq "on" ) {  gen_pmt();	}

gen_Define(); 	 # generates definitions at beginning of GDML
gen_Materials(); # generates materials to be used
gen_TPC();       # generate TPC for a given unit CRM
gen_Cryostat();  # 
gen_Enclosure(); # 
gen_World();	 # places the enclosure among DUSEL Rock


write_fragments(); # writes the XML input for make_gdml.pl
		   # which zips together the final GDML
exit;
