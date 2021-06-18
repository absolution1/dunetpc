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
use gdmlMaterials;
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
            "workspace|k:s" => \$workspace,
            "simpleStSu|u:s" => \$simpleStSu); 

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
        $simpleStSu = 1;
        print "\t\tCreating full geometry.\n";
    }
    elsif ( $protoDune == 1 )
    {
        print "\t\tCreating rough version of protoDUNE.\n";
        if ( ! defined $simpleStSu )
        {
           $simpleStSu = 0;
        }
        elsif ( $simpleStSu == 1 )      
        {
           print "\t\tUsing simple steel support geometry for protoDUNE.\n"; 
        }
    }
    if($protoDune == 0 && $simpleStSu == 0)
    {
       print "\t\tThis version of the steel support is only compatible with protoDUNE.\n";
       exit;
    }
}
elsif ( $workspace == 1)
{
    print "\t\tCreating smaller workspace geometry.\n";
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
$UWirePitch             =   0.4669;
$VWirePitch             =   0.4669;
$XWirePitch             =   0.4792;
$UAngle                 =   35.707;
$VAngle                 =   35.707;
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

$nAPAWide               =     3; 
$nAPAHigh               =     2;
$nAPALong               =     25;

if($protoDune==1){
    $nAPAWide           =     2; 
    $nAPAHigh           =     1;
    $nAPALong           =     3;    
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
    $nAPAWide           =     1; 
    $nAPALong           =     2;
}

$nAPAs                 =     $nAPAWide*$nAPAHigh*$nAPALong;


#$G10thickness = $inch/8;
$G10thickness = 0.335;  #1/8 inch plus epoxy
$WrapCover    = $inch/16;

#$SpaceAPAToCryoWall    =     15;
#$SpaceAPAToFloor       =     59.4  - 10.2; # 10.2cm just adjusts for how spaces are reported
#$SpaceAPAToTopLAr      =     50.9  - 10.2;
#$UpstreamLArPadding    =     311.4 - 10.2;
#$DownstreamLArPadding  =     59.4  - 10.2;

$Argon_x                =       854.8; #inside dimension of the cryostat
$Argon_y                =       790.0;
$Argon_z                =       854.8;

if($protoDune==1){
    $APAToSouthCrWall   =     27.7; #center APA to cryostat
    $CPAToWestCrWall    =     45.5; #center CPA to beam window side
    $DetCentToFloor     =     376.0; # center CPA to cryostat floor, y axis
    $CPACentToWestWall  =     396.2; # center CPA to beam window side, z axis
    $CPACentToEastWall  =     458.6;
    $CPACentToSaleWall  =     391.4; # center CPA to Saleve, x axis
    $CPACentToJuraWall  =     463.4; # center CPA to Jura, x axis
    $APAGap_y           =     0;

#Beam direction
#Distance to upstream wall  3962 mm
#Distance to downstream wall 4586 mm
#Distance wall to wall  8548 mm

#Drift direction
#Distance to Jura wall 4634 mm
#Distance to Saleve wall 3914 mm
#Distance wall to wall  8548 mm

#Vertical
#Distance to floor  3760.2 mm
#Distance to ceiling  4139.8 mm
#Floor to ceiling  7900
}

 #InnerDrift is the max distance form the edge of the CPA to the edge of the first wire plane
#$InnerDrift              = 359.4;
$InnerDrift              = 357.135;
$APAFrame_x              = 3*$inch; # ~2in -- this does not include the wire spacing

$TPCWireThickness        = 0.015;
$TPCWirePlaneThickness   = $TPCWireThickness;
#$APAWirePlaneSpacing     = 0.4730488 + $TPCWirePlaneThickness; # center to center spacing between all of the wire planes (g, u, v, and x)

if($Pitch3mmVersion==1){
    $APAWirePlaneSpacing     = 0.3;
} else {
    $APAWirePlaneSpacing     = 3/16*$inch;
}

# At creation of the plane volumes, the y and z boundaries will be increased
# by this much at each of the 4 edges. this is so the corners of the wire 
# tubes don't extrude. For all other purposes, the plane dimensions stay as originally defined
$UVPlaneBoundNudge = $TPCWireThickness;


# APA z-dimensions:
# The following are all widths about the z center of an APA
$APAphys_z    = 231.66;
$Zactive_z    = ($nZChannelsPerAPA/2-1)*$XWirePitch;# + $TPCWireThickness;
$APAFrame_z   = $APAphys_z - 2*(2*$G10thickness+$WrapCover);
$Vactive_z    = $APAFrame_z;
$Uactive_z    = $APAFrame_z + 2*$G10thickness;
$APAGap_z     =    0.4;  #separation between APAs (cover to cover) along the vertical axis

print "Zactive_z=".$Zactive_z."\n";
print "APAFrame_z=".$APAFrame_z."\n";
#print "=".$."\n";
#print "=".$."\n";


# APA y-dimensions:
$ReadoutBoardOverlap = 7.61; #board overlaps wires, chop this off of their active height
$APAFrame_y = 606.0;
$Zactive_y    = $APAFrame_y + 0*$G10thickness - $ReadoutBoardOverlap;
$Vactive_y    = $APAFrame_y + 1*$G10thickness - $ReadoutBoardOverlap; 
$Uactive_y    = $APAFrame_y + 2*$G10thickness - $ReadoutBoardOverlap;
  # the last G10 board for the grid, then a cover. This is not "covered" by the board
$APAphys_y    = $APAFrame_y + 4*$G10thickness + $WrapCover; 
#$APAGap_y     =    0.4;  #separation between APAs (cover to cover) along the incident beam axis

  # include APA spacing in y and z so volTPCs touch in y and z directions with correct APA
  # spacing - this makes for smoother event generation. 

print "Zactive_y=".$Zactive_y."\n";
print "APAphys_y=".$APAphys_y."\n";
#print "=".$."\n";
#print "=".$."\n";

$APA_UtoU_x   = $APAFrame_x + 6*$APAWirePlaneSpacing + (6-1)*$TPCWirePlaneThickness; # thickness of APA between center U wire to center U wire
#$SpaceAPAToCryoWall      = $APAToSouthCrWall - ($APA_UtoU_x + $TPCWireThickness)/2;
#$OuterDrift              = $SpaceAPAToCryoWall - 20.; #Subtract some to avoid overlaping.
                                                                               # outer wire planes (center to center)
#$TPCInner_x   = $InnerDrift + $APAWirePlaneSpacing;# + $TPCWirePlaneThickness;
$TPCInner_x   = $InnerDrift + 4*$APAWirePlaneSpacing + 4*$TPCWirePlaneThickness;
#$TPCOuter_x   = $OuterDrift + $APAWirePlaneSpacing;# + $TPCWirePlaneThickness; Making it smaller than the distance to the wall.
$TPCOuter_x   = 4*$APAWirePlaneSpacing + 4*$TPCWirePlaneThickness + 8;

print "TPCInner_x=".$TPCInner_x."\n";
print "TPCOuter_x=".$TPCOuter_x."\n";
print "APA_UtoU_x=".$APA_UtoU_x."\n";

$TPC_z    =   $APAphys_z;# + $APAGap_z;
$TPC_y    =   $APAphys_y;# + $APAGap_y;

#print "TPC_x=".$TPC_x."\n";
print "TPC_y=".$TPC_y."\n";
print "TPC_z=".$TPC_z."\n";

#$CPATube_OD = 5.08;
#$CPATube_ID = 4.747;

#$Cathode_x                 =    0.016;
$Cathode_x                 =    0.13*$inch;
#$Cathode_y                 =    $APAphys_y - $CPATube_OD;
#$Cathode_z                 =    $APAphys_z - $CPATube_OD;   
#$Cathode_y                 =    610.4;
$Cathode_z                 =    715.5;

$APAToAPA =   $APAFrame_x 
            + 2*$TPCInner_x
            + $Cathode_x; # center to center

$CPAToAPA =   $APAFrame_x/2
            + $TPCInner_x #2*$APAWirePlaneSpacing
            + $Cathode_x/2; # center to center

$SpaceAPAToCryoWall      =   $CPACentToSaleWall
                           - $CPAToAPA
                           - ($APA_UtoU_x + $TPCWireThickness)/2;

print "SpaceAPAToCryoWall=".$SpaceAPAToCryoWall."\n";
print "APAToAPA=".$APAToAPA."\n";
print "CPAToAPA=".$CPAToAPA."\n";

$SteelThickness         =       0.12; #half inch
$HeightGaseousAr        =       80;


#$Argon_x                =       854.8; #inside dimension of the cryostat 
#$Argon_x               =       ($nAPAWide-1)*$APAToAPA
#                              + $APA_UtoU_x + $TPCWirePlaneThickness
#                              + 2*$SpaceAPAToCryoWall;

if($workspace==1){ # this is arbitrary for the workspace, but size down a little
$Argon_x = 2*$CPAToAPA + $Cathode_x + 2*$SpaceAPAToCryoWall;
}

#$Argon_y               =       790.0; #inside dimension of the cryostat
#$Argon_y               =       $nAPAHigh*$APAphys_y 
#                              + ($nAPAHigh-1)*$APAGap_y
#                              + $SpaceAPAToFloor + $SpaceAPAToTopLAr 
#                              + $HeightGaseousAr; 
                                    # both liquid AND gaseous argon
                                    #
#$Argon_z                =      854.8; #inside dimension of the cryostat 
#$Argon_z               =       $nAPALong*$APAphys_z 
#                              + ($nAPALong-1)*$APAGap_z
#                              + $UpstreamLArPadding + $DownstreamLArPadding;

print "Argon_x=".$Argon_x."\n";
print "Argon_y=".$Argon_y."\n";
print "Argon_z=".$Argon_z."\n";
print "APAToAPA =".$APAToAPA."\n";
print "APAToSouthCrWall =".$APAToSouthCrWall."\n";

$APAToNorthCrWall       =       $Argon_x 
                              - $APAToAPA - $APAToSouthCrWall;

#$CPAToEastCrWall       =       $Argon_z 
#                             - ($Cathode_z - $CPATube_OD) - $CPAToWestCrWall;

$SpaceAPAToFloor        =       $DetCentToFloor - $APAphys_y/2;   

$SpaceAPAToTopLAr       =       $Argon_y
                              - $APAphys_y
                              - $SpaceAPAToFloor
                              - $HeightGaseousAr;

$UpstreamLArPadding     =       $CPACentToWestWall +
                              - ($nAPALong*$APAphys_z + ($nAPALong-1)*$APAGap_z)/2;

$DownstreamLArPadding   =       $CPACentToEastWall +
                              - ($nAPALong*$APAphys_z + ($nAPALong-1)*$APAGap_z)/2;

$Cryostat_x             =       $Argon_x + 2*$SteelThickness;
$Cryostat_y             =       $Argon_y + 2*$SteelThickness;
$Cryostat_z             =       $Argon_z + 2*$SteelThickness;

print "Cryo_x=".$Cryostat_x."\n";
print "Cryo_y=".$Cryostat_y."\n";
print "Cryo_z=".$Cryostat_z."\n";
##################################################################
############## DetEnc and World relevant parameters  #############


$SteelSupport_x        =        62.8;
$SteelSupport_y        =        62.8;
$SteelSupport_z        =        62.8;
$SteelPlate            =         1.0;
$FoamPadding           =        80. - $SteelThickness;
$FracVolOfSteel        =        0.08; #The steel support is not a solid block, but a mixture of air and steel
$FracMassOfSteel       =        $FracVolOfSteel/($FracVolOfSteel+1.205/7930*(1-$FracVolOfSteel)); #The steel support is not a solid block, but a mixture of air and steel

$SpaceSteelSupportToWall    = 900;
$SpaceSteelSupportToCeiling = 900;

$DetEncWidth            =    ($Cryostat_x
                              + 2*($SteelSupport_x + $FoamPadding) + 2*$SpaceSteelSupportToCeiling);
$DetEncHeight           =    ($Cryostat_y
                              + 2*($SteelSupport_y + $FoamPadding) + 2*$SpaceSteelSupportToWall);
$DetEncLength           =    ($Cryostat_z
                              + 2*($SteelSupport_z + $FoamPadding) + 2*$SpaceSteelSupportToWall);

$posCryoInDetEnc_x = 0;
#$posCryoInDetEnc_y = - $DetEncHeight/2 + $SteelSupport_y + $FoamPadding + $Cryostat_y/2;
$posCryoInDetEnc_y = 0;

$RockThickness          =       400;

  # We want the world origin to be at the very front of the fiducial volume.
  # move it to the front of the enclosure, then back it up through the concrete/foam, 
  # then through the Cryostat shell, then through the upstream dead LAr (including the
  # dead LAr on the edge of the TPC, but this is covered in $UpstreamLArPadding).
  # This is to be added to the z position of every volume in volWorld

$OriginZSet             =      $APAphys_z
                             + $APAGap_z
                             + $Uactive_z/2
                             + ($CPACentToEastWall - $CPACentToWestWall)/2; 
$ZAssym=($CPACentToEastWall - $CPACentToWestWall)/2;
print "Z assym ".$ZAssym."\n";
#$OriginZSet             =       $DetEncLength/2 
#                              - $SpaceSteelSupportToWall
#                              - $SteelSupport_z
#                              - $FoamPadding
#                              - $SteelThickness
#                              - $UpstreamLArPadding
#                              - ($APAphys_z - $Uactive_z)/2;

  # We want the world origin to be vertically centered between the stacked APAs.
  # This is to be added to the y position of every volume in volWorld

$OriginYSet             =       $APAphys_y/2
                              + $Argon_y/2
                              - $DetCentToFloor
#                              + $SpaceSteelSupportToWall;
                              - $posCryoInDetEnc_y;

print "DetCentToFloor = ".$DetCentToFloor."\n";                                          
#$OriginYSet             =       $DetEncHeight/2
#                              - $SteelSupport_y
#                              - $FoamPadding
#                              - $SteelThickness
#                              - $SpaceAPAToFloor
#                              - $APAphys_y
#                              - $APAGap_y/2;

#if($protoDune==1){ # bring the origin to the bottom of the APAs for protoDUNE
#    $OriginYSet = $OriginYSet + $APAphys_y + $APAGap_y/2;
#}

#$OriginXSet             =        ($APAToNorthCrWall - $APAToSouthCrWall)/2;
$OriginXSet             =        ($CPACentToJuraWall - $CPACentToSaleWall)/2;

print "OriginXSet =".$OriginXSet.", OriginYSet =".$OriginYSet.", OriginZSet =".$OriginZSet."\n";

#$OriginXSet             =     0; # centered for now

#$Cathode_z    =  $Cathode_z/3 - $CPATube_OD;

###################################################################
######################## Beam Window 2 Parameters ###################

$thetaYZ = 11.342;
$theta2XZ = 8.189;

$BeamTheta2 = atan (sqrt(tan(deg2rad($theta2XZ))**2 +tan(deg2rad($thetaYZ))**2));
$BeamPhi2   = atan (tan(deg2rad($thetaYZ))/tan(deg2rad($theta2XZ)));

$thetaYZprime = rad2deg(atan(sin($BeamTheta2)*sin($BeamPhi2+deg2rad(180))/sqrt(cos($BeamTheta2)**2 + sin($BeamTheta2)**2*cos($BeamPhi2)**2)));

print "thetaYZprime =".$thetaYZprime."\n";

$DeltaXZ2      = tan($BeamTheta2)*cos($BeamPhi2);
$DeltaYZ2      = tan($BeamTheta2)*sin($BeamPhi2);

print "DeltaXZ2 = ".$DeltaXZ2."\n";
print "DeltaYZ2 = ".$DeltaYZ2."\n";

$BeamTheta2Deg  = rad2deg($BeamTheta2);
$BeamPhi2Deg    = rad2deg($BeamPhi2);   

######################### Beam Window 3 Parameters ###################

$thetaYZ  = 11.342;
$theta3XZ = 11.844;

$BeamTheta3 = atan (sqrt(tan(deg2rad($theta3XZ))**2 +tan(deg2rad($thetaYZ))**2));
$BeamPhi3   = atan (tan(deg2rad($thetaYZ))/tan(deg2rad($theta3XZ)));

print "BeamTheta3 = ".$BeamTheta3."\n";
print "BeamPhi3 = ".$BeamPhi3."\n";


$thetaYZ3prime = rad2deg(atan(sin($BeamTheta3)*sin($BeamPhi3+deg2rad(180))/sqrt(cos($BeamTheta3)**2 + sin($BeamTheta3)**2*cos($BeamPhi3)**2)));

print "thetaYZ3prime =".$thetaYZ3prime."\n";

$DeltaXZ3      = tan($BeamTheta3)*cos($BeamPhi3);
$DeltaYZ3      = tan($BeamTheta3)*sin($BeamPhi3);

$BeamPipeRad      =   12.5;
$BeamVaPipeRad    =   $BeamPipeRad - 0.2;
$BeamPipeLe       =   900.0;
$BeamVaPipeLe     =   $BeamPipeLe;
$BeamWFoLe        =   52.0;
$BeamWGlLe        =   10.0;
$BeamWStPlateFF_x = - 41.0;
$BeamWStPlateFF_y =  125.7;
$BeamWStPlateFF_z =  - ($Cryostat_z/2 + $FoamPadding + $SteelPlate);
print "BeamWStPlateFF_z = ".$BeamWStPlateFF_z."\n";
#$BeamWStPlateFF_z = -508.4; 

$BeamWStPlateLe   = $SteelPlate/cos($BeamTheta3)+0.001;
$BeamWStPlate_x   = $BeamWStPlateFF_x - ($SteelPlate/2)*$DeltaXZ3;
$BeamWStPlate_y   = $BeamWStPlateFF_y - ($SteelPlate/2)*$DeltaYZ3;
$BeamWStPlate_z   = $BeamWStPlateFF_z + $SteelPlate/2;

$BeamWFoRemLe     = $FoamPadding/cos($BeamTheta3)+0.001;
$BeamWFoRemPosDZ  = $SteelPlate + $FoamPadding/2;
$BeamWFoRem_x     = $BeamWStPlateFF_x - $BeamWFoRemPosDZ*$DeltaXZ3;
$BeamWFoRem_y     = $BeamWStPlateFF_y - $BeamWFoRemPosDZ*$DeltaYZ3;
$BeamWFoRem_z     = $BeamWStPlateFF_z + $BeamWFoRemPosDZ;

$BeamWStSuLe      = ($SteelSupport_z - $SteelPlate)/cos($BeamTheta3)+0.001;
$BeamWStSuPosDZ   = - ($SteelSupport_z - $SteelPlate)/2;   # going upstream from the steel plate 
$BeamWStSu_x      = $BeamWStPlateFF_x - $BeamWStSuPosDZ*$DeltaXZ3;
$BeamWStSu_y      = $BeamWStPlateFF_y - $BeamWStSuPosDZ*$DeltaYZ3;
$BeamWStSu_z      = $BeamWStPlateFF_z + $BeamWStSuPosDZ;

$BeamWFoPosDZ     = $SteelPlate + $FoamPadding - $BeamWFoLe*cos($BeamTheta3)/2;
$BeamWFo_x        = $BeamWStPlateFF_x - $BeamWFoPosDZ*$DeltaXZ3;
$BeamWFo_y        = $BeamWStPlateFF_y - $BeamWFoPosDZ*$DeltaYZ3 + $posCryoInDetEnc_y;
$BeamWFo_z        = $BeamWStPlateFF_z + $BeamWFoPosDZ;

$BeamWGlPosDZ     = $SteelPlate + $FoamPadding - ($BeamWFoLe + $BeamWGlLe/2)*cos($BeamTheta3);
$BeamWGl_x        = $BeamWStPlateFF_x - $BeamWGlPosDZ*$DeltaXZ3;
$BeamWGl_y        = $BeamWStPlateFF_y - $BeamWGlPosDZ*$DeltaYZ3 + $posCryoInDetEnc_y;
$BeamWGl_z        = $BeamWStPlateFF_z + $BeamWGlPosDZ;

$BeamWVaPosDZ     = $SteelPlate + $FoamPadding - ($BeamWFoLe + $BeamWGlLe + $BeamPipeLe/2)*cos($BeamTheta3);
$BeamWVa_x        = $BeamWStPlateFF_x - $BeamWVaPosDZ*$DeltaXZ3; 
$BeamWVa_y        = $BeamWStPlateFF_y - $BeamWVaPosDZ*$DeltaYZ3 + $posCryoInDetEnc_y;
$BeamWVa_z        = $BeamWStPlateFF_z + $BeamWVaPosDZ;

$BeamPlugRad      = 10.48;
$BeamPlugNiRad    = 9.72;
$BeamPlugUSAr     = 1/cos($BeamTheta3); # 1 cm US LAr layer between beam plug and primary membrane
$BeamPlugLe       = (44.6)/cos($BeamTheta3) - $BeamPlugUSAr; #with current geometry and 49.22 Dz the flange's front face just touches the active volume.
$BeamPlugNiLe     = $BeamPlugLe-0.59/cos($BeamTheta3);

$BeamPlugPosDZ    = $SteelPlate + $FoamPadding + $SteelThickness + $BeamPlugUSAr + $BeamPlugLe*cos($BeamTheta3)/2;
$BeamPlug_x       = $BeamWStPlateFF_x - $BeamPlugPosDZ*$DeltaXZ3;
$BeamPlug_y       = $BeamWStPlateFF_y - $BeamPlugPosDZ*$DeltaYZ3;
$BeamPlug_z       = $BeamWStPlateFF_z + $BeamPlugPosDZ;

$BePlFlangePosDZ    = $SteelPlate + $FoamPadding + $SteelThickness + $BeamPlugUSAr + $BeamPlugLe*cos($BeamTheta3); #This is Dz to the end of the beam plug pipe needed for x,y position.
$BePlFlange_x       = $BeamWStPlateFF_x - $BePlFlangePosDZ*$DeltaXZ3;
$BePlFlange_y       = $BeamWStPlateFF_y - $BePlFlangePosDZ*$DeltaYZ3;
$BePlFlange_z       = $BeamWStPlateFF_z + $BePlFlangePosDZ + 1.8;  # Adding the half the thickness of the flange.

#$BeamPlugERingPosZ = -$BeamPlugLe/2 + 6.83;
$BeamPlugERingPosZ = -$BeamPlugLe/2 + 5.5;
$BeamPlugNiPos_z   = 0.59/2/cos($BeamTheta3);

print "BeamPlugLe = $BeamPlugLe"."\n";
print "BeamTheta3 = $BeamTheta3"."\n";
print "BeamWStPlate x=$BeamWStPlate_x".", y=$BeamWStPlate_y".", z=$BeamWStPlate_z"."\n";
print "BeamWStSu x=$BeamWStSu_x".", y=$BeamWStSu_y".", z=$BeamWStSu_z"."\n";
print "BeamWFoRem x=$BeamWFoRem_x".", y=$BeamWFoRem_y".", z=$BeamWFoRem_z"."\n";
print "BeamWFo x=$BeamWFo_x".", y=$BeamWFo_y".", z=$BeamWFo_z"."\n";
print "BeamWGl x=$BeamWGl_x".", y=$BeamWGl_y".", z=$BeamWGl_z"."\n";
print "BeamWVa x=$BeamWVa_x".", y=$BeamWVa_y".", z=$BeamWVa_z"."\n";
print "CPAToWestCrWall=$CPAToWestCrWall"."\n";
print "BeamPlug x=$BeamPlug_x".", y=$BeamPlug_y".", z=$BeamPlug_z"."\n";
print "BeamPlugFlange x=$BePlFlange_x".", y=$BePlFlange_y".", z=$BePlFlange_z"."\n";

$BeamTheta3Deg  = rad2deg($BeamTheta3);
$BeamPhi3Deg    = rad2deg($BeamPhi3);


$BWFFCoord3X   =       $BeamWStPlateFF_x - $BeamWStSuPosDZ*$DeltaXZ3*2 
                    +  $OriginXSet;
$BWFFCoord3Y   =       $BeamWStPlateFF_y - $BeamWStSuPosDZ*$DeltaYZ3*2
                    +  $OriginYSet + $posCryoInDetEnc_y;   
$BWFFCoord3Z   =    -  ($Cryostat_z/2 + $SteelSupport_z + $FoamPadding)
                    +  $OriginZSet;
print "BeamWStPlateFF_x=".$BeamWStPlateFF_x.", BeamWStPlateFF_y=".$BeamWStPlateFF_y."\n";
print "BeamWStSuPosDZ= $BeamWStSuPosDZ; "." DeltaYZ3=".$DeltaYZ3."\n";
print "BWFFCoord3X =".$BWFFCoord3X."\n";
print "BWFFCoord3Y =".$BWFFCoord3Y."\n";
print "BWFFCoord3Z =".$BWFFCoord3Z."\n";

$BW3StPlCoordX =  $BeamWStPlateFF_x + $OriginXSet;
$BW3StPlCoordY =  $BeamWStPlateFF_y + $OriginYSet + $posCryoInDetEnc_y;
$BW3StPlCoordZ =  $BeamWStPlateFF_z + $OriginZSet;

print "BW3StPlCoordX =".$BW3StPlCoordX."\n";
print "BW3StPlCoordY =".$BW3StPlCoordY."\n";
print "BW3StPlCoordZ =".$BW3StPlCoordZ."\n";


###################################################################
########################  CRT Dimensions ##########################

$CRTPaddleWidth      =  5.0;    
$CRTPaddleHeight     =  1.0;
$CRTPaddleLength     =  322.5;

$CRTModWidth         =  162.5;
$CRTModHeight        =  2.0;
$CRTModLength        =  322.5;

# SuperModule Centers as per Survey Document
$CRT_DSTopLeft_x     =  171.2;
$CRT_DSTopLeft_y     = -473.88;
$CRT_DSTopLeftFr_z   =  1042.13;
$CRT_DSTopLeftBa_z   =  1050.13;

$CRT_DSBotLeft_x     =  176.51;
$CRT_DSBotLeft_y     = -840.6;
$CRT_DSBotLeftFr_z   =  1041.74;
$CRT_DSBotLeftBa_z   =  1050.13;

$CRT_DSTopRight_x    = -176.23;
$CRT_DSTopRight_y    = -474.85;
$CRT_DSTopRightFr_z  =  1042.64;
$CRT_DSTopRightBa_z  =  1050.85;

$CRT_DSBotRight_x    = -169.6;
$CRT_DSBotRight_y    = -840.55;
$CRT_DSBotRightFr_z  =  1042.88;
$CRT_DSBotRightBa_z  =  1051.93;

$CRT_USTopLeft_x     =  393.6;
$CRT_USTopLeft_y     = -401.33;
$CRT_USTopLeftFr_z   = -295.05;
$CRT_USTopLeftBa_z   = -286.85;

$CRT_USBotLeft_x     =  394.14;
$CRT_USBotLeft_y     = -734.48;
$CRT_USBotLeftFr_z   = -320.24;
$CRT_USBotLeftBa_z   = -310.88;

$CRT_USTopRight_x    = -38.85;
$CRT_USTopRight_y    = -400.85;
$CRT_USTopRightFr_z  = -998.95;
$CRT_USTopRightBa_z  = -990.97;

$CRT_USBotRight_x    = -31.47;
$CRT_USBotRight_y    = -735.13;
$CRT_USBotRightFr_z  = -1022.25;
$CRT_USBotRightBa_z  = -1015.01;

# Coordinates of Survey origin w.r.t. Steel Plate box center as per penetrations steel plate drawing (flange 1.2)
$CRTSurveyOrigin_x   = -36.0;
$CRTSurveyOrigin_y   =  534.43;
$CRTSurveyOrigin_z   = -344.1;

# Distance between CRT module centers and CRT SuperModule Centers
$ModuleSMDist        = 85.6;
$ModuleOff_z         = 1;    # approx. correction for the center of a Module. Survey measures Z to the outside surface. Negative for the most US CRTs (surveyed from behind).
$ModuleLongCorr      = 5.6;  # allign the the modules at the frame's edge

# Beam Spot on the inside of the cryostat

$BeamSpotDSS_x          = -20.58;
$BeamSpotDSS_y          = -425.41;
$BeamSpotDSS_z          = -82.96;

$BeamSpot_x = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $BeamSpotDSS_x +  $OriginXSet;;
$BeamSpot_y = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $BeamSpotDSS_y +  $OriginYSet;;
$BeamSpot_z = $posCryoInDetEnc_z + $CRTSurveyOrigin_z + $BeamSpotDSS_z +  $OriginZSet;;

print "BeamSpot_x =".$BeamSpot_x.", BeamSpot_y =".$BeamSpot_y.", BeamSpot_z =".$BeamSpot_z."\n";

####################### End of Survey data ##########


my @posCRTDS_x = ();
my @posCRTDS_y = ();
my @posCRTDS_z = ();
my @posCRTDS_rot = ();

$posCRTDS_x[0] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSTopLeft_x - $ModuleSMDist;
$posCRTDS_y[0] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSTopLeft_y - $ModuleLongCorr;
$posCRTDS_z[0] = $CRTSurveyOrigin_z + $CRT_DSTopLeftBa_z + $ModuleOff_z;
$posCRTDS_rot[0] = "rPlus90AboutX";

$posCRTDS_x[1] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSTopLeft_x + $ModuleSMDist;
$posCRTDS_y[1] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSTopLeft_y - $ModuleLongCorr;
$posCRTDS_z[1] = $CRTSurveyOrigin_z + $CRT_DSTopLeftBa_z + $ModuleOff_z;
$posCRTDS_rot[1] = "rPlus90AboutX";

$posCRTDS_x[2] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSTopLeft_x - $ModuleLongCorr;
$posCRTDS_y[2] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSTopLeft_y - $ModuleSMDist;
$posCRTDS_z[2] = $CRTSurveyOrigin_z + $CRT_DSTopLeftFr_z + $ModuleOff_z;
$posCRTDS_rot[2] = "rMinus90AboutYMinus90AboutX";

$posCRTDS_x[3] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSTopLeft_x - $ModuleLongCorr;
$posCRTDS_y[3] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSTopLeft_y + $ModuleSMDist;
$posCRTDS_z[3] = $CRTSurveyOrigin_z + $CRT_DSTopLeftFr_z + $ModuleOff_z;
$posCRTDS_rot[3] = "rMinus90AboutYMinus90AboutX";

$posCRTDS_x[4] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSBotLeft_x - $ModuleSMDist;
$posCRTDS_y[4] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSBotLeft_y + $ModuleLongCorr;
$posCRTDS_z[4] = $CRTSurveyOrigin_z + $CRT_DSBotLeftFr_z + $ModuleOff_z;
$posCRTDS_rot[4] = "rPlus90AboutX";

$posCRTDS_x[5] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSBotLeft_x + $ModuleSMDist;
$posCRTDS_y[5] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSBotLeft_y + $ModuleLongCorr;
$posCRTDS_z[5] = $CRTSurveyOrigin_z + $CRT_DSBotLeftFr_z + $ModuleOff_z;
$posCRTDS_rot[5] = "rPlus90AboutX";

$posCRTDS_x[6] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSBotLeft_x - $ModuleLongCorr;
$posCRTDS_y[6] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSBotLeft_y - $ModuleSMDist;
$posCRTDS_z[6] = $CRTSurveyOrigin_z + $CRT_DSBotLeftBa_z + $ModuleOff_z;
$posCRTDS_rot[6] = "rMinus90AboutYMinus90AboutX";

$posCRTDS_x[7] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSBotLeft_x - $ModuleLongCorr;
$posCRTDS_y[7] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSBotLeft_y + $ModuleSMDist;
$posCRTDS_z[7] = $CRTSurveyOrigin_z + $CRT_DSBotLeftBa_z + $ModuleOff_z;
$posCRTDS_rot[7] = "rMinus90AboutYMinus90AboutX";

$posCRTDS_x[8] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSTopRight_x - $ModuleSMDist;
$posCRTDS_y[8] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSTopRight_y - $ModuleLongCorr;
$posCRTDS_z[8] = $CRTSurveyOrigin_z + $CRT_DSTopRightFr_z + $ModuleOff_z;
$posCRTDS_rot[8] = "rPlus90AboutX";

$posCRTDS_x[9] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSTopRight_x + $ModuleSMDist;
$posCRTDS_y[9] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSTopRight_y - $ModuleLongCorr;
$posCRTDS_z[9] = $CRTSurveyOrigin_z + $CRT_DSTopRightFr_z + $ModuleOff_z;
$posCRTDS_rot[9] = "rPlus90AboutX";

$posCRTDS_x[10] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSTopRight_x + $ModuleLongCorr;
$posCRTDS_y[10] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSTopRight_y - $ModuleSMDist;
$posCRTDS_z[10] = $CRTSurveyOrigin_z + $CRT_DSTopRightBa_z + $ModuleOff_z;
$posCRTDS_rot[10] = "rMinus90AboutYMinus90AboutX";

$posCRTDS_x[11] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSTopRight_x + $ModuleLongCorr;
$posCRTDS_y[11] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSTopRight_y + $ModuleSMDist;
$posCRTDS_z[11] = $CRTSurveyOrigin_z + $CRT_DSTopRightBa_z + $ModuleOff_z;
$posCRTDS_rot[11] = "rMinus90AboutYMinus90AboutX";

$posCRTDS_x[12] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSBotRight_x - $ModuleSMDist;
$posCRTDS_y[12] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSBotRight_y + $ModuleLongCorr;
$posCRTDS_z[12] = $CRTSurveyOrigin_z + $CRT_DSBotRightBa_z + $ModuleOff_z;
$posCRTDS_rot[12] = "rPlus90AboutX";

$posCRTDS_x[13] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSBotRight_x + $ModuleSMDist;
$posCRTDS_y[13] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSBotRight_y + $ModuleLongCorr;
$posCRTDS_z[13] = $CRTSurveyOrigin_z + $CRT_DSBotRightBa_z + $ModuleOff_z;
$posCRTDS_rot[13] = "rPlus90AboutX";

$posCRTDS_x[14] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSBotRight_x + $ModuleLongCorr;
$posCRTDS_y[14] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSBotRight_y - $ModuleSMDist;
$posCRTDS_z[14] = $CRTSurveyOrigin_z + $CRT_DSBotRightFr_z + $ModuleOff_z;
$posCRTDS_rot[14] = "rMinus90AboutYMinus90AboutX";

$posCRTDS_x[15] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_DSBotRight_x + $ModuleLongCorr;
$posCRTDS_y[15] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_DSBotRight_y + $ModuleSMDist;
$posCRTDS_z[15] = $CRTSurveyOrigin_z + $CRT_DSBotRightFr_z + $ModuleOff_z;
$posCRTDS_rot[15] = "rMinus90AboutYMinus90AboutX";


my @posCRTUS_x = ();
my @posCRTUS_y = ();
my @posCRTUS_z = ();
my @posCRTUS_rot = ();


$posCRTUS_x[0] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USTopLeft_x - $ModuleSMDist;
$posCRTUS_y[0] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USTopLeft_y - $ModuleLongCorr;
$posCRTUS_z[0] = $CRTSurveyOrigin_z + $CRT_USTopLeftBa_z + $ModuleOff_z;
$posCRTUS_rot[0] = "rPlus90AboutX";

$posCRTUS_x[1] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USTopLeft_x + $ModuleSMDist;
$posCRTUS_y[1] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USTopLeft_y - $ModuleLongCorr;
$posCRTUS_z[1] = $CRTSurveyOrigin_z + $CRT_USTopLeftBa_z + $ModuleOff_z;
$posCRTUS_rot[1] = "rPlus90AboutX";

$posCRTUS_x[2] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USTopLeft_x - $ModuleLongCorr;
$posCRTUS_y[2] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USTopLeft_y - $ModuleSMDist;
$posCRTUS_z[2] = $CRTSurveyOrigin_z + $CRT_USTopLeftFr_z + $ModuleOff_z;
$posCRTUS_rot[2] = "rMinus90AboutYMinus90AboutX";

$posCRTUS_x[3] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USTopLeft_x - $ModuleLongCorr;
$posCRTUS_y[3] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USTopLeft_y + $ModuleSMDist;
$posCRTUS_z[3] = $CRTSurveyOrigin_z + $CRT_USTopLeftFr_z + $ModuleOff_z;
$posCRTUS_rot[3] = "rMinus90AboutYMinus90AboutX";

$posCRTUS_x[4] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USBotLeft_x - $ModuleSMDist;
$posCRTUS_y[4] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USBotLeft_y + $ModuleLongCorr;
$posCRTUS_z[4] = $CRTSurveyOrigin_z + $CRT_USBotLeftFr_z + $ModuleOff_z;
$posCRTUS_rot[4] = "rPlus90AboutX";

$posCRTUS_x[5] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USBotLeft_x + $ModuleSMDist;
$posCRTUS_y[5] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USBotLeft_y + $ModuleLongCorr;
$posCRTUS_z[5] = $CRTSurveyOrigin_z + $CRT_USBotLeftFr_z + $ModuleOff_z;
$posCRTUS_rot[5] = "rPlus90AboutX";

$posCRTUS_x[6] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USBotLeft_x - $ModuleLongCorr;
$posCRTUS_y[6] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USBotLeft_y - $ModuleSMDist;
$posCRTUS_z[6] = $CRTSurveyOrigin_z + $CRT_USBotLeftBa_z + $ModuleOff_z;
$posCRTUS_rot[6] = "rMinus90AboutYMinus90AboutX";

$posCRTUS_x[7] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USBotLeft_x - $ModuleLongCorr;
$posCRTUS_y[7] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USBotLeft_y + $ModuleSMDist;
$posCRTUS_z[7] = $CRTSurveyOrigin_z + $CRT_USBotLeftBa_z + $ModuleOff_z;
$posCRTUS_rot[7] = "rMinus90AboutYMinus90AboutX";

$posCRTUS_x[8] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USTopRight_x - $ModuleSMDist;
$posCRTUS_y[8] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USTopRight_y - $ModuleLongCorr;
$posCRTUS_z[8] = $CRTSurveyOrigin_z + $CRT_USTopRightFr_z - $ModuleOff_z;
$posCRTUS_rot[8] = "rPlus90AboutX";

$posCRTUS_x[9] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USTopRight_x + $ModuleSMDist;
$posCRTUS_y[9] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USTopRight_y - $ModuleLongCorr;
$posCRTUS_z[9] = $CRTSurveyOrigin_z + $CRT_USTopRightFr_z - $ModuleOff_z;
$posCRTUS_rot[9] = "rPlus90AboutX";

$posCRTUS_x[10] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USTopRight_x + $ModuleLongCorr;
$posCRTUS_y[10] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USTopRight_y - $ModuleSMDist;
$posCRTUS_z[10] = $CRTSurveyOrigin_z + $CRT_USTopRightBa_z - $ModuleOff_z;
$posCRTUS_rot[10] = "rMinus90AboutYMinus90AboutX";

$posCRTUS_x[11] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USTopRight_x + $ModuleLongCorr;
$posCRTUS_y[11] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USTopRight_y + $ModuleSMDist;
$posCRTUS_z[11] = $CRTSurveyOrigin_z + $CRT_USTopRightBa_z - $ModuleOff_z;
$posCRTUS_rot[11] = "rMinus90AboutYMinus90AboutX";

$posCRTUS_x[12] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USBotRight_x - $ModuleSMDist;
$posCRTUS_y[12] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USBotRight_y + $ModuleLongCorr;
$posCRTUS_z[12] = $CRTSurveyOrigin_z + $CRT_USBotRightBa_z - $ModuleOff_z;
$posCRTUS_rot[12] = "rPlus90AboutX";

$posCRTUS_x[13] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USBotRight_x + $ModuleSMDist;
$posCRTUS_y[13] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USBotRight_y + $ModuleLongCorr;
$posCRTUS_z[13] = $CRTSurveyOrigin_z + $CRT_USBotRightBa_z - $ModuleOff_z;
$posCRTUS_rot[13] = "rPlus90AboutX";

$posCRTUS_x[14] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USBotRight_x + $ModuleLongCorr;
$posCRTUS_y[14] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USBotRight_y - $ModuleSMDist;
$posCRTUS_z[14] = $CRTSurveyOrigin_z + $CRT_USBotRightFr_z - $ModuleOff_z;
$posCRTUS_rot[14] = "rMinus90AboutYMinus90AboutX";

$posCRTUS_x[15] = $posCryoInDetEnc_x + $CRTSurveyOrigin_x + $CRT_USBotRight_x + $ModuleLongCorr;
$posCRTUS_y[15] = $posCryoInDetEnc_y + $CRTSurveyOrigin_y + $CRT_USBotRight_y + $ModuleSMDist;
$posCRTUS_z[15] = $CRTSurveyOrigin_z + $CRT_USBotRightFr_z - $ModuleOff_z;
$posCRTUS_rot[15] = "rMinus90AboutYMinus90AboutX";

###################################################################
########################  Paddle Dimensions ########################

$APAFrameZSide_y = 4*$inch;
$APAFrameYSide_z = 4*$inch;
$LightPaddleIU_x          =    0.863;
$LightPaddleIU_y          =    8.23;
$LightPaddleIU_z          =    200.074;
$LightPaddleMIT_x         =    0.551;
$LightPaddleMIT_y         =    8.37;
$LightPaddleMIT_z         =    203.925;
$LightPaddle_z            =    $APAFrame_z - 2*$APAFrameYSide_z; #Length of the APA Frame --> used for rails
$nLightPaddlesPerAPA      =    10;   # 10, or 20 for double coverage (for now)
$PaddleYInterval          =    59.2;
$FrameToPaddleSpace       =    30.412;
$SiPM_z                   =    0;

$LightPaddleIUOut_x       =    0.86;   #IU Rail
$LightPaddleMITOut_x      =    0.55;   #MIT Rail
$LightPaddleIUOut_y       =    8.47;   #IU Rail
$LightPaddleMITOut_y      =    8.61;   #MIT Rail
$LightPaddleOut_z         =    $LightPaddle_z;
$list_sub_paddle[0]=$APACenter_y-$APAphys_y/2+$FrameToPaddleSpace+$LightPaddle_y/2;
for($list=1; $list<10; $list++){
$list_sub_paddle[$list]=$list_pos[$i-1]+$LightPaddle_y/2+$PaddleYInterval;
}

# $PaddleYInterval is defined so that the center-to-center distance in the 
# y direction between paddles is uniform between vertically stacked APAs.
# $FrameToPaddleSpace is from the BOTTOM of the APA frame (4" in y direction)
# to the CENTER of a paddle, including the 4" part of the frame. This variable's
# primary purpose is to position the lowest paddle in each APA.


####################################################################
######################## ARAPUCA Dimensions ########################

$ArapucaOut_x = 16.0;
$ArapucaOut_y = 96.0;
$ArapucaOut_z = 10.0*$LightPaddle_z; 
$ArapucaIn_x = 15.0;
$ArapucaIn_y = 80.0;
$ArapucaIn_z = 100.0;
$ArapucaAcceptanceWindow_x =1.0;
$ArapucaAcceptanceWindow_y =80.0;
$ArapucaAcceptanceWindow_z =100.0;
$pos_subtraction_arapuca_x = 2.0;
$gapCenterLeft_arapuca_z = 50.0;#if not simmetrical, positioning of windows for odd APAs needs to change
$gapCenterRight_arapuca_z = 100.0-$gapCenterLeft_arapuca_z;
$nAraAPA1 = 0;
$nAraAPA2 = 3;
$nSlotAra1 = 6;
$nSlotAra2 = 4;

$list_pos[0]=-$gapCenterLeft_arapuca_z-0.5*$ArapucaIn_z;
for($list=1; $list<8; $list++){
$list_pos[$list]=$list_pos[$i-1]-1.0*$ArapucaIn_z-2.0;
}
$list_pos[8]=$gapCenterRight_arapuca_z+0.5*$ArapucaIn_z;
for($list=9; $list<16; $list++){
$list_pos[$list]=$list_pos[$i-1]+$ArapucaIn_z+2.0;
}

#nAraAPAi is the number of the APA that the ith arapuca bar goes into and nSlotArai, its slot in that APA. These numbers are valid for the protodune-sp geometry.




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
#++++++++++++++++++++++++++++++++++++++ gen_Extend +++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Extend()
{

# Create the <define> fragment file name, 
# add file to list of fragments,
# and open it
    $DEF = "subfile_protodune_v7_Ext" . $suffix . ".gdml";
    push (@gdmlFiles, $DEF);
    $DEF = ">" . $DEF;
    open(DEF) or die("Could not open file $DEF for writing");

print DEF <<EOF;
<?xml version='1.0'?>
<gdml>
<extension>
   <color name="magenta"     R="0.0"  G="1.0"  B="0.0"  A="1.0" />
   <color name="green"       R="0.0"  G="1.0"  B="0.0"  A="1.0" />
   <color name="red"         R="1.0"  G="0.0"  B="0.0"  A="1.0" />
   <color name="blue"        R="0.0"  G="0.0"  B="1.0"  A="1.0" />
   <color name="yellow"      R="1.0"  G="1.0"  B="0.0"  A="1.0" />    
</extension>
</gdml>
EOF
    close (DEF);
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_Define +++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_Define()
{

# Create the <define> fragment file name, 
# add file to list of fragments,
# and open it
    $DEF = "subfile_protodune_v7_Def" . $suffix . ".gdml";
    push (@gdmlFiles, $DEF);
    $DEF = ">" . $DEF;
    open(DEF) or die("Could not open file $DEF for writing");

$UWireAngle = 90-$UAngle;
$VWireAngle = 90+$VAngle;
print DEF <<EOF;
<?xml version='1.0'?>
<gdml>
<define>

<!--



-->

   <position name="posCryoInDetEnc"     unit="cm" x="0" y="$posCryoInDetEnc_y" z="0"/>
   <position name="posCenter"           unit="cm" x="0" y="0" z="0"/>
   <position name="posBeamWAr3"         unit="cm" x="$BeamPlug_x" y="$BeamPlug_y" z="$BeamPlug_z"/>
   <rotation name="rPlus90AboutX"       unit="deg" x="90" y="0" z="0"/>
   <rotation name="rMinus90AboutY"      unit="deg" x="0" y="270" z="0"/>
   <rotation name="rMinus90AboutYMinus90AboutX"       unit="deg" x="270" y="270" z="0"/>
   <rotation name="rPlusUAngleAboutX"   unit="deg" x="$UWireAngle" y="0"   z="0"/>
   <rotation name="rPlusVAngleAboutX"   unit="deg" x="$VWireAngle" y="0"   z="0"/>
   <rotation name="rPlus180AboutX"      unit="deg" x="180" y="0" z="0"/>
   <rotation name="rPlus180AboutY"      unit="deg" x="0" y="180" z="0"/>
   <rotation name="rPlus180AboutZ"      unit="deg" x="0" y="0" z="180"/>
   <rotation name="rPlus180AboutXPlus180AboutY" unit="deg" x="180" y="180"   z="0"/>
   <rotation name="rPlus90AboutXPlus90AboutZ"       unit="deg" x="90" y="0" z="90"/>
   <rotation name="rPlus180AboutXPlus180AboutZ"   unit="deg" x="180" y="0" z="180"/>
   <rotation name="rIdentity"           unit="deg" x="0" y="0"   z="0"/>
   <rotation name="rBeamW2"             unit="deg" x="0" y="-$BeamTheta2Deg" z="$BeamPhi2Deg"/>
   <rotation name="rBeamWRev2"          unit="deg" x="-11.342" y="8.03118195044" z="-55.1415281209"/>
   <rotation name="rBeamW3"             unit="deg" x="0" y="-$BeamTheta3Deg" z="$BeamPhi3Deg"/>
   <rotation name="rBeamWRev3"          unit="deg" x="-11.342" y="11.6190450403" z="-44.8829268772"/>
   <position name="posArapucaSub0"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[0]]}"/>
   <position name="posArapucaSub1"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[1]]}"/>
   <position name="posArapucaSub2"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[2]]}"/>
   <position name="posArapucaSub3"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[3]]}"/>
   <position name="posArapucaSub4"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[4]]}"/>
   <position name="posArapucaSub5"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[5]]}"/>
   <position name="posArapucaSub6"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[6]]}"/>
   <position name="posArapucaSub7"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[7]]}"/>
   <position name="posArapucaSub8"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[8]]}"/>
   <position name="posArapucaSub9"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[9]]}"/>
   <position name="posArapucaSub10"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[10]]}"/>
   <position name="posArapucaSub11"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[11]]}"/>
   <position name="posArapucaSub12"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[12]]}"/>
   <position name="posArapucaSub13"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[13]]}"/>
   <position name="posArapucaSub14"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[14]]}"/>
   <position name="posArapucaSub15"    unit="mm" x="@{[$pos_subtraction_arapuca_x]}" y="0" z="@{[$list_pos[15]]}"/>
   <position name="posPaddleSub0"    unit="mm" x="0" y="@{[$list_sub_paddle[0]]}" z="0"/>
   <position name="posPaddleSub1"    unit="mm" x="0" y="@{[$list_sub_paddle[1]]}" z="0"/>
   <position name="posPaddleSub2"    unit="mm" x="0" y="@{[$list_sub_paddle[2]]}" z="0"/>
   <position name="posPaddleSub3"    unit="mm" x="0" y="@{[$list_sub_paddle[3]]}" z="0"/>
   <position name="posPaddleSub4"    unit="mm" x="0" y="@{[$list_sub_paddle[4]]}" z="0"/>
   <position name="posPaddleSub5"    unit="mm" x="0" y="@{[$list_sub_paddle[5]]}" z="0"/>
   <position name="posPaddleSub6"    unit="mm" x="0" y="@{[$list_sub_paddle[6]]}" z="0"/>
   <position name="posPaddleSub7"    unit="mm" x="0" y="@{[$list_sub_paddle[7]]}" z="0"/>
   <position name="posPaddleSub8"    unit="mm" x="0" y="@{[$list_sub_paddle[8]]}" z="0"/>
   <position name="posPaddleSub9"    unit="mm" x="0" y="@{[$list_sub_paddle[9]]}" z="0"/>
   <position name="posSubRailIU"    unit="cm" x="0" y="0" z="-4.80425"/>
   <position name="posSubRailMIT"    unit="cm" x="0" y="0" z="-2.87875"/>

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
    $MAT = "subfile_protodune_v7_Materials" . $suffix . ".gdml";
    push (@gdmlFiles, $MAT);
    $MAT = ">" . $MAT;

    open(MAT) or die("Could not open file $MAT for writing");

    # Add any materials special to this geometry by defining a mulitline string
    # and passing it to the gdmlMaterials::gen_Materials() function.
    $spdensity = 0.001205*(1-$FracVolOfSteel) + 7.9300*$FracVolOfSteel;
    $ssairfrac = (1-$FracMassOfSteel);
my $asmix = <<EOF;
  <!-- preliminary values -->
  <material name="AirSteelMixture" formula="AirSteelMixture">
   <D value="$spdensity" unit="g/cm3"/>
   <fraction n="$FracMassOfSteel" ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
   <fraction n="$ssairfrac"   ref="Air"/>
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

    my $TPCActive_x   =  $_[0] - 3*$APAWirePlaneSpacing -3*$TPCWirePlaneThickness;
    if( ($protoDune == 1) && ($_[3] eq 'Outer') ){
       $TPCActive_x   =  3*($APAWirePlaneSpacing + $TPCWirePlaneThickness);
    }
    my $TPCActive_y   =  $_[1] - $APAGap_y/2 - $ReadoutBoardOverlap + $G10thickness; #TODO: make the Active height more accurate
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
    $TPC = "subfile_protodune_v7_TPC_${_[3]}" . $suffix . ".gdml";
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
     $uplact_y = $Uactive_y + $UVPlaneBoundNudge;
     $uplact_z = $Uactive_z + $UVPlaneBoundNudge;
     $vplact_y = $Vactive_y + $UVPlaneBoundNudge;
     $vplact_z = $Vactive_z + $UVPlaneBoundNudge;
print TPC <<EOF;
<solids>
    <box name="$_[3]" lunit="cm" 
      x="$_[0]" 
      y="$_[1]" 
      z="$_[2]"/>
    <box name="${_[3]}UPlane" lunit="cm" 
      x="$TPCWirePlaneThickness" 
      y="$uplact_y" 
      z="$uplact_z"/>
    <box name="${_[3]}VPlane" lunit="cm" 
      x="$TPCWirePlaneThickness" 
      y="$vplact_y" 
      z="$vplact_z"/>
    <box name="${_[3]}ZPlane" lunit="cm" 
      x="$TPCWirePlaneThickness" 
      y="$Zactive_y" 
      z="$Zactive_z"/>
EOF
#if(${_[3]} eq "Inner"){
#print TPC <<EOF;
#    <box name="${_[3]}Active_NoCPA" lunit="cm"
#      x="$TPCActive_x"
#      y="$TPCActive_y"
#      z="$TPCActive_z"/>
#EOF
#}
#else{
print TPC <<EOF;
    <box name="${_[3]}Active" lunit="cm"
      x="$TPCActive_x"
      y="$TPCActive_y"
      z="$TPCActive_z"/>
EOF
#}


#++++++++++++++++++++++++++++ Wire Solids ++++++++++++++++++++++++++++++
$wirerad = 0.5*$TPCWireThickness;
print TPC <<EOF;

    <tube name="${_[3]}WireVert"
      rmax="$wirerad"
      z="$Zactive_y"               
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
      rmax="$wirerad"
      z="$CornerUWireLength[$i]"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

    }

    $CommonUWireLength = $Uactive_z/$SinUAngle;

   print TPC <<EOF;
    <tube name="${_[3]}WireUCommon"
      rmax="$wirerad"
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
      rmax="$wirerad"
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
        $CornerVWireLength[$i] = ($FirstVWireOffset + $i*$VWire_zint)/$SinVAngle;

        print TPC <<EOF;

    <tube name="${_[3]}WireV$i"
      rmax="$wirerad"
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
      rmax="$wirerad"
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
      rmax="$wirerad"
      z="$TopCornerVWireLength[$i]"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>
EOF

    }


}

# The active inner and active outer volumes are treated slightly differently 
my $subname = '';
my $fieldStrength = '500';
if ( ${_[3]} ne 'Inner' ){ 
  $fieldStrength = '0'; #assumption: no field in outer volume. This is currently unused in larg4
  $subname = ${_[3]}; #The volume name has to be "volTPCActive" for Inner
                    #Can be anything for the outer volumes
}

# Begin structure and create the vertical wire logical volume
# D.R. Add auxtypes utilized by larg4
print TPC <<EOF;
</solids>
<structure>
    <volume name="volTPCActive${subname}">
      <materialref ref="LAr"/>
      <solidref ref="${_[3]}Active"/>
      <auxiliary auxtype="SensDet" auxvalue="SimEnergyDeposit"/>
      <auxiliary auxtype="StepLimit" auxunit="cm" auxvalue="0.47625*cm"/>
      <auxiliary auxtype="Efield" auxunit="V/cm" auxvalue="${fieldStrength}*V/cm"/>
      <colorref ref="blue"/>
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

$xrotUAang = 90-$UAngle;
print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireU$i${_[3]}"/>
        <position name="pos${_[3]}WireU$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotation name="rUAngle$i"   unit="deg" x="$xrotUAang" y="0"   z="0"/>
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

$xrotUAang = 90-$UAngle;
print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireUCommon${_[3]}"/>
        <position name="pos${_[3]}WireU$i" unit="cm" x="0" y="$ypos " z="0"/>
        <rotation name="rUAngle$i"   unit="deg" x="$xrotUAang" y="0"   z="0"/>
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

$xrotUAang = 90-$UAngle;
print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireU$j${_[3]}"/>
        <position name="pos${_[3]}WireU$j" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotation name="rUAngle$j"   unit="deg" x="$xrotUAang" y="0"   z="0"/>
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

$xrotVAang = 90+$VAngle; 
print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireV$i${_[3]}"/>
        <position name="pos${_[3]}WireV$i" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotation name="rVAngle$i"   unit="deg" x="$xrotVAang" y="0"   z="0"/>
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

$xrotVAang = 90+$VAngle;
print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireVCommon${_[3]}"/>
        <position name="pos${_[3]}WireV$i" unit="cm" x="0" y="$ypos " z="0"/>
        <rotation name="rVAngle$i"   unit="deg" x="$xrotVAang" y="0"   z="0"/>
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

$xrotVAang = 90+$VAngle;
print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireV$j${_[3]}"/>
        <position name="pos${_[3]}WireV$j" unit="cm" x="0" y="$ypos " z="$zpos"/>
        <rotation name="rVAngle$j"   unit="deg" x="$xrotVAang" y="0"   z="0"/>
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

    my $BottomOfAPA = - $TPC_y/2 + $APAGap_y/2;


    $posZplane[0]   = -$_[0]/2 + $APAWirePlaneSpacing + $TPCWirePlaneThickness/2;
    $posZplane[1]   = $BottomOfAPA + $WrapCover + 4*$G10thickness + $Zactive_y/2;
    $posZplane[2]   = 0;

    $posVplane[0]   = $posZplane[0] + $APAWirePlaneSpacing + $TPCWirePlaneThickness;
    $posVplane[1]   = $BottomOfAPA + $WrapCover + 3*$G10thickness + $Vactive_y/2;
    $posVplane[2]   = $posZplane[2];

    $posUplane[0]   = $posVplane[0] + $APAWirePlaneSpacing + $TPCWirePlaneThickness;
    $posUplane[1]   = $BottomOfAPA + $WrapCover + 2*$G10thickness + $Uactive_y/2;
    $posUplane[2]   = $posZplane[2];

    $posTPCActive[0] = $posUplane[0] + $TPCWirePlaneThickness/2 + $TPCActive_x/2;
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
         x="$posZplane[0]" y="$posZplane[1]" z="$posZplane[2]"/>
       <rotationref ref="rIdentity"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneV${_[3]}"/>
       <position name="pos${_[3]}PlaneV" unit="cm" 
         x="$posVplane[0]" y="$posVplane[1]" z="$posVplane[2]"/>
       <rotationref ref="rIdentity"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneU${_[3]}"/>
       <position name="pos${_[3]}PlaneU" unit="cm" 
         x="$posUplane[0]" y="$posUplane[1]" z="$posUplane[2]"/>
       <rotationref ref="rIdentity"/>
     </physvol>
EOF

# Option to omit the active volume on the outside for the prototype.
if( !($killOuterActive==1 && $_[3] eq 'Outer') ){
print TPC <<EOF;
     <physvol>
       <volumeref ref="volTPCActive${subname}"/>
       <position name="pos${_[3]}Active" unit="cm" 
         x="$posTPCActive[0]" y="$posTPCActive[1]" z="$posTPCActive[2]"/>
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
    $CRYO = "subfile_protodune_v7_Cryostat" . $suffix . ".gdml";
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
$APAFrameYCenter_x = 2*$inch;
$APAFrameYCenter_z = 3*$inch;

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

#Cathode dimensions
$CathFrame_x     = 2.5*$inch;
$CathFrame_y     = 237.51*$inch;
$CathFrame_z     = 91.14*$inch;

$CathTop_x   = 1.185*$inch;
$CathTop_y   = 76.42*$inch;
$CathTop_z   = 43.47*$inch;

$CathMid_x   = 1.185*$inch;
$CathMid_y   = 76.67*$inch;
$CathMid_z   = 43.47*$inch;

$CathPos_z   =  $CathFrame_z/2 - $CathTop_z/2 - 1*$inch;
$CathPos_x   =  $CathFrame_x/2 - $CathTop_x/2 + 0.001;  # 0.001 is for clean cut
$CathTopPos_y=  $CathFrame_y/2 - $CathTop_y/2 - 3*$inch;
$CathBotPos_y= -$CathFrame_y/2 + $CathTop_y/2 + 3*$inch;
$CathMidPos_y=  0.; 
print "CathMidPos_y = $CathMidPos_y"."\n";

# All the cryostat solids.
print CRYO <<EOF;
<solids>
    <box name="Cryostat" lunit="cm" 
      x="$Cryostat_x" 
      y="$Cryostat_y" 
      z="$Cryostat_z"/>

    <cutTube name="BeamPlugMod" rmin="0" rmax="12.0" z="$BeamPlugLe" startphi="0" deltaphi="360" lowX="-0.278696241193" lowY="0" lowZ="-0.960379302747" highX="0.278696241193" highY="0" highZ="0.960379302747" aunit="deg" lunit="cm"/>
  
    <cutTube name="BeamPlug_0" rmin="0" rmax="$BeamPlugRad" z="$BeamPlugLe" startphi="0" deltaphi="360" lowX="-0.278696241193" lowY="0" lowZ="-0.960379302747" highX="0.278696241193" highY="0" highZ="0.960379302747" aunit="deg" lunit="cm"/>

    <cutTube name="BeamPlugNi_0" rmin="0" rmax="$BeamPlugNiRad" z="$BeamPlugNiLe" startphi="0" deltaphi="360" lowX="-0.278696241193" lowY="0" lowZ="-0.960379302747" highX="0.278696241193" highY="0" highZ="0.960379302747" aunit="deg" lunit="cm"/>

    <torus   name="TorRing" rmin="0" rmax="1.7" rtor="10.1" deltaphi="360." startphi="0" aunit="deg" lunit="cm"/>
    <tube   name="TubeRing1" rmin="$BeamPlugNiRad" rmax="$BeamPlugRad" z="0.5" deltaphi="360." startphi="0" aunit="deg" lunit="cm"/>
    <tube   name="TubeRing2" rmin="$BeamPlugNiRad" rmax="$BeamPlugRad" z="3.4" deltaphi="360." startphi="0" aunit="deg" lunit="cm"/>

    <tube   name="TubeBePlFlangeOut" rmin="0" rmax="15" z="3.6" deltaphi="360." startphi="0" aunit="deg" lunit="cm"/>
    <tube   name="TubeBePlFlangeIn" rmin="0" rmax="8.25" z="3.2" deltaphi="360." startphi="0" aunit="deg" lunit="cm"/>
    <tube   name="TubeBePlFlangeCut" rmin="12.5" rmax="15.01" z="2.302" deltaphi="360." startphi="0" aunit="deg" lunit="cm"/>
    <tube   name="TubeBePlFlangeTPCRem" rmin="0." rmax="12.5" z="0.91" deltaphi="360." startphi="0" aunit="deg" lunit="cm"/>

    <subtraction name="BPRing">
      <first ref="TorRing"/>
      <second ref="TubeRing2"/>
    </subtraction>
    <union name="BeamPlugERing">
      <first ref="BPRing"/>
      <second ref="TubeRing1"/>
    </union>

    <subtraction name="BePlFlangeSolid">
      <first ref="TubeBePlFlangeOut"/>
      <second ref="TubeBePlFlangeCut"/>
      <position name="posSubBePlFlangeSolid" x="0" y="0" z="0.651" unit="cm"/>
    </subtraction>
    <subtraction name="BePlFlange">
      <first ref="BePlFlangeSolid"/>
      <second ref="TubeBePlFlangeIn"/>
      <position name="posSubBePlFlange" x="0" y="0" z="-0.20" unit="cm"/>
    </subtraction>

EOF

for($i=0 ; $i<7 ; $i++){

$old_i = $i;
$new_i = $i+1;
$BeamPlugERing_z   = $BeamPlugERingPosZ + $i*5.7;
$BeamPlugNiERing_z = $BeamPlugERingPosZ + $i*5.7 - $BeamPlugNiPos_z;
 
print CRYO <<EOF;
    <subtraction name="BeamPlug_$new_i">
      <first ref="BeamPlug_$old_i"/>
      <second ref="BeamPlugERing"/>
      <position name="posBeamPlugSubtr_$new_i" x="0" y="0" z="$BeamPlugERing_z" unit="cm"/>
    </subtraction>
    <subtraction name="BeamPlugNi_$new_i">
      <first ref="BeamPlugNi_$old_i"/>
      <second ref="BeamPlugERing"/>
      <position name="posBeamPlugNiSubtr_$new_i" x="0" y="0" z="$BeamPlugNiERing_z" unit="cm"/>
    </subtraction>
EOF
}

print CRYO <<EOF;
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

    <box name="Cathode" lunit="cm"
      x="$Cathode_x"
      y="$CathFrame_y"
      z="$CathFrame_z"/>


EOF
#    <box name="CathFrame" lunit="cm"
#      x="$CathFrame_x"
#      y="$CathFrame_y"
#      z="$CathFrame_z"/>

#    <box name="CathTop" lunit="cm"
#      x="$CathTop_x"
#      y="$CathTop_y"
#      z="$CathTop_z"/>

#    <box name="CathMid" lunit="cm"
#      x="$CathMid_x"
#      y="$CathMid_y"
#      z="$CathMid_z"/>

#    <subtraction name="Cath-1">
#      <first ref="CathFrame"/>
#      <second ref="CathTop"/>
#      <position name="posCath-1" x="$CathPos_x" y="$CathTopPos_y" z="$CathPos_z" unit="cm"/>
#    </subtraction>
#    <subtraction name="Cath-2">
#      <first ref="Cath-1"/>
#      <second ref="CathTop"/>
#      <position name="posCath-2" x="-$CathPos_x" y="$CathTopPos_y" z="$CathPos_z" unit="cm"/>
#    </subtraction>
#    <subtraction name="Cath-3">
#      <first ref="Cath-2"/>
#      <second ref="CathTop"/>
#      <position name="posCath-3" x="$CathPos_x" y="$CathTopPos_y" z="-$CathPos_z" unit="cm"/>
#    </subtraction>
#    <subtraction name="Cath-4">
#      <first ref="Cath-3"/>
#      <second ref="CathTop"/>
#      <position name="posCath-4" x="-$CathPos_x" y="$CathTopPos_y" z="-$CathPos_z" unit="cm"/>
#    </subtraction>

#    <subtraction name="Cath-5">
#      <first ref="Cath-4"/>
#      <second ref="CathMid"/>
#      <position name="posCath-5" x="$CathPos_x" y="$CathMidPos_y" z="$CathPos_z" unit="cm"/>
#    </subtraction>
#    <subtraction name="Cath-6">
#      <first ref="Cath-5"/>
#      <second ref="CathMid"/>
#      <position name="posCath-6" x="-$CathPos_x" y="$CathMidPos_y" z="$CathPos_z" unit="cm"/>
#    </subtraction>
#    <subtraction name="Cath-7">
#      <first ref="Cath-6"/>
#      <second ref="CathMid"/>
#      <position name="posCath-7" x="$CathPos_x" y="$CathMidPos_y" z="-$CathPos_z" unit="cm"/>
#    </subtraction>
#    <subtraction name="Cath-8">
#      <first ref="Cath-7"/>
#      <second ref="CathMid"/>
#      <position name="posCath-8" x="-$CathPos_x" y="$CathMidPos_y" z="-$CathPos_z" unit="cm"/>
#    </subtraction>


#    <subtraction name="Cath-9">
#      <first ref="Cath-8"/>
#      <second ref="CathTop"/>
#      <position name="posCath-9" x="$CathPos_x" y="$CathBotPos_y" z="$CathPos_z" unit="cm"/>
#    </subtraction>
#    <subtraction name="Cath-10">
#      <first ref="Cath-9"/>
#      <second ref="CathTop"/>
#      <position name="posCath-10" x="-$CathPos_x" y="$CathBotPos_y" z="$CathPos_z" unit="cm"/>
#    </subtraction>
#    <subtraction name="Cath-11">
#      <first ref="Cath-10"/>
#      <second ref="CathTop"/>
#      <position name="posCath-11" x="$CathPos_x" y="$CathBotPos_y" z="-$CathPos_z" unit="cm"/>
#    </subtraction>
#    <subtraction name="Cath-12">
#      <first ref="Cath-11"/>
#      <second ref="CathTop"/>
#      <position name="posCath-12" x="-$CathPos_x" y="$CathBotPos_y" z="-$CathPos_z" unit="cm"/>
#    </subtraction>
#EOF

#    <subtraction name="InnerActive">
#      <first ref="InnerActive_NoCPA"/>
#      <second ref="Cath-12"/>
#      <position name="posCPAInnerRem" x="179.23175" y="-3.805" z="0." unit="cm"/>
#    </subtraction>
#
#################################
print CRYO <<EOF;
    <box name="boxFCEWmod" x="354.3" y="151.2" z="12.8" lunit="cm"/>
    <box name="boxFCEWHorSuppOut" x="354.3" y="15.2" z="10.2" lunit="cm"/>
    <box name="boxFCEWHorSuppInn" x="354.301" y="14.0" z="9.0" lunit="cm"/>
    <box name="boxFCEWRem" x="354.301" y="15.201" z="1.3" lunit="cm"/>
    <box name="boxFCEWFieldShaping" x="354.3" y="151.2" z="1.3" lunit="cm"/>
    <box name="boxFCEWVerSupp" x="15.2" y="151.2" z="1.3" lunit="cm"/>

    <subtraction name="FCEWHorSupp">
      <first ref="boxFCEWHorSuppOut"/>
      <second ref="boxFCEWHorSuppInn"/>
    </subtraction>
    <subtraction name="FCEWHorSuppBar">
      <first ref="FCEWHorSupp"/>
      <second ref="boxFCEWRem"/>
      <position name="posFCEWHorSuppBar" x="0" y="0" z="2." unit="cm"/>
    </subtraction>

    <box name="boxFCEW-BP-Subtract" x="50.8" y="151.2" z="1.301" lunit="cm"/>
    <tube name="tubBePlHolemod" rmin="0" rmax="15.25" z="12.801" deltaphi="360." startphi="0" aunit="deg" lunit="cm"/>
    <tube name="tubBePlHoleFiSh" rmin="0" rmax="15.25" z="1.301" deltaphi="360." startphi="0" aunit="deg" lunit="cm"/>
    <box name="boxFCEW-BPFla-VerSupp" x="50.8" y="151.2" z="1.3" lunit="cm"/>
    <box name="boxFCEW-BPFla-Subtr" x="25.4" y="64" z="1.301" lunit="cm"/>
    <tube name="tubFCEW-BPFla-Subtr" rmin="0" rmax="14.0" z="1.301" deltaphi="360." startphi="0" aunit="deg" lunit="cm"/>
    <tube name="tubFCEW-BPVer-Subtr" rmin="0" rmax="12.5" z="1.301" deltaphi="360." startphi="0" aunit="deg" lunit="cm"/>

    <subtraction name="boxFCEW-BP-NoHole">
      <first ref="boxFCEWmod"/>
      <second ref="boxFCEW-BP-Subtract"/>
      <position name="posFCEW-BP-NoHole" x="147.8" y="0" z="5.7505" unit="cm"/>
    </subtraction>
    <subtraction name="boxFCEW-BP-mod">
      <first ref="boxFCEW-BP-NoHole"/>
      <second ref="tubBePlHolemod"/>
      <position name="posboxFCEW-BP-mod" x="146" y="41.6" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="FCEW-BP-FieldShaping">
      <first ref="boxFCEWFieldShaping"/>
      <second ref="tubBePlHoleFiSh"/>
      <position name="posFCEW-BP-FieldShapingSol" x="146" y="41.6" z="0" unit="cm"/>
    </subtraction>

    <subtraction name="FCEW-BPFla-Ver">
      <first ref="boxFCEW-BPFla-VerSupp"/>
      <second ref="boxFCEW-BPFla-Subtr"/>
      <position name="posFCEW-BPFla-Ver" x="0." y="-11.6" z="0." unit="cm"/>
    </subtraction>
    <subtraction name="FCEW-BPFla">
      <first ref="FCEW-BPFla-Ver"/>
      <second ref="tubFCEW-BPFla-Subtr"/>
      <position name="posFCEW-BPFlaSolid" x="-1.35" y="41.4" z="0." unit="cm"/>
    </subtraction>

    <subtraction name="FCEW-BPVer">
      <first ref="boxFCEWVerSupp"/>
      <second ref="tubFCEW-BPVer-Subtr"/>
      <position name="posFCEW-BP-Ver" x="-14.9" y="43.3" z="0." unit="cm"/>
    </subtraction>

    <box name="boxFCMod" x="352" y="27.4" z="230" lunit="cm"/>
    <box name="boxShieldPlate" x="264" y="2" z="230" lunit="cm"/>
    <box name="boxFCXBarOut" x="352" y="15.2" z="7.6" lunit="cm"/>
    <box name="boxFCXBarRem" x="352.001" y="13.9" z="3.501" lunit="cm"/>
    <box name="boxFCXFieldBarRem" x="352.001" y="1.3" z="0.61" lunit="cm"/>
    <box name="boxFCZBarOut" x="3.8" y="7.6" z="115.4" lunit="cm"/>
    <box name="boxFCZBarRem" x="1.58" y="6.3" z="115.401" lunit="cm"/>
    <box name="boxFCFieldBarPlane" x="352" y="1.3" z="230" lunit="cm"/>

    <subtraction name="FCXBarOut-1">
      <first ref="boxFCXBarOut"/>
      <second ref="boxFCXBarRem"/>
      <position name="posFCXBarOutSol-1" x="0" y="0" z="2.05" unit="cm"/>
    </subtraction>
    <subtraction name="FCXBarOut-2">
      <first ref="FCXBarOut-1"/>
      <second ref="boxFCXBarRem"/>
      <position name="posFCXBarOutSol-2" x="0" y="0" z="-2.05" unit="cm"/>
    </subtraction>
    <subtraction name="FCXBarOut">
      <first ref="FCXBarOut-2"/>
      <second ref="boxFCXFieldBarRem"/>
      <position name="posFCXBarOut" x="0" y="-3.3" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="FCZBarOut-1">
      <first ref="boxFCZBarOut"/>
      <second ref="boxFCZBarRem"/>
      <position name="posFCZBarOut-1" x="1.111" y="0" z="0" unit="cm"/>
    </subtraction>
    <subtraction name="FCZBar">
      <first ref="FCZBarOut-1"/>
      <second ref="boxFCZBarRem"/>
      <position name="posFCZBar" x="-1.111" y="0" z="0" unit="cm"/>
    </subtraction>


EOF

$LightPaddleIUSiPM_z = $LightPaddleIU_z + $SiPM_z;
$LightPaddleMITSiPM_z = $LightPaddleMIT_z + $SiPM_z;
$LightPaddleSiPM_z = $LightPaddle_z + $SiPM_z;
$APAFrameYSideHollow_x = $APAFrameYSide_x-2*$EdgeFrameSteelThickness;
$APAFrameYSideHollow_y = $APAFrameYSide_y-2*$EdgeFrameSteelThickness;
$APAFrameZSideHollow_x = $APAFrameZSide_x-2*$EdgeFrameSteelThickness;
$APAFrameZSideHollow_y = $APAFrameZSide_y-2*$EdgeFrameSteelThickness;
$APAFrameYCenterHollow_x = $APAFrameYCenter_x-2*$EdgeFrameSteelThickness;
print CRYO <<EOF;
    <box name="LightPaddleIU" lunit="cm"
      x="$LightPaddleIU_x"
      y="$LightPaddleIU_y"
      z="$LightPaddleIUSiPM_z"/>

    <box name="LightPaddleIUOut" lunit="cm"
      x="$LightPaddleIUOut_x"
      y="$LightPaddleIUOut_y"
      z="$LightPaddleSiPM_z"/>
     <subtraction name="PaddleRailIU">
      <first  ref="LightPaddleIUOut"/>
      <second ref="LightPaddleIU"/>
      <positionref ref="posSubRailIU"/>
      </subtraction>

    <box name="LightPaddleMIT" lunit="cm"
      x="$LightPaddleMIT_x"
      y="$LightPaddleMIT_y"
      z="$LightPaddleMITSiPM_z"/>

    <box name="LightPaddleMITOut" lunit="cm"
      x="$LightPaddleMITOut_x"
      y="$LightPaddleMITOut_y"
      z="$LightPaddleSiPM_z"/>
     <subtraction name="PaddleRailMIT">
      <first  ref="LightPaddleMITOut"/>
      <second ref="LightPaddleMIT"/>
      <positionref ref="posSubRailMIT"/>
      </subtraction>

    <box name="ArapucaOut" lunit="mm"
      x="@{[$ArapucaOut_x]}"
      y="@{[$ArapucaOut_y]}"
      z="@{[$ArapucaOut_z]}"/>

    <box name="ArapucaIn" lunit="mm"
      x="@{[$ArapucaIn_x]}"
      y="@{[$ArapucaIn_y]}"
      z="@{[$ArapucaIn_z]}"/>

     <subtraction name="ArapucaWalls0">
      <first  ref="ArapucaOut"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub0"/>
      </subtraction>
     <subtraction name="ArapucaWalls1">
      <first  ref="ArapucaWalls0"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub1"/>
      </subtraction>
     <subtraction name="ArapucaWalls2">
      <first  ref="ArapucaWalls1"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub2"/>
      </subtraction>
     <subtraction name="ArapucaWalls3">
      <first  ref="ArapucaWalls2"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub3"/>
      </subtraction>
     <subtraction name="ArapucaWalls4">
      <first  ref="ArapucaWalls3"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub4"/>
      </subtraction>
     <subtraction name="ArapucaWalls5">
      <first  ref="ArapucaWalls4"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub5"/>
      </subtraction>
     <subtraction name="ArapucaWalls6">
      <first  ref="ArapucaWalls5"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub6"/>
      </subtraction>
     <subtraction name="ArapucaWalls7">
      <first  ref="ArapucaWalls6"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub7"/>
      </subtraction>
     <subtraction name="ArapucaWalls8">
      <first  ref="ArapucaWalls7"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub8"/>
      </subtraction>
     <subtraction name="ArapucaWalls9">
      <first  ref="ArapucaWalls8"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub9"/>
      </subtraction>
     <subtraction name="ArapucaWalls10">
      <first  ref="ArapucaWalls9"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub10"/>
      </subtraction>
     <subtraction name="ArapucaWalls11">
      <first  ref="ArapucaWalls10"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub11"/>
      </subtraction>
     <subtraction name="ArapucaWalls12">
      <first  ref="ArapucaWalls11"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub12"/>
      </subtraction>
     <subtraction name="ArapucaWalls13">
      <first  ref="ArapucaWalls12"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub13"/>
      </subtraction>
     <subtraction name="ArapucaWalls14">
      <first  ref="ArapucaWalls13"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub14"/>
      </subtraction>
     <subtraction name="ArapucaWalls">
      <first  ref="ArapucaWalls14"/>
      <second ref="ArapucaIn"/>
      <positionref ref="posArapucaSub15"/>
      </subtraction>


    <box name="ArapucaAcceptanceWindow" lunit="mm"
      x="@{[$ArapucaAcceptanceWindow_x]}"
      y="@{[$ArapucaAcceptanceWindow_y]}"
      z="@{[$ArapucaAcceptanceWindow_z]}"/>

     <box name="APAFrameYSideHollow" lunit="cm"
      x="$APAFrameYSideHollow_x"
      y="$APAFrameYSideHollow_y"
      z="$APAFrameYSide_z"/>
     <box name="APAFrameYSideShell" lunit="cm"
      x="$APAFrameYSide_x"
      y="$APAFrameYSide_y"
      z="$APAFrameYSide_z"/>

     <subtraction name="APAFrameYSide">
      <first  ref="APAFrameYSideShell"/>
      <second ref="APAFrameYSideHollow"/>
      <positionref ref="posCenter"/>
      <rotationref ref="rIdentity"/>
      </subtraction>

     <box name="APAFrameYCenterHollow" lunit="cm"
      x="$APAFrameYCenterHollow_x"
      y="$APAFrameYSideHollow_y"
      z="$APAFrameYCenter_z"/>
     <box name="APAFrameYCenterShell" lunit="cm"
      x="$APAFrameYCenter_x"
      y="$APAFrameYSide_y"
      z="$APAFrameYCenter_z"/>

     <subtraction name="APAFrameYCenter_v1">
      <first  ref="APAFrameYCenterShell"/>
      <second ref="APAFrameYCenterHollow"/>
      <positionref ref="posCenter"/>
      <rotationref ref="rIdentity"/>
      </subtraction>

    <subtraction name="APAFrameYCenter0">
      <first  ref="APAFrameYCenter_v1"/>
      <second ref="ArapucaOut"/>
      <positionref ref="posPaddleSub0"/>
      </subtraction>
    <subtraction name="APAFrameYCenter1">
      <first  ref="APAFrameYCenter0"/>
      <second ref="ArapucaOut"/>
      <positionref ref="posPaddleSub1"/>
      </subtraction>
    <subtraction name="APAFrameYCenter2">
      <first  ref="APAFrameYCenter1"/>
      <second ref="ArapucaOut"/>
      <positionref ref="posPaddleSub2"/>
      </subtraction>
    <subtraction name="APAFrameYCenter3">
      <first  ref="APAFrameYCenter2"/>
      <second ref="ArapucaOut"/>
      <positionref ref="posPaddleSub3"/>
      </subtraction>
    <subtraction name="APAFrameYCenter4">
      <first  ref="APAFrameYCenter3"/>
      <second ref="ArapucaOut"/>
      <positionref ref="posPaddleSub4"/>
      </subtraction>
    <subtraction name="APAFrameYCenter5">
      <first  ref="APAFrameYCenter4"/>
      <second ref="ArapucaOut"/>
      <positionref ref="posPaddleSub5"/>
      </subtraction>
    <subtraction name="APAFrameYCenter6">
      <first  ref="APAFrameYCenter5"/>
      <second ref="ArapucaOut"/>
      <positionref ref="posPaddleSub6"/>
      </subtraction>
    <subtraction name="APAFrameYCenter7">
      <first  ref="APAFrameYCenter6"/>
      <second ref="ArapucaOut"/>
      <positionref ref="posPaddleSub7"/>
      </subtraction>
    <subtraction name="APAFrameYCenter8">
      <first  ref="APAFrameYCenter7"/>
      <second ref="ArapucaOut"/>
      <positionref ref="posPaddleSub8"/>
      </subtraction>
    <subtraction name="APAFrameYCenter">
      <first  ref="APAFrameYCenter8"/>
      <second ref="ArapucaOut"/>
      <positionref ref="posPaddleSub9"/>
      </subtraction>

     <box name="APAFrameZSideHollow" lunit="cm"
      x="$APAFrameZSideHollow_x"
      y="$APAFrameZSideHollow_y"
      z="$APAFrameZSide_z"/>
     <box name="APAFrameZSideShell" lunit="cm"
      x="$APAFrameZSide_x"
      y="$APAFrameZSide_y"
      z="$APAFrameZSide_z"/>
     <subtraction name="APAFrameZSide">
      <first  ref="APAFrameZSideShell"/>
      <second ref="APAFrameZSideHollow"/>
      <positionref ref="posCenter"/>
      <rotationref ref="rIdentity"/>
      </subtraction>

     <box name="APAFrameYOuterSupport" lunit="cm"
      x="$EdgeFrameSteelThickness"
      y="$APAFrameYOuterSupport_y"
      z="$APAFrameYOuterSupport_z"/>


     <box name="G10BoardYSideCenterSeg" lunit="cm"
      x="$G10BoardYSide_x"
      y="$G10BoardYSide_y"
      z="$G10BoardYSide_z"/>

     <box name="G10BoardZSideCenterSeg" lunit="cm"
      x="$G10BoardZSide_x"
      y="$G10BoardZSide_y"
      z="$G10BoardZSide_z"/>

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

    <volume name="volBeamPlugNi">
      <materialref ref="NiGas1atm80K"/>
      <solidref ref="BeamPlugNi_7"/>
    </volume>
    <volume name="volBeamPlugG10">
      <materialref ref="G10"/>
      <solidref ref="BeamPlug_7"/>
      <physvol>
        <volumeref ref="volBeamPlugNi"/>
        <position name="posBeamPlugNi" x="0" y="0" z="$BeamPlugNiPos_z" unit="cm"/>
      </physvol>      
    </volume>

    <volume name="volBeamPlugERing">
      <materialref ref="ALUMINUM_Al"/>
      <solidref ref="BeamPlugERing"/>
    </volume>

    <volume name="volBeamPlugMod">
      <materialref ref="LAr"/>
      <solidref ref="BeamPlugMod"/>
      <physvol>
        <volumeref ref="volBeamPlugG10"/>
        <positionref ref="posCenter"/>
      </physvol>
EOF

for($i=0 ; $i<7 ; $i++){
$BeamPlugERing_z = $BeamPlugERingPosZ + $i*5.7; 
print CRYO <<EOF;
      <physvol name="volBeamPlugERing_$i">
        <volumeref ref="volBeamPlugERing"/>
        <position name="posBeamPlugERing_$i" x="0" y="0" z="$BeamPlugERing_z" unit="cm"/>
      </physvol>
EOF
}
print CRYO <<EOF;
    </volume>

    <volume name="volBePlFlangeNi">
      <materialref ref="NiGas1atm80K"/>
      <solidref ref="TubeBePlFlangeIn"/>
    </volume>
    <volume name="volBePlFlange">
      <materialref ref="G10"/>
      <solidref ref="BePlFlange"/>
    </volume>
    <volume name="volCathode">
      <materialref ref="G10" />
      <solidref ref="Cathode" />
    </volume>
EOF

for($i=0 ; $i<$nAPAs ; $i++){
for($p=0 ; $p<10 ; $p++){
if($p % 2 == 0){
print CRYO <<EOF;
    <volume name="volOpDetSensitive_$i\-$p">
      <materialref ref="LAr"/>
      <solidref ref="LightPaddleIU"/>
    </volume>
EOF
print CRYO <<EOF;
    <volume name="volPaddleRail_$i\-$p">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="PaddleRailIU"/>
    </volume>
EOF
}else{
print CRYO <<EOF;
    <volume name="volOpDetSensitive_$i\-$p">
      <materialref ref="LAr"/>
      <solidref ref="LightPaddleMIT"/>
    </volume>
EOF
print CRYO <<EOF;
    <volume name="volPaddleRail_$i\-$p">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="PaddleRailMIT"/>
    </volume>
EOF
}
}
}

for($i=0 ; $i<$nAPAs ; $i++){
if($i==$nAraAPA1 || $i==$nAraAPA2){
print CRYO <<EOF;
    <volume name="volArapuca_$i">
      <materialref ref="G10" />
      <solidref ref="ArapucaWalls" />
    </volume>
EOF
for($p=0 ; $p<16 ; $p++){
print CRYO <<EOF;
    <volume name="volOpDetSensitive_Arapuca_$i\-$p">
      <materialref ref="LAr"/>
      <solidref ref="ArapucaAcceptanceWindow"/>
    </volume>
EOF
}
}
}

#### End wall field cage module
print CRYO <<EOF;
    <volume name="volFCEWHorSuppBar">
      <materialref ref="FR4"/>
      <solidref ref="FCEWHorSuppBar"/>
    </volume>
    <volume name="volFCEWFieldBarPlane">
      <materialref ref="ALUMINUM_Al"/>
      <solidref ref="boxFCEWFieldShaping"/>
    </volume>
    <volume name="volFCEWVerSuppBar">
      <materialref ref="FR4"/>
      <solidref ref="boxFCEWVerSupp"/>
    </volume>
    <volume name="volFCEWmod">
      <materialref ref="LAr"/>
      <solidref ref="boxFCEWmod"/>
      <physvol name="volFCEWHorSuppBar-1">
        <volumeref ref="volFCEWHorSuppBar"/>
        <position name="posFCEWHorSuppBar-1" x="0" y="54.5" z="0" unit="cm"/>
      </physvol>
      <physvol name="volFCEWHorSuppBar-2">
        <volumeref ref="volFCEWHorSuppBar"/>
        <position name="posFCEWHorSuppBar-2" x="0" y="-54.5" z="0" unit="cm"/>
      </physvol>
      <physvol>
        <volumeref ref="volFCEWFieldBarPlane"/>
        <position name="posFCEWFieldBarPlane" x="0" y="0" z="2" unit="cm"/>
      </physvol>
      <physvol name="volFCEWVerSuppBar-1">
        <volumeref ref="volFCEWVerSuppBar"/>
        <position name="posFCEWVerSuppBar-1" x="-165.75" y="0" z="5.75" unit="cm"/>
      </physvol>
      <physvol name="volFCEWVerSuppBar-2">
        <volumeref ref="volFCEWVerSuppBar"/>
        <position name="posFCEWVerSuppBar-2" x="74.2" y="0" z="5.75" unit="cm"/>
      </physvol>
      <physvol name="volFCEWVerSuppBar-3">
        <volumeref ref="volFCEWVerSuppBar"/>
        <position name="posFCEWVerSuppBar-3" x="164.2" y="0" z="5.75" unit="cm"/>
      </physvol>
      <physvol name="volFCEWVerSuppBar-4">
        <volumeref ref="volFCEWVerSuppBar"/>
        <position name="posFCEWVerSuppBar-4" x="-165.75" y="0" z="-5.75" unit="cm"/>
      </physvol>
      <physvol name="volFCEWVerSuppBar-5">
        <volumeref ref="volFCEWVerSuppBar"/>
        <position name="posFCEWVerSuppBar-5" x="74.2" y="0" z="-5.75" unit="cm"/>
      </physvol>
      <physvol name="volFCEWVerSuppBar-6">
        <volumeref ref="volFCEWVerSuppBar"/>
        <position name="posFCEWVerSuppBar-6" x="164.2" y="0" z="-5.75" unit="cm"/>
      </physvol>
    </volume>


EOF
##### Beam plug end wall module

print CRYO <<EOF;

    <volume name="volFCEW-BP-FieldShaping">
      <materialref ref="ALUMINUM_Al"/>
      <solidref ref="FCEW-BP-FieldShaping"/>
    </volume>
    <volume name="volFCEW-BPFla">
      <materialref ref="FR4"/>
      <solidref ref="FCEW-BPFla"/>
    </volume>
    <volume name="volFCEW-BPVer">
      <materialref ref="FR4"/>
      <solidref ref="FCEW-BPVer"/>
    </volume>

    <volume name="volFCEW-BP-mod">
      <materialref ref="LAr"/>
      <solidref ref="boxFCEW-BP-mod"/>
      <physvol name="volFCEW-BP-HorSuppBar-1">
        <volumeref ref="volFCEWHorSuppBar"/>
        <position name="posFCEW-BP-HorSuppBar-1" x="0" y="64.5" z="0" unit="cm"/>
      </physvol>
      <physvol name="volFCEW-BP-HorSuppBar-2">
        <volumeref ref="volFCEWHorSuppBar"/>
        <position name="posFCEW-BP-HorSuppBar-2" x="0" y="-54.5" z="0" unit="cm"/>
      </physvol>
      <physvol>
        <volumeref ref="volFCEW-BP-FieldShaping"/>
        <position name="posFCEW-BP-FieldShaping" x="0" y="0" z="2" unit="cm"/>
      </physvol>
      <physvol name="volFCEW-BP-VerSuppBar-1">
        <volumeref ref="volFCEWVerSuppBar"/>
        <position name="posFCEW-BP-VerSuppBar-1" x="-74.2" y="0" z="5.75" unit="cm"/>
      </physvol>
      <physvol name="volFCEW-BP-VerSuppBar-2">
        <volumeref ref="volFCEWVerSuppBar"/>
        <position name="posFCEW-BP-VerSuppBar-2" x="-164.2" y="0" z="5.75" unit="cm"/>
      </physvol>
      <physvol name="volFCEW-BP-VerSuppBar-3">
        <volumeref ref="volFCEWVerSuppBar"/>
        <position name="posFCEW-BP-VerSuppBar-3" x="-74.2" y="0" z="-5.75" unit="cm"/>
      </physvol>
      <physvol name="volFCEW-BP-VerSuppBar-4">
        <volumeref ref="volFCEWVerSuppBar"/>
        <position name="posFCEW-BP-VerSuppBar-4" x="-164.2" y="0" z="-5.75" unit="cm"/>
      </physvol>
      <physvol name="volFCEW-BP-VerSuppBar-5">
        <volumeref ref="volFCEWVerSuppBar"/>
        <position name="posFCEW-BP-VerSuppBar-5" x="122.5" y="0" z="-5.75" unit="cm"/>
      </physvol>
      <physvol>
        <volumeref ref="volFCEW-BPVer"/>
        <position name="posFCEW-BPVer" x="169.4" y="0" z="-5.75" unit="cm"/>
      </physvol>
    </volume>

EOF

##### Top and bottom field cage module volume

print CRYO <<EOF;

    <volume name="volFCXBarOut">
      <materialref ref="FR4"/>
      <solidref ref="FCXBarOut"/>
    </volume>
    <volume name="volFCZBar">
      <materialref ref="FR4"/>
      <solidref ref="FCZBar"/>
    </volume>
    <volume name="volFCFieldBarPlane">
      <materialref ref="ALUMINUM_Al"/>
      <solidref ref="boxFCFieldBarPlane"/>
    </volume>
    <volume name="volShieldPlate">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="boxShieldPlate"/>
    </volume>
    <volume name="volFCMod">
      <materialref ref="LAr"/>
      <solidref ref="boxFCMod"/>
      <physvol name="volFCXBarOut-1">
        <volumeref ref="volFCXBarOut"/>
        <position name="posFCXBarOut-1" x="0" y="-6.1" z="58" unit="cm"/>
      </physvol>
      <physvol name="volFCXBarOut-2">
        <volumeref ref="volFCXBarOut"/>
        <position name="posFCXBarOut-2" x="0" y="-6.1" z="-58" unit="cm"/>
      </physvol>
      <physvol>
        <volumeref ref="volFCFieldBarPlane"/>
        <position name="posFCFieldBarPlane" x="0" y="-9.4" z="0" unit="cm"/>
      </physvol>
      <physvol name="volFCZBar-1">
        <volumeref ref="volFCZBar"/>
        <position name="posFCZBar-1" x="141.6" y="-4.1" z="0" unit="cm"/>
      </physvol>
      <physvol name="volFCZBar-2">
        <volumeref ref="volFCZBar"/>
        <position name="posFCZBar-2" x="-41" y="-4.1" z="0" unit="cm"/>
      </physvol>
      <physvol name="volFCZBar-3">
        <volumeref ref="volFCZBar"/>
        <position name="posFCZBar-3" x="-143.5" y="-4.1" z="0" unit="cm"/>
      </physvol>
      <physvol>
        <volumeref ref="volShieldPlate"/>
        <position name="posShieldPlate" x="44" y="12.7" z="0" unit="cm"/>
      </physvol>
    </volume>
EOF


$posGaseousArgon_y = $Argon_y/2-$HeightGaseousAr/2;
$BePlFlangeNi_z    = $BePlFlange_z - 0.2;
print CRYO <<EOF;

    <volume name="volAPAFrameYSide">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="APAFrameYSide"/>
    </volume>

    <volume name="volAPAFrameYCenter">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="APAFrameYCenter"/>
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
        <position name="posGaseousArgon" unit="cm" x="0" y="$posGaseousArgon_y" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volSteelShell"/>
        <position name="posSteelShell" unit="cm" x="0" y="0" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volBeamPlugMod"/>
        <positionref ref="posBeamWAr3"/>       
        <rotationref ref="rBeamWRev3"/> 
      </physvol>
      <physvol>
        <volumeref ref="volBePlFlange"/>
        <position name="posBePlFlange" unit="cm" x="$BePlFlange_x" y="$BePlFlange_y" z="$BePlFlange_z"/>       
      </physvol>
      <physvol>
         <volumeref ref="volBePlFlangeNi"/>
         <position name="posBePlFlangeNi" x="$BePlFlange_x" y="$BePlFlange_y" z="$BePlFlangeNi_z" unit="cm"/>
      </physvol>

EOF
print "Beam Plug position x=".$BePlFlange_x." y=".$BePlFlange_y." z=".$BePlFlange_z."\n";
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
    {     # $j=0 for bottom APAs
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
#                print "$i".", $nAPAWide".", CPA_x = $CPA_x"."\n";

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


$posvolTPC0_x = $APACenter_x - $APAFrame_x/2 - $TPC_0_x/2;
$posvolTPC1_x = $APACenter_x + $APAFrame_x/2 + $TPC_1_x/2;
print CRYO <<EOF;

      <physvol>
        <volumeref ref="volTPC$SelectTPC_0"/>
        <position name="posTPC\-$tpc_0" unit="cm" 
          x="$posvolTPC0_x" 
          y="$APACenter_y" 
          z="$APACenter_z"/> 
        <rotationref ref="$rot_0"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPC$SelectTPC_1"/>
        <position name="posTPC\-$tpc_1" unit="cm" 
          x="$posvolTPC1_x" 
          y="$APACenter_y" 
          z="$APACenter_z"/> 
        <rotationref ref="$rot_1"/>
      </physvol>

EOF


if( $i < $nAPAWide - 1 ){ # avoid placeing the last row of CPAs since the APAs are on the outside

print "$cpa_i".". CPA_x = $CPA_x".", APACenter_y = $APACenter_y".", APACenter_z = $APACenter_z"."\n";
$CPAmodbarGap = 1.8*$inch;
$FCmodRight_x[$cpa_i]  = $CPA_x + 3 + 352/2; # CPAmodCenter_x to FCmodCenter_x 
$FCmodLeft_x[$cpa_i]   = $CPA_x - 3 - 352/2;
$FCmodTop_y[$cpa_i]    = $APACenter_y + 603.3/2 + $CPAmodbarGap + 9.4 + 1.83; # CPAmodCenter_y to FCmodCenter_y with 1.5 inch distance between CPA and module bars
$FCmodBottom_y[$cpa_i] = $APACenter_y - 603.3/2 - $CPAmodbarGap - 9.4 - 1.83;
$FCmod_z[$cpa_i]       = $APACenter_z;

$FCEWmodRight_x = $CPA_x + 3 + 354.3/2;
$FCEWmodLeft_x = $CPA_x - 3 - 354.3/2;
if($cpa_i == 0) {
  $FCEWmodUS_z = $APACenter_z - 231.5/2 - $CPAmodbarGap - 2; #CPAmodCenter_z to FCEWmod 
  for($fcewy_i = 0; $fcewy_i < 4; $fcewy_i++){

    $FCEWmod_y[$fcewy_i] = $APACenter_y + (-3 + 2*$fcewy_i)*77; #distance between FCEW module is 2.8 cm (originally 0.8 cm)

  }
}
if($cpa_i == 2) {
  $FCEWmodDS_z = $APACenter_z + 231.5/2 + $CPAmodbarGap + 2;
}

print CRYO <<EOF;
      <physvol>
        <volumeref ref="volCathode"/>
        <position name="posCathode\-$cpa_i" unit="cm" 
        x="$CPA_x" 
        y="$APACenter_y" 
        z="$APACenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol name="volFCModTopRight\-$cpa_i">
        <volumeref ref="volFCMod"/>
        <position name="posFCModTopRight\-$cpa_i" unit="cm"
         x="$FCmodRight_x[$cpa_i]"
         y="$FCmodTop_y[$cpa_i]"
         z="$FCmod_z[$cpa_i]"/>
        <rotationref ref="rPlus180AboutY"/>
      </physvol>
      <physvol name="volFCModTopLeft\-$cpa_i">
        <volumeref ref="volFCMod"/>
        <position name="posFCModTopLeft\-$cpa_i" unit="cm"
         x="$FCmodLeft_x[$cpa_i]"
         y="$FCmodTop_y[$cpa_i]"
         z="$FCmod_z[$cpa_i]"/>
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol name="volFCModBottomRight\-$cpa_i">
        <volumeref ref="volFCMod"/>
        <position name="posFCModBottomRight\-$cpa_i" unit="cm"
         x="$FCmodRight_x[$cpa_i]"
         y="$FCmodBottom_y[$cpa_i]"
         z="$FCmod_z[$cpa_i]"/>
        <rotationref ref="rPlus180AboutXPlus180AboutY"/>
      </physvol>
      <physvol name="volFCModBottomLeft\-$cpa_i">
        <volumeref ref="volFCMod"/>
        <position name="posFCModBottomLeft\-$cpa_i" unit="cm"
         x="$FCmodLeft_x[$cpa_i]"
         y="$FCmodBottom_y[$cpa_i]"
         z="$FCmod_z[$cpa_i]"/>
        <rotationref ref="rPlus180AboutX"/>
      </physvol>



EOF
          $cpa_i++;
}

################################## End Wall Field Cage Placement


            } #wide



        } # end if not the smaller workspace


        # Make the workspace have only one center APA with CPAs and the full drift on either side
        elsif( $workspace == 1 ){

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

$posvolTPCInner0_x = $APACenter_x - $APAFrame_x/2 - $TPCInner_x/2;              
$posvolTPCInner1_x = $APACenter_x + $APAFrame_x/2 + $TPCInner_x/2;
print CRYO <<EOF;

      <physvol>
        <volumeref ref="volTPCInner"/>
        <position name="posTPC\-$tpc_0" unit="cm" 
          x="$posvolTPCInner0_x" 
          y="$APACenter_y" 
          z="$APACenter_z"/> 
        <rotationref ref="$rot_0"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPCInner"/>
        <position name="posTPC\-$tpc_1" unit="cm"
          x="$posvolTPCInner1_x" 
          y="$APACenter_y"
          z="$APACenter_z"/>
        <rotationref ref="$rot_1"/>
      </physvol>

      <physvol>
        <volumeref ref="volCathode"/>
        <position name="posCathode\-$apa_i-0" unit="cm" 
        x="$CPA_0_x" 
        y="$APACenter_y" 
        z="$APACenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volCathode"/>
        <position name="posCathode\-$apa_i-1" unit="cm" 
        x="$CPA_1_x" 
        y="$APACenter_y" 
        z="$APACenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>

EOF

        } # if the smaller workspace


    } #high
  } #long

}# end if tpc

for($fcewy_i = 0; $fcewy_i < 4; $fcewy_i++){
if($fcewy_i !=2){
print CRYO <<EOF;
      <physvol name="volFCEWModUSLeft\-$fcewy_i">
        <volumeref ref="volFCEWmod"/>
        <position name="posFCEWModUSLeft\-$fcewy_i" unit="cm"
         x="$FCEWmodLeft_x"
         y="$FCEWmod_y[$fcewy_i]"
         z="$FCEWmodUS_z"/>
        <rotationref ref="rPlus180AboutZ"/>
      </physvol>
EOF
}
print CRYO <<EOF;
      <physvol name="volFCEWModUSRight\-$fcewy_i">
        <volumeref ref="volFCEWmod"/>
        <position name="posFCEWModUSRight\-$fcewy_i" unit="cm"
         x="$FCEWmodRight_x"
         y="$FCEWmod_y[$fcewy_i]"
         z="$FCEWmodUS_z"/>
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol name="volFCEWModDSRight\-$fcewy_i">
        <volumeref ref="volFCEWmod"/>
        <position name="posFCEWModDSRight\-$fcewy_i" unit="cm"
         x="$FCEWmodRight_x"
         y="$FCEWmod_y[$fcewy_i]"
         z="$FCEWmodDS_z"/>
        <rotationref ref="rPlus180AboutX"/>
      </physvol>
      <physvol name="volFCEWModDSLeft\-$fcewy_i">
        <volumeref ref="volFCEWmod"/>
        <position name="posFCEWModDSLeft\-$fcewy_i" unit="cm"
         x="$FCEWmodLeft_x"
         y="$FCEWmod_y[$fcewy_i]"
         z="$FCEWmodDS_z"/>
        <rotationref ref="rPlus180AboutXPlus180AboutZ"/>
      </physvol>
EOF
}

$FCEWBPmod_x = $FCEWmodLeft_x + 2.8;
$FCEWBPmod_y = $FCEWmod_y[2]  + 0.9;
$FCEWBPmod_z = $FCEWmodUS_z;

$FCEWBPFla_x = $FCEWBPmod_x + 147.8;
$FCEWBPFla_y = $FCEWBPmod_y  + 0;
$FCEWBPFla_z = $FCEWBPmod_z   + 5.75;

print "volFCEW-BPFla pos x=".$FCEWBPFla_x." y=".$FCEWBPFla_y." z=".$FCEWBPFla_z."\n";

print CRYO <<EOF;
      <physvol>
        <volumeref ref="volFCEW-BP-mod"/>
        <position name="posFCEW-BP-mod" unit="cm"
         x="$FCEWBPmod_x"
         y="$FCEWBPmod_y"
         z="$FCEWBPmod_z"/>
      </physvol>

      <physvol>
        <volumeref ref="volFCEW-BPFla"/>
        <position name="posFCEW-BPFla" unit="cm"
         x="$FCEWBPFla_x"
         y="$FCEWBPFla_y"
         z="$FCEWBPFla_z"/>
      </physvol>
EOF

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

    $test = 0;

for ($paddle = 0; $paddle<$nLightPaddlesPerAPA; $paddle++)
{

             # All Light Paddle centers will have the same
             #  X coordinate as the center of the current APA
             #  Z coordinate as the current TPC pair
             # The Y coordinate must be looped over:

             #the multiplication by j here is a temporary dirty trick to get around some strange behavior

             $Paddle_Y   =    $APACenter_y 
                            - $APAphys_y/2 
                            + $FrameToPaddleSpace 
                            + (1-$j)*($LightPaddleMIT_y/2) 
                            + $PaddleYInterval*$paddle;

     if($apa_i==$nAraAPA1){
        if($paddle==$nSlotAra1) {$test=1;}
    } #where the arapuca bar goes for pDUNE
    if($apa_i==$nAraAPA2){
        if($paddle==$nSlotAra2) {$test=1;}
    }

    if(($test!=1)){
             # Alternate the paddle orientations
             if ( $paddle % 2 == 0 ) { $rot = "rIdentity"; }
             else                    { $rot = "rPlus180AboutY"; }

             $posvolRail_z = $APACenter_z + $SiPM_z/2;

             if ( $paddle % 2 == 0 ) { $posvolOpDetSensitive_z = $APACenter_z + $SiPM_z/2 - 4.80425; } #IU bar
             else                    { $posvolOpDetSensitive_z = $APACenter_z + $SiPM_z/2 + 2.87875; }

             print CRYO <<EOF;
     <physvol>
       <volumeref ref="volOpDetSensitive_$apa_i\-$paddle"/>
       <position name="posPaddle\-$paddle\-TPC\-$i\-$j\-$k" unit="cm" 
         x="$APACenter_x" 
         y="$Paddle_Y" 
         z="$posvolOpDetSensitive_z"/>
       <rotationref ref="$rot"/>
     </physvol>
EOF

             print CRYO <<EOF;
     <physvol>
       <volumeref ref="volPaddleRail_$apa_i\-$paddle"/>
       <position name="posPaddleRail\-$paddle\-TPC\-$i\-$j\-$k" unit="cm" 
         x="$APACenter_x" 
	 y="$Paddle_Y" 
         z="$posvolRail_z"/>
       <rotationref ref="$rot"/>
     </physvol>
EOF

$test = 0;
    }else{
        if ( $apa_i % 2 == 0 ) { 
            $rot = "rIdentity"; 
            $posArax = ($APACenter_x+0.05*$ArapucaOut_x-0.05*$ArapucaAcceptanceWindow_x-0.01);
        }else{ 
            $rot = "rPlus180AboutY"; 
            $posArax = ($APACenter_x-0.05*$ArapucaOut_x+0.05*$ArapucaAcceptanceWindow_x+0.01);
        }
        
        print CRYO <<EOF;
     <physvol>
       <volumeref ref="volArapuca_$apa_i"/>
       <position name="posArapuca$apa_i-TPC\-$i\-$j\-$k" unit="cm" 
         x="@{[$APACenter_x]}"
         y="@{[$Paddle_Y]}" 
         z="@{[$APACenter_z]}"/>
       <rotationref ref="$rot"/>
     </physvol>
EOF
        for ($ara = 0; $ara<16; $ara++){
            print CRYO <<EOF;
     <physvol>
       <volumeref ref="volOpDetSensitive_Arapuca_$apa_i\-$ara"/>
       <position name="posOpArapuca$apa_i\-$ara\-TPC\-$i\-$j\-$k" unit="cm" 
         x="@{[$posArax]}"
         y="@{[$Paddle_Y]}" 
         z="@{[(0.1*$list_pos[$ara]+$APACenter_z)]}"/>
     </physvol>
EOF
}
$test = 0;
}
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
#       x="$APAFrameCenter_x - ($APAFrameYOuterSupport_x + $APAFrameYInnerSupport_x/2 - $EdgeFrameSteelThickness/2)" 
#       y="$APAFrameCenter_y" 
#       z="$APAFrameCenter_z"/> 
#        <rotationref ref="rIdentity"/>
#      </physvol>
#      <physvol>
#        <volumeref ref="volAPAFrameYOuterSupport\-$APA_i"/>
#        <position name="posAPAFrameYOuterSupportPos\-$APA_i" unit="cm" 
#       x="$APAFrameCenter_x + ($APAFrameYOuterSupport_x + $APAFrameYInnerSupport_x/2 - $EdgeFrameSteelThickness/2)" 
#       y="$APAFrameCenter_y" 
#       z="$APAFrameCenter_z"/> 
#        <rotationref ref="rIdentity"/>
#      </physvol>
#-->
#EOF

  $posvolAPAFrameYSideNeg_z = $APAFrameCenter_z - $APAFrame_z/2 + $APAFrameYSide_z/2;
  $posvolAPAFrameYSidePos_z = $APAFrameCenter_z + $APAFrame_z/2 - $APAFrameYSide_z/2;
  $posvolAPAFrameZSidePos_y = $APAFrameCenter_y + $APAFrame_y/2 - $APAFrameZSide_y/2;
  $posvolAPAFrameZSideNeg_y = $APAFrameCenter_y - $APAFrame_y/2 + $APAFrameZSide_y/2;
  print CRYO <<EOF;
      <physvol>
        <volumeref ref="volAPAFrameYSide"/>
        <position name="posAPAFrameYSideNeg\-$APA_i" unit="cm" 
        x="$APAFrameCenter_x" 
        y="$APAFrameCenter_y" 
        z="$posvolAPAFrameYSideNeg_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameYSide"/>
        <position name="posAPAFrameYSidePos\-$APA_i" unit="cm" 
        x="$APAFrameCenter_x" 
        y="$APAFrameCenter_y" 
        z="$posvolAPAFrameYSidePos_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
     <physvol>
       <volumeref ref="volAPAFrameYCenter"/>
       <position name="posFrameYCenter$apa_i" unit="cm" 
         x="$APAFrameCenter_x"
	 y="$APAFrameCenter_y" 
	 z="$APAFrameCenter_z"/>
     </physvol>
      <physvol>
        <volumeref ref="volAPAFrameZSide"/>
        <position name="posAPAFrameZSidePos\-$APA_i" unit="cm" 
        x="$APAFrameCenter_x" 
        y="$posvolAPAFrameZSidePos_y" 
        z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volAPAFrameZSide"/>
        <position name="posAPAFrameZSideNeg\-$APA_i" unit="cm" 
        x="$APAFrameCenter_x" 
        y="$posvolAPAFrameZSideNeg_y" 
        z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>


      <!-- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        - Add the *parts* of the G10 boards that exist directly in volCryostat.
        - There are two boards on each the up and downstream end, 
            one each to wrap the U and V views around the APA frame
        - There are 4 on the bottom which anchor the U V and Z wires and the grid plane
        - The rest of the parts of the G10 boards must be placed directly in volTPC   -->
EOF

  $posvolG10BoardYSideCenterSegVup_z = $APAFrameCenter_z - $APAFrame_z/2 - (0+.5)*($G10BoardYSide_z); 
  $posvolG10BoardYSideCenterSegUup_z = $APAFrameCenter_z - $APAFrame_z/2 - (1+.5)*($G10BoardYSide_z);
  $posvolG10BoardYSideCenterSegVdo_z = $APAFrameCenter_z + $APAFrame_z/2 + (0+.5)*($G10BoardYSide_z);
  $posvolG10BoardYSideCenterSegUdo_z = $APAFrameCenter_z + $APAFrame_z/2 + (1+.5)*($G10BoardYSide_z);
  print CRYO <<EOF;
      <physvol>
        <volumeref ref="volG10BoardYSideCenterSeg"/>
        <position name="posG10BoardYSideCenterSeg\-Vup\-$APA_i" unit="cm" 
        x="$APAFrameCenter_x" 
        y="$APAFrameCenter_y" 
        z="$posvolG10BoardYSideCenterSegVup_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardYSideCenterSeg"/>
        <position name="posG10BoardYSideCenterSeg\-Uup\-$APA_i" unit="cm" 
        x="$APAFrameCenter_x" 
        y="$APAFrameCenter_y" 
        z="$posvolG10BoardYSideCenterSegUup_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>

      <physvol>
        <volumeref ref="volG10BoardYSideCenterSeg"/>
        <position name="posG10BoardYSideCenterSeg\-Vdown\-$APA_i" unit="cm" 
        x="$APAFrameCenter_x" 
        y="$APAFrameCenter_y" 
        z="$posvolG10BoardYSideCenterSegVdo_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardYSideCenterSeg"/>
        <position name="posG10BoardYSideCenterSeg\-Udown\-$APA_i" unit="cm" 
        x="$APAFrameCenter_x" 
        y="$APAFrameCenter_y" 
        z="$posvolG10BoardYSideCenterSegUdo_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>

      <physvol>
        <volumeref ref="volG10BoardZSideCenterSeg"/>
        <position name="posG10BoardZSideCenterSeg\-Z\-$APA_i" unit="cm" 
        x="$APAFrameCenter_x" 
        y="$posG10ZSideZ_y" 
        z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardZSideCenterSeg"/>
        <position name="posG10BoardZSideCenterSeg\-V\-$APA_i" unit="cm" 
        x="$APAFrameCenter_x" 
        y="$posG10ZSideV_y" 
        z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardZSideCenterSeg"/>
        <position name="posG10BoardZSideCenterSeg\-U\-$APA_i" unit="cm" 
        x="$APAFrameCenter_x" 
        y="$posG10ZSideU_y" 
        z="$APAFrameCenter_z"/> 
        <rotationref ref="rIdentity"/>
      </physvol>
      <physvol>
        <volumeref ref="volG10BoardZSideCenterSeg"/>
        <position name="posG10BoardZSideCenterSeg\-Grid\-$APA_i" unit="cm" 
        x="$APAFrameCenter_x" 
        y="$posG10ZSideGrid_y"
        z="$APAFrameCenter_z"/> 
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
    $ENCL = "subfile_protodune_v7_DetEnclosure" . $suffix . ".gdml";
    push (@gdmlFiles, $ENCL);
    $ENCL = ">" . $ENCL;
    open(ENCL) or die("Could not open file $ENCL for writing");


# The standard XML prefix and starting the gdml
    print ENCL <<EOF;
<?xml version='1.0'?>
<gdml>
EOF


# All the detector enclosure solids.
$dimFoamPadBlock_x = $Cryostat_x + 2*$FoamPadding;
$dimFoamPadBlock_y = $Cryostat_y + 2*$FoamPadding;
$dimFoamPadBlock_z = $Cryostat_z + 2*$FoamPadding;
$dimSteelSupportBlock_x = $Cryostat_x + 2*$FoamPadding + 2*$SteelSupport_x;
$dimSteelSupportBlock_y = $Cryostat_y + 2*$FoamPadding + 2*$SteelSupport_y;
$dimSteelSupportBlock_z = $Cryostat_z + 2*$FoamPadding + 2*$SteelSupport_z;
$dimSteelPlateBlock_x = $Cryostat_x + 2*$FoamPadding + 2*$SteelPlate;
$dimSteelPlateBlock_y = $Cryostat_y + 2*$FoamPadding + 2*$SteelPlate;
$dimSteelPlateBlock_z = $Cryostat_z + 2*$FoamPadding + 2*$SteelPlate;

$dimPitDepth = 900;
$dimPitWidth = 1650;
$dimStSuToDSPitWall = 160;
$dimPitLength = $DetEncLength/2 + $dimSteelSupportBlock_z/2 + $dimStSuToDSPitWall;
$dimEnclConcrete_y  = $SpaceSteelSupportToWall + $dimPitDepth; 
$DetCentToJuraPit   = 958;

$dimPit_y = $dimPitDepth + 0.1; # try to avoid overlap
$dimPit_z = $dimPitLength + 0.1;   
$posPitInConcrete_x =   $DetCentToJuraPit - $dimPitWidth/2;
$posPitInConcrete_y =   $dimEnclConcrete_y/2 - $dimPitDepth/2 ;
$posPitInConcrete_z = - $DetEncLength/2 + $dimPitLength/2 ;
$posSurrConcrete_y  = - $DetEncHeight/2 + $dimEnclConcrete_y/2 ;

print "FoamPad x=$dimFoamPadBlock_x".", y=$dimFoamPadBlock_y".", z=$dimFoamPadBlock_z"."\n";
print "SteelPlate x=$dimSteelPlateBlock_x".", y=$dimSteelPlateBlock_y".", z=$dimSteelPlateBlock_z"."\n";
print "SteelSupport x=$dimSteelSupportBlock_x".", y=$dimSteelSupportBlock_y".", z=$dimSteelSupportBlock_z"."\n";
print ENCL <<EOF;
<solids>
    <box name="FoamPadBlock" lunit="cm"
      x="$dimFoamPadBlock_x"
      y="$dimFoamPadBlock_y"
      z="$dimFoamPadBlock_z" />

    <box name="SteelPlateBlock" lunit="cm"
      x="$dimSteelPlateBlock_x"
      y="$dimSteelPlateBlock_y"
      z="$dimSteelPlateBlock_z" />

    <subtraction name="SteelPlateNoBW">
      <first ref="SteelPlateBlock"/>
      <second ref="FoamPadBlock"/>
      <positionref ref="posCenter"/>
    </subtraction>

    <subtraction name="FoamPaddingNoBW">
      <first ref="FoamPadBlock"/>
      <second ref="Cryostat"/>
      <positionref ref="posCenter"/>
    </subtraction>
EOF

#    <cutTube name="BeamWindowStPlate" rmin="0" rmax="$BeamPipeRad" z="$BeamWStPlateLe" startphi="0" deltaphi="360" lowX="-0.239669734044" lowY="0" lowZ="-0.970854478582" highX="0.239669734044" highY="0" highZ="0.970854478582" aunit="deg" lunit="cm"/>

#    <cutTube name="BeamWindowFoam" rmin="0" rmax="$BeamPipeRad" z="$BeamWFoLe" startphi="0" deltaphi="360" lowX="-0.239669734044" lowY="0" lowZ="-0.970854478582" highX="0.239669734044" highY="0" highZ="0.970854478582" aunit="deg" lunit="cm"/>

    print ENCL <<EOF;

    <cutTube name="BeamWindowStPlate" rmin="0" rmax="$BeamPipeRad" z="$BeamWStPlateLe" startphi="0" deltaphi="360" lowX="-0.278696241193" lowY="3.41304459712e-17" lowZ="-0.960379302747" highX="0.278696241193" highY="0" highZ="0.960379302747" aunit="deg" lunit="cm"/>

    <cutTube name="BeamWindowFoam" rmin="0" rmax="$BeamPipeRad" z="$BeamWFoLe" startphi="0" deltaphi="360" lowX="-0.278696241193" lowY="3.41304459712e-17" lowZ="-0.960379302747" highX="0.278696241193" highY="0" highZ="0.960379302747"  aunit="deg" lunit="cm"/>

    <cutTube name="BeamWindowFoamRem" rmin="0" rmax="$BeamPipeRad" z="$BeamWFoRemLe" startphi="0" deltaphi="360" lowX="-0.278696241193" lowY="3.41304459712e-17" lowZ="-0.960379302747" highX="0.278696241193" highY="0" highZ="0.960379302747"  aunit="deg" lunit="cm"/>

    <cutTube name="BeamWindowStSu" rmin="0" rmax="$BeamPipeRad" z="$BeamWStSuLe" startphi="0" deltaphi="360" lowX="-0.278696241193" lowY="3.41304459712e-17" lowZ="-0.960379302747" highX="0.278696241193" highY="0" highZ="0.960379302747"  aunit="deg" lunit="cm"/>

    <cutTube name="BeamWindowGlassWool" rmin="0" rmax="$BeamPipeRad" z="$BeamWGlLe" startphi="0" deltaphi="360" lowX="-0.278696241193" lowY="3.41304459712e-17" lowZ="-0.960379302747" highX="0.278696241193" highY="0" highZ="0.960379302747" aunit="deg" lunit="cm"/>

    <cutTube name="BeamPipe" rmin="0" rmax="$BeamPipeRad" z="$BeamPipeLe" startphi="0" deltaphi="360" lowX="-0.278696241193" lowY="3.41304459712e-17" lowZ="-0.960379302747" highX="0.278696241193" highY="0" highZ="0.960379302747" aunit="deg" lunit="cm"/>

    <cutTube name="BeamPipeVacuum" rmin="0" rmax="$BeamVaPipeRad" z="$BeamVaPipeLe" startphi="0" deltaphi="360" lowX="-0.278696241193" lowY="3.41304459712e-17" lowZ="-0.960379302747" highX="0.278696241193" highY="0" highZ="0.960379302747" aunit="deg" lunit="cm"/>

    <subtraction name="SteelPlate">
      <first ref="SteelPlateNoBW"/>
      <second ref="BeamWindowStPlate"/>
      <position name="posBWStPl" x="$BeamWStPlate_x" y="$BeamWStPlate_y" z="$BeamWStPlate_z" unit="cm"/>
      <rotationref ref="rBeamW3"/>
    </subtraction>

    <subtraction name="FoamPadding">
      <first ref="FoamPaddingNoBW"/>
      <second ref="BeamWindowFoamRem"/>
      <position name="posBWFoPa" x="$BeamWFoRem_x" y="$BeamWFoRem_y" z="$BeamWFoRem_z" unit="cm"/>
      <rotationref ref="rBeamW3"/>
    </subtraction>

    <box name="CRTPaddle" lunit="cm" 
      x="$CRTPaddleWidth"
      y="$CRTPaddleHeight"
      z="$CRTPaddleLength"/>

    <box name="CRTModule" lunit="cm" 
      x="$CRTModWidth"
      y="$CRTModHeight"
      z="$CRTModLength"/>

    <box name="DetEnclosure" lunit="cm" 
      x="$DetEncWidth"
      y="$DetEncHeight"
      z="$DetEncLength"/>

    <box name="SurrConcreteBox" lunit="cm" 
      x="$DetEncWidth"
      y="$dimEnclConcrete_y"
      z="$DetEncLength"/>

    <box name="PitBox" lunit="cm" 
      x="$dimPitWidth"
      y="$dimPit_y"
      z="$dimPit_z"/>

    <subtraction name="SurrConcrete">
      <first ref="SurrConcreteBox"/>
      <second ref="PitBox"/>
      <position name="posPitInConcrete" x="$posPitInConcrete_x" y="$posPitInConcrete_y" z="$posPitInConcrete_z" unit="cm"/>
    </subtraction>


EOF
if($simpleStSu == 1)
{
    print ENCL <<EOF;

    <box name="SteelSupportBlock" lunit="cm"
      x="$dimSteelSupportBlock_x"
      y="$dimSteelSupportBlock_y"
      z="$dimSteelSupportBlock_z" />

    <subtraction name="SteelSupportNoBW">
      <first ref="SteelSupportBlock"/>
      <second ref="SteelPlateBlock"/>
      <positionref ref="posCenter"/>
    </subtraction>

    <subtraction name="SteelSupport">
      <first ref="SteelSupportNoBW"/>
      <second ref="BeamWindowStSu"/>
      <position name="posBWStSu" x="$BeamWStSu_x" y="$BeamWStSu_y" z="$BeamWStSu_z" unit="cm"/>
      <rotationref ref="rBeamW3"/>
    </subtraction>

EOF
}
elsif($simpleStSu == 0)
{
    print ENCL <<EOF;

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

    <subtraction name="boxCryoWallSmUS">
      <first ref="boxCryoWallSm"/>
      <second ref="BeamWindowStSu"/>
      <position name="posBWCryoWallUS" x="$BeamWStSu_x" y="$BeamWStSu_y" z="0." unit="cm"/>
      <rotationref ref="rBeamW3"/>
    </subtraction>

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
EOF
}

    print ENCL <<EOF;
</solids>
EOF



# Detector enclosure structure
    print ENCL <<EOF;
<structure>
    <volume name="volBeamWinFoam">
      <materialref ref="ProtoDUNEBWFoam"/>
      <solidref ref="BeamWindowFoam"/>
    </volume>
    <volume name="volBeamWinGlassWool">
      <materialref ref="GlassWool"/>
      <solidref ref="BeamWindowGlassWool"/>
    </volume>
    <volume name="volBeamPipeVac">
      <materialref ref="Vacuum"/>
      <solidref ref="BeamPipeVacuum"/>
    </volume>
    <volume name="volBeamPipe">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="BeamPipe"/>
      <physvol>
         <volumeref ref="volBeamPipeVac"/>
         <positionref ref="posCenter"/> 
      </physvol>
    </volume>
    <volume name="volFoamPadding">
      <materialref ref="ProtoDUNEFoam"/>
      <solidref ref="FoamPadding"/>
    </volume>
    <volume name="volSteelPlate">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="SteelPlate"/>
    </volume>
    <volume name="volSurrConcrete">
      <materialref ref="Concrete"/>
      <solidref ref="SurrConcrete"/>
    </volume>
EOF


    print ENCL <<EOF;
    <volume name="volAuxDetSensitiveCRTPaddle">
      <materialref ref="G4_POLYSTYRENE"/>
      <solidref ref="CRTPaddle"/>
      <auxiliary auxtype="SensDet" auxvalue="AuxDet"/>
      <auxiliary auxtype="Solid" auxvalue="True"/>
    </volume>
EOF


for($imod=0 ; $imod<16 ; $imod++){
    $modnum = $imod + 1;
    print ENCL <<EOF;
    <volume name="volAuxDetCRTModuleU$modnum">
      <materialref ref="Air"/>
      <solidref ref="CRTModule"/>
EOF
for($i=0 ; $i<32 ; $i++){
    $paddle_x1 = - $CRTModWidth/2 + $CRTPaddleWidth*($i + 0.5);
    $paddle_x2 = - $CRTModWidth/2 + $CRTPaddleWidth*($i + 1);
    $paddle_y1 =   $CRTPaddleHeight/2;
    $paddle_y2 = - $CRTPaddleHeight/2;
    $paddle_z = 0;
    $padnum1 = 2*$i+(($modnum-1)*64);
    $padnum2 = (2*$i+1)+(($modnum-1)*64);
    $padnum1id = 2*$i+1;
    $padnum2id = 2*$i+2;
    $paddleid1 = "$padnum1id"."U$modnum";
    $paddleid2 = "$padnum2id"."U$modnum";
    print ENCL <<EOF;
       <physvol  name="CRTPaddle$paddleid1" copynumber="$padnum1">
           <volumeref ref="volAuxDetSensitiveCRTPaddle"/>
           <position name="posCRTPaddleSensitive$paddleid1" unit="cm" 
           x="$paddle_x1" 
           y="$paddle_y1" 
           z="$paddle_z"/> 
          <rotationref ref="rIdentity"/>
       </physvol>
       <physvol  name="CRTPaddle$paddleid2" copynumber="$padnum2">
           <volumeref ref="volAuxDetSensitiveCRTPaddle"/>
           <position name="posCRTPaddleSensitive$paddleid2" unit="cm" 
           x="$paddle_x2" 
           y="$paddle_y2" 
           z="$paddle_z"/> 
          <rotationref ref="rIdentity"/>
       </physvol>
EOF
}
      print ENCL <<EOF;
    </volume>
EOF
}




for($imod=0 ; $imod<16 ; $imod++){
    $modnum = $imod + 1;

    print ENCL <<EOF;
    <volume name="volAuxDetCRTModuleD$modnum">
      <materialref ref="Air"/>
      <solidref ref="CRTModule"/>
EOF
for($i=0 ; $i<32 ; $i++){
    $paddle_x1 = - $CRTModWidth/2 + $CRTPaddleWidth*($i + 0.5);
    $paddle_x2 = - $CRTModWidth/2 + $CRTPaddleWidth*($i + 1);
    $paddle_y1 =   $CRTPaddleHeight/2;
    $paddle_y2 = - $CRTPaddleHeight/2;
    $paddle_z = 0;
    $padnum1 = 2*$i+64*(($modnum-1)+16);
    $padnum2 = (2*$i+1)+64*(($modnum-1)+16);
    $padnum1id = 2*$i+1;
    $padnum2id = 2*$i+2;
    $paddleid1 = "$padnum1id"."D$modnum";
    $paddleid2 = "$padnum2id"."D$modnum";
    print ENCL <<EOF;
       <physvol  name="CRTPaddle$paddleid1" copynumber="$padnum1">
           <volumeref ref="volAuxDetSensitiveCRTPaddle"/>
           <position name="posCRTPaddleSensitive$paddleid1" unit="cm" 
           x="$paddle_x1" 
           y="$paddle_y1" 
           z="$paddle_z"/> 
          <rotationref ref="rIdentity"/>
       </physvol>
       <physvol  name="CRTPaddle$paddleid2" copynumber="$padnum2">
           <volumeref ref="volAuxDetSensitiveCRTPaddle"/>
           <position name="posCRTPaddleSensitive$paddleid2" unit="cm" 
           x="$paddle_x2" 
           y="$paddle_y2" 
           z="$paddle_z"/> 
          <rotationref ref="rIdentity"/>
       </physvol>
EOF
}
      print ENCL <<EOF;
    </volume>
EOF
}
#  for($i=0 ; $i<5 ; $i++){
#  for($j=0 ; $j<5 ; $j++){
#  $suffix = "$i"."-$j";
#

if($simpleStSu == 0)
{
print ENCL <<EOF;
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
EOF

#}}

print ENCL <<EOF;
    <volume name="volSteelSupport_TB">
      <materialref ref="Air"/>
      <solidref ref="boxCryoTop"/>
EOF

  for($i=0 ; $i<5 ; $i++){
  for($j=0 ; $j<5 ; $j++){
  $xupos = -320+160*$i;
  $yupos = -320+160*$j;
  $stsuind = "$i"."-$j";
print ENCL <<EOF;
      <physvol name="volUnitTBCent_$stsuind">
        <volumeref ref="volUnitCent"/>
        <position name="posUnitTBCent_$stsuind" x="$xupos" y="$yupos" z="0" unit="cm"/>
      </physvol>
EOF
}
print ENCL <<EOF;
      <physvol name="volUnitTBE_$i">
        <volumeref ref="volUnitTop"/>
        <position name="posUnitTBE_$i" x="454.2" y="$xupos" z="0" unit="cm"/>
      </physvol>
      <physvol name="volUnitTBS_$i">
        <volumeref ref="volUnitTop"/>
        <position name="posUnitTBS_$i" x="$xupos" y="454.2" z="0" unit="cm"/>
        <rotation name="rotUnitTBS_$i" x="0" y="0" z="-90" unit="deg"/>
      </physvol>
      <physvol name="volUnitTBW_$i">
        <volumeref ref="volUnitTop"/>
        <position name="posUnitTBW_$i" x="-454.2" y="$xupos" z="0" unit="cm"/>
        <rotation name="rotUnitTBW_$i" x="0" y="0" z="-180" unit="deg"/>
      </physvol>
      <physvol name="volUnitTBN_$i">
        <volumeref ref="volUnitTop"/>
        <position name="posUnitTBN_$i" x="$xupos" y="-454.2" z="0" unit="cm"/>
        <rotation name="rotUnitTBN_$i" x="0" y="0" z="-270" unit="deg"/>
      </physvol>
EOF
   
}

      print ENCL <<EOF;
    </volume>
EOF

print ENCL <<EOF;
    <volume name="volSteelSupport_WS">
      <materialref ref="Air"/>
      <solidref ref="boxCryoWallSm"/>
EOF

  for($i=0 ; $i<5 ; $i++){
  for($j=0 ; $j<5 ; $j++){
  $xupos = -320+160*$i;
  $yupos = -320+160*$j;
  $stsuind = "$i"."-$j";
print ENCL <<EOF;
      <physvol name="volUnitWSCent_$stsuind">
        <volumeref ref="volUnitCent"/>
        <position name="posUnitWSCent_$stsuind" x="$xupos" y="$yupos" z="0" unit="cm"/>
      </physvol>
EOF
}
print ENCL <<EOF;
      <physvol name="volUnitWSE_$i">
        <volumeref ref="volUnitTop"/>
        <position name="posUnitWSE_$i" x="454.2" y="$xupos" z="0" unit="cm"/>
      </physvol>
      <physvol name="volUnitWSS_$i">
        <volumeref ref="volUnitWallS"/>
        <position name="posUnitWSS_$i" x="$xupos" y="468.9" z="0" unit="cm"/>
        <rotation name="rotUnitWSS_$i" x="0" y="0" z="-90" unit="deg"/>
      </physvol>
      <physvol name="volUnitWSW_$i">
        <volumeref ref="volUnitTop"/>
        <position name="posUnitWSW_$i" x="-454.2" y="$xupos" z="0" unit="cm"/>
        <rotation  name="rotUnitWSW_$i" x="0" y="0" z="-180" unit="deg"/>
      </physvol>
      <physvol name="volUnitWSN_$i">
        <volumeref ref="volUnitWallS"/>
        <position name="posUnitWSN_$i" x="$xupos" y="-468.9" z="0" unit="cm"/>
        <rotation  name="rotUnitWSN_$i" x="0" y="0" z="-270" unit="deg"/>
      </physvol>
EOF

}
      print ENCL <<EOF;
    </volume>
EOF

print ENCL <<EOF;
    <volume name="volSteelSupport_US">
      <materialref ref="Air"/>
      <solidref ref="boxCryoWallSmUS"/>
EOF

  for($i=0 ; $i<5 ; $i++){
  for($j=0 ; $j<5 ; $j++){
  $xupos = -320+160*$i;
  $yupos = -320+160*$j;
  $stsuind = "$i"."-$j";
print ENCL <<EOF;

      <physvol name="volUnitUSCent_$stsuind">
        <volumeref ref="volUnitCent"/>
        <position name="posUnitUSCent_$stsuind" x="$xupos" y="$yupos" z="0" unit="cm"/>
        <rotationref ref="rPlus180AboutY"/>
      </physvol>
EOF
}
print ENCL <<EOF;
      <physvol name="volUnitUSE_$i">
        <volumeref ref="volUnitTop"/>
        <position name="posUnitUSE_$i" x="454.2" y="$xupos" z="0" unit="cm"/>
        <rotationref ref="rPlus180AboutX"/>
      </physvol>
      <physvol name="volUnitUSS_$i">
        <volumeref ref="volUnitWallS"/>
        <position name="posUnitUSS_$i" x="$xupos" y="468.9" z="0" unit="cm"/>
        <rotation name="rotUnitUSS_$i" x="0" y="180" z="-90" unit="deg"/>
      </physvol>
      <physvol name="volUnitUSW_$i">
        <volumeref ref="volUnitTop"/>
        <position name="posUnitUSW_$i" x="-454.2" y="$xupos" z="0" unit="cm"/>
        <rotation  name="rotUnitUSW_$i" x="180" y="0" z="-180" unit="deg"/>
      </physvol>
      <physvol name="volUnitUSN_$i">
        <volumeref ref="volUnitWallS"/>
        <position name="posUnitUSN_$i" x="$xupos" y="-468.9" z="0" unit="cm"/>
        <rotation  name="rotUnitUSN_$i" x="0" y="180" z="-270" unit="deg"/>
      </physvol>

EOF

}
      print ENCL <<EOF;
    </volume>
EOF


print ENCL <<EOF;
    <volume name="volSteelSupport_LR">
      <materialref ref="Air"/>
      <solidref ref="boxCryoWallLg"/>
EOF

  for($i=0 ; $i<5 ; $i++){
  for($j=0 ; $j<5 ; $j++){
  $xupos = -320+160*$i;
  $yupos = -320+160*$j;
  $stsuind = "$i"."-$j";
print ENCL <<EOF;
      <physvol name="volUnitLRCent_$stsuind">
        <volumeref ref="volUnitCent"/>
        <position name="posUnitLRCent_$stsuind" x="$xupos" y="$yupos" z="0" unit="cm"/>
      </physvol>
EOF
}
print ENCL <<EOF;
      <physvol name="volUnitLRE_$i">
        <volumeref ref="volUnitWallL"/>
        <position name="posUnitLRE_$i" x="485.1" y="$xupos" z="0" unit="cm"/>
      </physvol>
      <physvol name="volUnitLRS_$i">
        <volumeref ref="volUnitWallS"/>
        <position name="posUnitLRS_$i" x="$xupos" y="468.9" z="0" unit="cm"/>
        <rotation name="rotUnitLRS_$i" x="0" y="0" z="-90" unit="deg"/>
      </physvol>
      <physvol name="volUnitLRW_$i">
        <volumeref ref="volUnitWallL"/>
        <position name="posUnitLRW_$i" x="-485.1" y="$xupos" z="0" unit="cm"/>
        <rotation  name="rotUnitLRW_$i" x="0" y="0" z="-180" unit="deg"/>
      </physvol>
      <physvol name="volUnitLRN_$i">
        <volumeref ref="volUnitWallS"/>
        <position name="posUnitLRN_$i" x="$xupos" y="-468.9" z="0" unit="cm"/>
        <rotation  name="rotUnitLRN_$i" x="0" y="0" z="-270" unit="deg"/>
      </physvol>
EOF

}
      print ENCL <<EOF;
    </volume>
EOF
} # if(simpleStSu == 0)
elsif($simpleStSu==1){
print ENCL <<EOF;
    <volume name="volSteelSupportMod">
      <materialref ref="AirSteelMixture"/>
      <solidref ref="SteelSupport"/>
    </volume>
EOF
}

    print ENCL <<EOF;

    <volume name="volDetEnclosure">
      <materialref ref="Air"/>
      <solidref ref="DetEnclosure"/>

EOF
if($simpleStSu==1){
    print ENCL <<EOF;
       <physvol>
           <volumeref ref="volSteelSupportMod"/>
           <positionref ref="posCryoInDetEnc"/>
       </physvol>
EOF
}
elsif ($simpleStSu==0)
{
    $stsuTop_y =  506.9 + $posCryoInDetEnc_y;
    $stsuBot_y = -506.9 + $posCryoInDetEnc_y;
    $stsuOther_y  =  $posCryoInDetEnc_y;
    
    print ENCL <<EOF;

       <physvol name="volSteelSupport_Top">
         <volumeref ref="volSteelSupport_TB"/>
         <position name="posSteelSupport_Top" x="0" y="$stsuTop_y" z="0" unit="cm"/>
         <rotation name="rotSteelSupport_Top" x="90" y="0" z="0" unit="deg"/>
       </physvol>
       <physvol name="volSteelSupport_Bottom">
         <volumeref ref="volSteelSupport_TB"/>
         <position name="posSteelSupport_Bottom" x="0" y="$stsuBot_y" z="0" unit="cm"/>
         <rotation name="rotSteelSupport_Bottom" x="-90" y="0" z="0" unit="deg"/>
       </physvol>
       <physvol>
         <volumeref ref="volSteelSupport_US"/>
         <position name="posSteelSupport_US" x="0" y="$stsuOther_y" z="-539.3" unit="cm"/>
       </physvol>
       <physvol name="volSteelSupport_DS">
         <volumeref ref="volSteelSupport_WS"/>
         <position name="posSteelSupport_DS" x="0" y="$stsuOther_y" z="539.3" unit="cm"/>
       </physvol>
       <physvol name="volSteelSupport_LS">
         <volumeref ref="volSteelSupport_LR"/>
         <position name="posSteelSupport_LS" x="539.3" y="$stsuOther_y" z="0" unit="cm"/>
         <rotation name="rotSteelSupport_LS" x="0" y="-90" z="0" unit="deg"/>
       </physvol>
       <physvol name="volSteelSupport_RS">
         <volumeref ref="volSteelSupport_LR"/>
         <position name="posSteelSupport_RS" x="-539.3" y="$stsuOther_y" z="0" unit="cm"/>
         <rotation name="rotSteelSupport_RS" x="0" y="90" z="0" unit="deg"/>
       </physvol>

EOF
}
    print ENCL <<EOF;

       <physvol>
           <volumeref ref="volFoamPadding"/>
           <positionref ref="posCryoInDetEnc"/>
       </physvol>
       <physvol>
         <volumeref ref="volSteelPlate"/>
         <positionref ref="posCryoInDetEnc"/>
       </physvol>
       <physvol>
         <volumeref ref="volBeamWinFoam"/>
         <position name="posBeamWinFoam" unit="cm" 
           x="$BeamWFo_x" 
           y="$BeamWFo_y" 
           z="$BeamWFo_z"/> 
         <rotationref ref="rBeamWRev3"/>
      </physvol>
      <physvol>
         <volumeref ref="volBeamWinGlassWool"/>
         <position name="posBeamWinGlassWool" unit="cm" 
           x="$BeamWGl_x" 
           y="$BeamWGl_y" 
           z="$BeamWGl_z"/> 
         <rotationref ref="rBeamWRev3"/>
      </physvol>
      <physvol>
         <volumeref ref="volBeamPipe"/>
         <position name="posBeamPipe" unit="cm" 
           x="$BeamWVa_x" 
           y="$BeamWVa_y" 
           z="$BeamWVa_z"/> 
         <rotationref ref="rBeamWRev3"/>
      </physvol>
       <physvol>
           <volumeref ref="volCryostat"/>
           <positionref ref="posCryoInDetEnc"/>
       </physvol>
       <physvol>
           <volumeref ref="volSurrConcrete"/>
           <position name="posSurrConcrete" unit="cm" 
           x="0" 
           y="$posSurrConcrete_y" 
           z="0"/>
       </physvol>
EOF

  for($i=0 ; $i<16 ; $i++){
    $modnum = $i + 1;
    $modid  = "U$modnum";

print ENCL <<EOF;
       <physvol>
           <volumeref ref="volAuxDetCRTModule$modid"/>
           <position name="posvolAuxDetCRTModule$modid" unit="cm"
            x="$posCRTUS_x[$i]" 
            y="$posCRTUS_y[$i]" 
            z="$posCRTUS_z[$i]"/>
           <rotationref ref="$posCRTUS_rot[$i]"/>
       </physvol>
EOF
  }

  for($i=0 ; $i<16 ; $i++){
    $modnum = $i + 1; 
    $modid  = "D$modnum";

print ENCL <<EOF;
       <physvol>
           <volumeref ref="volAuxDetCRTModule$modid"/>
           <position name="posvolAuxDetCRTModule$modid" unit="cm"
            x="$posCRTDS_x[$i]" 
            y="$posCRTDS_y[$i]" 
            z="$posCRTDS_z[$i]"/>
           <rotationref ref="$posCRTDS_rot[$i]"/>
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
    $WORLD = "subfile_protodune_v7_World" . $suffix . ".gdml";
    push (@gdmlFiles, $WORLD);
    $WORLD = ">" . $WORLD;
    open(WORLD) or die("Could not open file $WORLD for writing");


# The standard XML prefix and starting the gdml
    print WORLD <<EOF;
<?xml version='1.0'?>
<gdml>
EOF


# All the World solids.
$dimWorld_x = $DetEncWidth+2*$RockThickness;
$dimWorld_y = $DetEncHeight+2*$RockThickness;
$dimWorld_z = $DetEncLength+2*$RockThickness;

print WORLD <<EOF;
<solids>
    <box name="World" lunit="cm" 
      x="$dimWorld_x" 
      y="$dimWorld_y" 
      z="$dimWorld_z"/>
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



# run the sub routines that generate the fragments

gen_Extend(); #generates the GDML color extension for the refactored geometry 

gen_Define();    # generates definitions at beginning of GDML
gen_Materials(); # generates materials to be used


    gen_TPC( $TPCInner_x,  $TPC_y,  $TPC_z,  'Inner');
    gen_TPC( $TPCOuter_x,  $TPC_y,  $TPC_z,  'Outer');

close $wout;

gen_Cryostat();  # places (2*nAPAWide x nAPAHigh x nAPALong) volTPC,
                 # half rotated 180 about Y
gen_Enclosure(); # places two cryostats and concrete volumes

gen_World();     # places the enclosure among DUSEL Rock


write_fragments(); # writes the XML input for make_gdml.pl
                        # which zips together the final GDML
exit;
