#!/usr/bin/perl

#
#
#  First attempt to make a GDML fragment generator for the DUNE vertical drift 
#  10kt detector geometry with 3 views (2 orthogonal + induction at angle)
#  The lower chamber is not added yet. 
#  !!!NOTE!!!: the readout is on a positive Y plane (drift along horizontal X)
#              due to current reco limitations)
#  No photon detectors declared
#  Simplified treatment of inter-module dead spaces
#
#  Created: Thu Oct  1 16:45:27 CEST 2020
#           Vyacheslav Galymov <vgalymov@ipnl.in2p3.fr>
#
#  Modified:
#           VG: Added defs to enable use in the refactored sim framework
#           VG: 23.02.21 Adjust plane dimensions to fit a given number of ch per side
#           VG: 23.02.21 Group CRUs in CRPs
#           Laura Paulucci (lpaulucc@fnal.gov): Sept 2021 PDS added. 
#             Use option -pds=1 for backup design (membrane only coverage).
#             Default (pds=0) is the reference design (~4-pi).
#             This is linked with a larger geometry to account for photon propagation, generate it with -k=4.
#             Field Cage is turned on with reference and backup designs to match PDS option.
#	      For not including the pds, please use option -pds=-1
#
#
#################################################################################

# Each subroutine generates a fragment GDML file, and the last subroutine
# creates an XML file that make_gdml.pl will use to appropriately arrange
# the fragment GDML files to create the final desired DUNE GDML file, 
# to be named by make_gdml output command

##################################################################################


#use warnings;
use gdmlMaterials;
use Math::Trig;
use Getopt::Long;
use Math::BigFloat;
Math::BigFloat->precision(-16);

###
GetOptions( "help|h" => \$help,
	    "suffix|s:s" => \$suffix,
	    "output|o:s" => \$output,
	    "wires|w:s" => \$wires,
	    "workspace|k:s" => \$wkspc,
            "pdsconfig|pds:s" => \$pdsconfig);

my $FieldCage_switch="on";
my $Cathode_switch="off";

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


$workspace = 0;
if(defined $wkspc ) 
{
    $workspace = $wkspc;
}
elsif ( $workspace != 0 )
{
    print "\t\tCreating smaller workspace geometry.\n";
}

if ( ! defined $pdsconfig )
{
    $pdsconfig = 0;
    print "\t\tCreating reference design: 4-pi PDS converage.\n";
}
elsif ( $pdsconfig == 1 )
{
    print "\t\tCreating backup design: membrane-only PDS coverage.\n";
}

# set wires on to be the default, unless given an input by the user
$wires_on = 1; # 1=on, 0=off
if (defined $wires)
{
    $wires_on = $wires;
}

$tpc_on = 1;
$basename="_";


##################################################################
############## Parameters for One Readout Panel ##################

# parameters for 1.5 x 1.7 sub-unit Charge Readout Module / Unit
#$widthPCBActive   = 169.0; # cm 
#$lengthPCBActive  = 150.0; # cm

# views and channel counts
%nChans = ('Ind1', 256, 'Ind1Bot', 128, 'Ind2', 320, 'Col', 288);
$nViews = keys %nChans;
#print "$nViews %nChans\n";

# first induction view
$wirePitchU      = 0.8695;  # cm
$wireAngleU      = 131.63;  #-48.37;  # deg

# second induction view
$wirePitchY      = 0.525;
$widthPCBActive  = 168.00;   #$wirePitchY * $nChans{'Ind2'};

# last collection view
$wirePitchZ      = 0.517;
$lengthPCBActive = 148.9009; #$wirePitchZ * $nChans{'Col'};

#
$borderCRM       = 0.0;      # border space aroud each CRM 

$widthCRM_active  = $widthPCBActive;  
$lengthCRM_active = $lengthPCBActive; 

$widthCRM  = $widthPCBActive  + 2 * $borderCRM;
$lengthCRM = $lengthPCBActive + 2 * $borderCRM;

$borderCRP = 0.5; # cm

# number of CRMs in y and z
$nCRM_x   = 4 * 2;
$nCRM_z   = 20 * 2;

# create a smaller geometry
if( $workspace == 1 )
{
    $nCRM_x = 1 * 2;
    $nCRM_z = 1 * 2;
}

# create a smaller geometry
if( $workspace == 2 )
{
    $nCRM_x = 2 * 2;
    $nCRM_z = 2 * 2;
}

# create a smaller geometry
if( $workspace == 3 )
{
    $nCRM_x = 3 * 2;
    $nCRM_z = 3 * 2;
}

# create pds geometry
if( $workspace == 4 )
{
    $nCRM_x = 4 * 2;
    $nCRM_z = 7 * 2;
}


# calculate tpc area based on number of CRMs and their dimensions 
# each CRP should have a 2x2 CRMs
$widthTPCActive  = $nCRM_x * $widthCRM + $nCRM_x * $borderCRP;  # around 1200 for full module
$lengthTPCActive = $nCRM_z * $lengthCRM + $nCRM_z * $borderCRP; # around 6000 for full module

# active volume dimensions 
$driftTPCActive  = 650.0;

# model anode strips as wires of some diameter
$padWidth          = 0.02;
$ReadoutPlane      = $nViews * $padWidth; # 3 readout planes (no space b/w)!

##################################################################
############## Parameters for TPC and inner volume ###############

# inner volume dimensions of the cryostat
$Argon_x = 1510;
$Argon_y = 1510;
$Argon_z = 6200;

# width of gas argon layer on top
$HeightGaseousAr = 100;

if( $workspace != 0 )
{
    #active tpc + 1.0 m buffer on each side
    $Argon_x = $driftTPCActive + $HeightGaseousAr + $ReadoutPlane + 100;
    $Argon_y = $widthTPCActive + 200;
    $Argon_z = $lengthTPCActive + 200;
}


# size of liquid argon buffer
$xLArBuffer = $Argon_x - $driftTPCActive - $HeightGaseousAr - $ReadoutPlane;
$yLArBuffer = 0.5 * ($Argon_y - $widthTPCActive);
$zLArBuffer = 0.5 * ($Argon_z - $lengthTPCActive);

# cryostat 
$SteelThickness = 0.12; # membrane

$Cryostat_x = $Argon_x + 2*$SteelThickness;
$Cryostat_y = $Argon_y + 2*$SteelThickness;
$Cryostat_z = $Argon_z + 2*$SteelThickness;

##################################################################
############## DetEnc and World relevant parameters  #############

$SteelSupport_x  =  100;
$SteelSupport_y  =  100;
$SteelSupport_z  =  100; 
$FoamPadding     =  80;  
$FracMassOfSteel =  0.5; #The steel support is not a solid block, but a mixture of air and steel
$FracMassOfAir   =  1 - $FracMassOfSteel;


$SpaceSteelSupportToWall    = 100;
$SpaceSteelSupportToCeiling = 100;

$DetEncX  =    $Cryostat_x
                  + 2*($SteelSupport_x + $FoamPadding) + $SpaceSteelSupportToCeiling;

$DetEncY  =    $Cryostat_y
                  + 2*($SteelSupport_y + $FoamPadding) + 2*$SpaceSteelSupportToWall;

$DetEncZ  =    $Cryostat_z
                  + 2*($SteelSupport_z + $FoamPadding) + 2*$SpaceSteelSupportToWall;

$posCryoInDetEnc_x = - $DetEncX/2 + $SteelSupport_x + $FoamPadding + $Cryostat_x/2;


$RockThickness = 4000;

  # We want the world origin to be vertically centered on active TPC
  # This is to be added to the x and y position of every volume in volWorld

$OriginXSet =  $DetEncX/2.0
             - $SteelSupport_x
             - $FoamPadding
             - $SteelThickness
             - $xLArBuffer
             - $driftTPCActive/2.0;

$OriginYSet =   $DetEncY/2.0
              - $SpaceSteelSupportToWall
              - $SteelSupport_y
              - $FoamPadding
              - $SteelThickness
              - $yLArBuffer
              - $widthTPCActive/2.0;

  # We want the world origin to be at the very front of the fiducial volume.
  # move it to the front of the enclosure, then back it up through the concrete/foam, 
  # then through the Cryostat shell, then through the upstream dead LAr (including the
  # dead LAr on the edge of the TPC)
  # This is to be added to the z position of every volume in volWorld

$OriginZSet =   $DetEncZ/2.0 
              - $SpaceSteelSupportToWall
              - $SteelSupport_z
              - $FoamPadding
              - $SteelThickness
              - $zLArBuffer
              - $borderCRM;

##################################################################
############## Field Cage Parameters ###############
$FieldShaperLongTubeLength  =  $lengthTPCActive;
$FieldShaperShortTubeLength =  $widthTPCActive;
#$FieldShaperInnerRadius = 1.485;
#$FieldShaperOuterRadius = 1.685;
#$FieldShaperTorRad = 1.69;
$FieldShaperInnerRadius = 0.5; #cm
$FieldShaperOuterRadius = 2.285; #cm
$FieldShaperOuterRadiusSlim = 0.75; #cm
$FieldShaperTorRad = 2.3; #cm

$FieldShaperLength = $FieldShaperLongTubeLength + 2*$FieldShaperOuterRadius+ 2*$FieldShaperTorRad;
$FieldShaperWidth =  $FieldShaperShortTubeLength + 2*$FieldShaperOuterRadius+ 2*$FieldShaperTorRad;

$FieldShaperSeparation = 6.0; #cm
$NFieldShapers = ($driftTPCActive/$FieldShaperSeparation) - 1;

$FieldCageSizeX = $FieldShaperSeparation*$NFieldShapers+2;
$FieldCageSizeY = $FieldShaperWidth+2;
$FieldCageSizeZ = $FieldShaperLength+2;


####################################################################
######################## ARAPUCA Dimensions ########################
## in cm

$ArapucaOut_x = 65.0; 
$ArapucaOut_y = 2.5;
$ArapucaOut_z = 65.0; 
$ArapucaIn_x = 60.0;
$ArapucaIn_y = 2.0;
$ArapucaIn_z = 60.0;
$ArapucaAcceptanceWindow_x = 60.0;
$ArapucaAcceptanceWindow_y = 1.0;
$ArapucaAcceptanceWindow_z = 60.0;
$BlockPD_x = 0.5*$widthCRM_active - 6.35; #Sub-division of Frame (4 n x, 4 in z)
$BlockPD_z = 0.5*$lengthCRM_active - 6.25;
$GapPD = 0.5; #Size of supporting bars in Frame
$FrameToArapucaSpace       =    1.0; #At this moment, should cover the thickness of Frame + small gap to prevent overlap. VALUE NEEDS TO BE CHECKED!!!
$FrameToArapucaSpaceLat    =   $yLArBuffer - 60.0; #Arapucas 60 cm behind FC. At this moment, should cover the thickness of Frame + small gap to prevent overlap. VALUE NEEDS TO BE CHECKED!!!
$VerticalPDdist = 75.0; #distance of arapucas (center to center) in the y direction 
$HorizontalPDdist = 150.0; #distance of arapucas (center to center) in the x direction

#Positions of the 4 arapucas with respect to the Frame center --> arapucas over the cathode
$list_posx_bot[0]=-2*$BlockPD_x - 1.5*$GapPD + 0.5*$ArapucaOut_x;
$list_posz_bot[0]= 0.5*$BlockPD_z + $GapPD;
$list_posx_bot[1]=- $GapPD - 0.5*$ArapucaOut_x;
$list_posz_bot[1]=-1.5*$BlockPD_z - 1.5*$GapPD;
$list_posx_bot[2]=-$list_posx_bot[1];
$list_posz_bot[2]=-$list_posz_bot[1];
$list_posx_bot[3]=-$list_posx_bot[0];
$list_posz_bot[3]=-$list_posz_bot[0];



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


sub gen_Extend()
{

# Create the <define> fragment file name, 
# add file to list of fragments,
# and open it
    $DEF = $basename."_Ext" . $suffix . ".gdml";
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
   <rotation name="rUWireAboutX"        unit="deg" x="$wireAngleU" y="0" z="0"/>
   <rotation name="rPlus90AboutX"       unit="deg" x="90" y="0" z="0"/>
   <rotation name="rPlus90AboutY"       unit="deg" x="0" y="90" z="0"/>
   <rotation name="rPlus90AboutXPlus90AboutY" unit="deg" x="90" y="90" z="0"/>
   <rotation name="rMinus90AboutX"      unit="deg" x="270" y="0" z="0"/>
   <rotation name="rMinus90AboutY"      unit="deg" x="0" y="270" z="0"/>
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
   <D value=" 0.001205*(1-$FracMassOfSteel) + 7.9300*$FracMassOfSteel " unit="g/cm3"/>
   <fraction n="$FracMassOfSteel" ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
   <fraction n="$FracMassOfAir"   ref="Air"/>
  </material>
  <material name="vm2000" formula="vm2000">
    <D value="1.2" unit="g/cm3"/>
    <composite n="2" ref="carbon"/>
    <composite n="4" ref="hydrogen"/>
  </material>
EOF

    # add the general materials used anywere
    print MAT gdmlMaterials::gen_Materials( $asmix );

    close(MAT);
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++ gen_TPC ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# line clip on the rectangle boundary
sub lineClip {
    my $x0  = $_[0];
    my $y0  = $_[1];
    my $nx  = $_[2];
    my $ny  = $_[3];
    my $rcl = $_[4];
    my $rcw = $_[5];

    my $tol = 1.0E-4;
    my @endpts = ();
    if( abs( nx ) < tol ){
	push( @endpts, ($x0, 0) );
	push( @endpts, ($x0, $rcw) );
	return @endpts;
    }
    if( abs( ny ) < tol ){
	push( @endpts, (0, $y0) );
	push( @endpts, ($rcl, $y0) );
	return @endpts;
    }
    
    # left border at x = 0
    my $y = $y0 - $x0 * $ny/$nx;
    if( $y >= 0 && $y <= $rcw ){
	push( @endpts, (0, $y) );
    }

    # right border at x = l
    $y = $y0 + ($rcl-$x0) * $ny/$nx;
    if( $y >= 0 && $y <= $rcw ){
	push( @endpts, ($rcl, $y) );
	if( scalar(@endpts) == 4 ){
	    return @endpts;
	}
    }

    # bottom border at y = 0
    my $x = $x0 - $y0 * $nx/$ny;
    if( $x >= 0 && $x <= $rcl ){
	push( @endpts, ($x, 0) );
	if( scalar(@endpts) == 4 ){
	    return @endpts;
	}
    }
    
    # top border at y = w
    $x = $x0 + ($rcw-$y0)* $nx/$ny;
    if( $x >= 0 && $x <= $rcl ){
	push( @endpts, ($x, $rcw) );
    }
    
    return @endpts;
}

sub gen_Wires
{
    my $length = $_[0];  # 
    my $width  = $_[1];  # 
    my $nch    = $_[2];  # 
    my $nchb   = $_[3];  # nch per bottom side
    my $pitch  = $_[4];  # 
    my $theta  = $_[5];  # deg
    my $dia    = $_[6];  #
    
    $theta  = $theta * pi()/180.0;
    my @dirw   = (cos($theta), sin($theta));
    my @dirp   = (cos($theta - pi()/2), sin($theta - pi()/2));

    # calculate
    my $alpha = $theta;
    if( $alpha > pi()/2 ){
	$alpha = pi() - $alpha;
    }
    my $dX = $pitch / sin( $alpha );
    my $dY = $pitch / sin( pi()/2 - $alpha );
    if( $length <= 0 ){
        $length = $dX * $nchb;
    }
    if( $width <= 0 ){
	$width = $dY * ($nch - $nchb);
    }

    my @orig   = (0, 0);
    if( $dirp[0] < 0 ){
	$orig[0] = $length;
    }
    if( $dirp[1] < 0 ){
	$orig[1] = $width;
    }
  
    #print "origin    : @orig\n";
    #print "pitch dir : @dirp\n";
    #print "wire dir  : @dirw\n";
    #print "$length x $width cm2\n";
    
    # gen wires
    my @winfo  = ();
    my $offset = $pitch/2;
    foreach my $ch (0..$nch-1){
	#print "Processing $ch\n";

	# calculate reference point for this strip
	my @wcn = (0, 0);
	$wcn[0] = $orig[0] + $offset * $dirp[0];
	$wcn[1] = $orig[1] + $offset * $dirp[1];

	# line clip on the rectangle boundary
	@endpts = lineClip( $wcn[0], $wcn[1], $dirw[0], $dirw[1], $length, $width );

	if( scalar(@endpts) != 4 ){
	    print "Could not find end points for wire $ch : @endpts\n";
	    $offset = $offset + $pitch;
	    next;
	}

	# re-center on the mid-point
	$endpts[0] -= $length/2;
	$endpts[2] -= $length/2;
	$endpts[1] -= $width/2;
	$endpts[3] -= $width/2;

	# calculate the strip center in the rectangle of CRU
	$wcn[0] = ($endpts[0] + $endpts[2])/2;
	$wcn[1] = ($endpts[1] + $endpts[3])/2;

	# calculate the length
	my $dx = $endpts[0] - $endpts[2];
	my $dy = $endpts[1] - $endpts[3];
	my $wlen = sqrt($dx**2 + $dy**2);

	# put all info together
	my @wire = ($ch, $wcn[0], $wcn[1], $wlen);
	push( @wire, @endpts );
	push( @winfo, \@wire);
	$offset = $offset + $pitch;
	#last;
    }
    return @winfo;
}

#
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

    print " TPC dimensions     : $TPC_x x $TPC_y x $TPC_z\n";
    
    $TPC = $basename."_TPC" . $suffix . ".gdml";
    push (@gdmlFiles, $TPC);
    $TPC = ">" . $TPC;
    open(TPC) or die("Could not open file $TPC for writing");

    # The standard XML prefix and starting the gdml
print TPC <<EOF;
    <?xml version='1.0'?>
	<gdml>
EOF

    # compute wires for 1st induction
    my @winfo = ();
    if( $wires_on == 1 ){
	@winfo = gen_Wires( 0, 0, # $TPCActive_y,
			    $nChans{'Ind1'}, $nChans{'Ind1Bot'}, 
			    $wirePitchU, $wireAngleU, $padWidth );
    }

    # All the TPC solids save the wires.
print TPC <<EOF;
    <solids>
EOF

print TPC <<EOF;
   <box name="CRM"
      x="$TPC_x" 
      y="$TPC_y" 
      z="$TPC_z"
      lunit="cm"/>
   <box name="CRMUPlane" 
      x="$padWidth" 
      y="$TPCActive_y" 
      z="$TPCActive_z"
      lunit="cm"/>
   <box name="CRMYPlane" 
      x="$padWidth" 
      y="$TPCActive_y" 
      z="$TPCActive_z"
      lunit="cm"/>
   <box name="CRMZPlane" 
      x="$padWidth"
      y="$TPCActive_y"
      z="$TPCActive_z"
      lunit="cm"/>
   <box name="CRMActive" 
      x="$TPCActive_x"
      y="$TPCActive_y"
      z="$TPCActive_z"
      lunit="cm"/>
EOF

#++++++++++++++++++++++++++++ Wire Solids ++++++++++++++++++++++++++++++
if($wires_on==1){
	    
    foreach my $wire (@winfo) {
	my $wid = $wire->[0];
	my $wln = $wire->[3];
print TPC <<EOF;
   <tube name="CRMWireU$wid"
      rmax="0.5*$padWidth"
      z="$wln"               
      deltaphi="360"
      aunit="deg" lunit="cm"/>
EOF
    }
    
print TPC <<EOF;
   <tube name="CRMWireY"
      rmax="0.5*$padWidth"
      z="$TPCActive_z"               
      deltaphi="360" 
      aunit="deg" lunit="cm"/>
   <tube name="CRMWireZ"
      rmax="0.5*$padWidth"
      z="$TPCActive_y"               
      deltaphi="360"
      aunit="deg" lunit="cm"/>
EOF
}
print TPC <<EOF;
</solids>
EOF


# Begin structure and create wire logical volumes
print TPC <<EOF;
<structure>
    <volume name="volTPCActive">
      <materialref ref="LAr"/>
      <solidref ref="CRMActive"/>
      <auxiliary auxtype="SensDet" auxvalue="SimEnergyDeposit"/>
      <auxiliary auxtype="StepLimit" auxunit="cm" auxvalue="0.5208*cm"/>
      <auxiliary auxtype="Efield" auxunit="V/cm" auxvalue="500*V/cm"/>
      <colorref ref="blue"/>
    </volume>
EOF

if($wires_on==1) 
{
    foreach my $wire (@winfo) 
    {
	my $wid = $wire->[0];
print TPC <<EOF;
    <volume name="volTPCWireU$wid">
      <materialref ref="Copper_Beryllium_alloy25"/>
      <solidref ref="CRMWireU$wid"/>
    </volume>
EOF
    }

print TPC <<EOF;
    <volume name="volTPCWireY">
      <materialref ref="Copper_Beryllium_alloy25"/>
      <solidref ref="CRMWireY"/>
    </volume>
    <volume name="volTPCWireZ">
      <materialref ref="Copper_Beryllium_alloy25"/>
      <solidref ref="CRMWireZ"/>
    </volume>
EOF
}
    # 1st induction plane
print TPC <<EOF;
   <volume name="volTPCPlaneU">
     <materialref ref="LAr"/>
     <solidref ref="CRMUPlane"/>
EOF
if ($wires_on==1) # add wires to U plane 
{
    # the coordinates were computed with a corner at (0,0)
    # so we need to move to plane coordinates
    my $offsetZ = 0; #-0.5 * $TPCActive_z;
    my $offsetY = 0; #-0.5 * $TPCActive_y;

    foreach my $wire (@winfo) {
	my $wid  = $wire->[0];
	my $zpos = $wire->[1] + $offsetZ;
	my $ypos = $wire->[2] + $offsetY;
print TPC <<EOF;
     <physvol>
       <volumeref ref="volTPCWireU$wid"/> 
       <position name="posWireU$wid" unit="cm" x="0" y="$ypos" z="$zpos"/>
       <rotationref ref="rUWireAboutX"/> 
     </physvol>
EOF
    }
}
print TPC <<EOF;
   </volume>
EOF

# 2nd induction plane
print TPC <<EOF;
  <volume name="volTPCPlaneY">
    <materialref ref="LAr"/>
    <solidref ref="CRMYPlane"/>
EOF

if ($wires_on==1) # add wires to Y plane (plane with wires reading y position)
{
    for(my $i=0;$i<$nChans{'Ind2'};++$i)
    {
	#my $ypos = -0.5 * $TPCActive_y + ($i+0.5)*$wirePitchY + 0.5*$padWidth;
	my $ypos = ($i + 0.5 - $nChans{'Ind2'}/2)*$wirePitchY;
	if( (0.5 * $TPCActive_y - abs($ypos)) < 0 ){
	    die "Cannot place wire $i in view Y, as plane is too small\n";
	}
print TPC <<EOF;
      <physvol>
        <volumeref ref="volTPCWireY"/> 
        <position name="posWireY$i" unit="cm" x="0" y="$ypos" z="0"/>
	<rotationref ref="rIdentity"/> 
      </physvol>
EOF
   }
}
print TPC <<EOF;
  </volume>
EOF

# collection plane
print TPC <<EOF;
  <volume name="volTPCPlaneZ">
    <materialref ref="LAr"/>
    <solidref ref="CRMZPlane"/>
EOF
if ($wires_on==1) # add wires to Z plane (plane with wires reading z position)
   {
       for(my $i=0;$i<$nChans{'Col'};++$i)
       {
	   #my $zpos = -0.5 * $TPCActive_z + ($i+0.5)*$wirePitchZ + 0.5*$padWidth;
	   my $zpos = ($i + 0.5 - $nChans{'Col'}/2)*$wirePitchZ;
	   if( (0.5 * $TPCActive_z - abs($zpos)) < 0 ){
	       die "Cannot place wire $i in view Z, as plane is too small\n";
	   }
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

       
$posUplane[0] = 0.5*$TPC_x - 2.5*$padWidth;
$posUplane[1] = 0;
$posUplane[2] = 0;

$posYplane[0] = 0.5*$TPC_x - 1.5*$padWidth;
$posYplane[1] = 0;
$posYplane[2] = 0;

$posZplane[0] = 0.5*$TPC_x - 0.5*$padWidth;
$posZplane[1] = 0; 
$posZplane[2] = 0;

$posTPCActive[0] = -$ReadoutPlane/2;
$posTPCActive[1] = 0;
$posTPCActive[2] = 0;


#wrap up the TPC file
print TPC <<EOF;
   <volume name="volTPC">
     <materialref ref="LAr"/>
       <solidref ref="CRM"/>
       <physvol>
       <volumeref ref="volTPCPlaneU"/>
       <position name="posPlaneU" unit="cm" 
         x="$posUplane[0]" y="$posUplane[1]" z="$posUplane[2]"/>
       <rotationref ref="rIdentity"/>
     </physvol>
     <physvol>
       <volumeref ref="volTPCPlaneY"/>
       <position name="posPlaneY" unit="cm" 
         x="$posYplane[0]" y="$posYplane[1]" z="$posYplane[2]"/>
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
        x="$posTPCActive[0]" y="$posTPCAtive[1]" z="$posTPCActive[2]"/>
       <rotationref ref="rIdentity"/>
     </physvol>
   </volume>
EOF
## x="@{[$posTPCActive[0]+$padWidth]}" y="$posTPCActive[1]" z="$posTPCActive[2]"/>

print TPC <<EOF;
 </structure>
 </gdml>
EOF

    close(TPC);
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++ gen_FieldCage ++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub gen_FieldCage {

    $FieldCage = $basename."_FieldCage" . $suffix . ".gdml";
    push (@gdmlFiles, $FieldCage);
    $FieldCage = ">" . $FieldCage;
    open(FieldCage) or die("Could not open file $FieldCage for writing");

# The standard XML prefix and starting the gdml
print FieldCage <<EOF;
   <?xml version='1.0'?>
   <gdml>
EOF
# The printing solids used in the Field Cage
#print "lengthTPCActive      : $lengthTPCActive \n";
#print "widthTPCActive       : $widthTPCActive \n";


print FieldCage <<EOF;
<solids>
     <torus name="FieldShaperCorner" rmin="$FieldShaperInnerRadius" rmax="$FieldShaperOuterRadius" rtor="$FieldShaperTorRad" deltaphi="90" startphi="0" aunit="deg" lunit="cm"/>
     <tube name="FieldShaperLongtube" rmin="$FieldShaperInnerRadius" rmax="$FieldShaperOuterRadius" z="$FieldShaperLongTubeLength" deltaphi="360" startphi="0" aunit="deg" lunit="cm"/>
     <tube name="FieldShaperLongtubeSlim" rmin="$FieldShaperInnerRadius" rmax="$FieldShaperOuterRadiusSlim" z="$FieldShaperLongTubeLength" deltaphi="360" startphi="0" aunit="deg" lunit="cm"/>
     <tube name="FieldShaperShorttube" rmin="$FieldShaperInnerRadius" rmax="$FieldShaperOuterRadius" z="$FieldShaperShortTubeLength" deltaphi="360" startphi="0" aunit="deg" lunit="cm"/>

    <union name="FSunion1">
      <first ref="FieldShaperLongtube"/>
      <second ref="FieldShaperCorner"/>
   		<position name="esquinapos1" unit="cm" x="@{[-$FieldShaperTorRad]}" y="0" z="@{[0.5*$FieldShaperLongTubeLength]}"/>
		<rotation name="rot1" unit="deg" x="90" y="0" z="0" />
    </union>

    <union name="FSunion2">
      <first ref="FSunion1"/>
      <second ref="FieldShaperShorttube"/>
   		<position name="esquinapos2" unit="cm" x="@{[-0.5*$FieldShaperShortTubeLength-$FieldShaperTorRad]}" y="0" z="@{[+0.5*$FieldShaperLongTubeLength+$FieldShaperTorRad]}"/>
   		<rotation name="rot2" unit="deg" x="0" y="90" z="0" />
    </union>

    <union name="FSunion3">
      <first ref="FSunion2"/>
      <second ref="FieldShaperCorner"/>
   		<position name="esquinapos3" unit="cm" x="@{[-$FieldShaperShortTubeLength-$FieldShaperTorRad]}" y="0" z="@{[0.5*$FieldShaperLongTubeLength]}"/>
		<rotation name="rot3" unit="deg" x="90" y="270" z="0" />
    </union>

    <union name="FSunion4">
      <first ref="FSunion3"/>
      <second ref="FieldShaperLongtube"/>
   		<position name="esquinapos4" unit="cm" x="@{[-$FieldShaperShortTubeLength-2*$FieldShaperTorRad]}" y="0" z="0"/>
    </union>

    <union name="FSunion5">
      <first ref="FSunion4"/>
      <second ref="FieldShaperCorner"/>
   		<position name="esquinapos5" unit="cm" x="@{[-$FieldShaperShortTubeLength-$FieldShaperTorRad]}" y="0" z="@{[-0.5*$FieldShaperLongTubeLength]}"/>
		<rotation name="rot5" unit="deg" x="90" y="180" z="0" />
    </union>

    <union name="FSunion6">
      <first ref="FSunion5"/>
      <second ref="FieldShaperShorttube"/>
   		<position name="esquinapos6" unit="cm" x="@{[-0.5*$FieldShaperShortTubeLength-$FieldShaperTorRad]}" y="0" z="@{[-0.5*$FieldShaperLongTubeLength-$FieldShaperTorRad]}"/>
		<rotation name="rot6" unit="deg" x="0" y="90" z="0" />
    </union>

    <union name="FieldShaperSolid">
      <first ref="FSunion6"/>
      <second ref="FieldShaperCorner"/>
   		<position name="esquinapos7" unit="cm" x="@{[-$FieldShaperTorRad]}" y="0" z="@{[-0.5*$FieldShaperLongTubeLength]}"/>
		<rotation name="rot7" unit="deg" x="90" y="90" z="0" />
    </union>
    
    <union name="FSunionSlim1">
      <first ref="FieldShaperLongtubeSlim"/>
      <second ref="FieldShaperCorner"/>
   		<position name="esquinapos1" unit="cm" x="@{[-$FieldShaperTorRad]}" y="0" z="@{[0.5*$FieldShaperLongTubeLength]}"/>
		<rotation name="rot1" unit="deg" x="90" y="0" z="0" />
    </union>

    <union name="FSunionSlim2">
      <first ref="FSunionSlim1"/>
      <second ref="FieldShaperShorttube"/>
   		<position name="esquinapos2" unit="cm" x="@{[-0.5*$FieldShaperShortTubeLength-$FieldShaperTorRad]}" y="0" z="@{[+0.5*$FieldShaperLongTubeLength+$FieldShaperTorRad]}"/>
   		<rotation name="rot2" unit="deg" x="0" y="90" z="0" />
    </union>

    <union name="FSunionSlim3">
      <first ref="FSunionSlim2"/>
      <second ref="FieldShaperCorner"/>
   		<position name="esquinapos3" unit="cm" x="@{[-$FieldShaperShortTubeLength-$FieldShaperTorRad]}" y="0" z="@{[0.5*$FieldShaperLongTubeLength]}"/>
		<rotation name="rot3" unit="deg" x="90" y="270" z="0" />
    </union>

    <union name="FSunionSlim4">
      <first ref="FSunionSlim3"/>
      <second ref="FieldShaperLongtubeSlim"/>
   		<position name="esquinapos4" unit="cm" x="@{[-$FieldShaperShortTubeLength-2*$FieldShaperTorRad]}" y="0" z="0"/>
    </union>

    <union name="FSunionSlim5">
      <first ref="FSunionSlim4"/>
      <second ref="FieldShaperCorner"/>
   		<position name="esquinapos5" unit="cm" x="@{[-$FieldShaperShortTubeLength-$FieldShaperTorRad]}" y="0" z="@{[-0.5*$FieldShaperLongTubeLength]}"/>
		<rotation name="rot5" unit="deg" x="90" y="180" z="0" />
    </union>

    <union name="FSunionSlim6">
      <first ref="FSunionSlim5"/>
      <second ref="FieldShaperShorttube"/>
   		<position name="esquinapos6" unit="cm" x="@{[-0.5*$FieldShaperShortTubeLength-$FieldShaperTorRad]}" y="0" z="@{[-0.5*$FieldShaperLongTubeLength-$FieldShaperTorRad]}"/>
		<rotation name="rot6" unit="deg" x="0" y="90" z="0" />
    </union>

    <union name="FieldShaperSolidSlim">
      <first ref="FSunionSlim6"/>
      <second ref="FieldShaperCorner"/>
   		<position name="esquinapos7" unit="cm" x="@{[-$FieldShaperTorRad]}" y="0" z="@{[-0.5*$FieldShaperLongTubeLength]}"/>
		<rotation name="rot7" unit="deg" x="90" y="90" z="0" />
    </union>
    
</solids>

EOF

print FieldCage <<EOF;

<structure>
<volume name="volFieldShaper">
  <materialref ref="Al2O3"/>
  <solidref ref="FieldShaperSolid"/>
</volume>
<volume name="volFieldShaperSlim">
  <materialref ref="Al2O3"/>
  <solidref ref="FieldShaperSolidSlim"/>
</volume>

</structure>

EOF

print FieldCage <<EOF;

</gdml>
EOF
close(FieldCage);
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
# External active are two side volumes for generating light outside the field cage (no top or bottom buffers included)
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

    <box name="ExternalAuxOut" lunit="cm" 
      x="$Argon_x - $xLArBuffer"
      y="$Argon_y - 2*$ArapucaOut_y - 2*$FrameToArapucaSpaceLat"
      z="$Argon_z"/>

    <box name="ExternalAuxIn" lunit="cm" 
      x="$Argon_x"
      y="$FieldCageSizeY"
      z="$Argon_z + 1"/>

   <subtraction name="ExternalActive">
      <first ref="ExternalAuxOut"/>
      <second ref="ExternalAuxIn"/>
    </subtraction>

    <subtraction name="SteelShell">
      <first ref="Cryostat"/>
      <second ref="ArgonInterior"/>
    </subtraction>

</solids>
EOF
    
#PDS
#Double sided detectors should only be included when both top and bottom volumes become available
#Optical sensitive volumes cannot be rotated because Larsoft cannot pick up the rotation when obtinaing the lengths needed for the semi-analytic model --> two acceptance windows for single sided lateral and cathode
print CRYO <<EOF;
<solids>
    <box name="ArapucaOut" lunit="cm"
      x="@{[$ArapucaOut_x]}"
      y="@{[$ArapucaOut_y]}"
      z="@{[$ArapucaOut_z]}"/>

    <box name="ArapucaIn" lunit="cm"
      x="@{[$ArapucaIn_x]}"
      y="@{[$ArapucaOut_y]}"
      z="@{[$ArapucaIn_z]}"/>

     <subtraction name="ArapucaWalls">
      <first  ref="ArapucaOut"/>
      <second ref="ArapucaIn"/>
      <position name="posArapucaSub" x="0" y="@{[$ArapucaOut_y/2.0]}" z="0." unit="cm"/>
      </subtraction>

    <box name="ArapucaAcceptanceWindow" lunit="cm"
      x="@{[$ArapucaAcceptanceWindow_x]}"
      y="@{[$ArapucaAcceptanceWindow_y]}"
      z="@{[$ArapucaAcceptanceWindow_z]}"/>

    <box name="ArapucaDoubleIn" lunit="cm"
      x="@{[$ArapucaIn_x]}"
      y="@{[$ArapucaOut_y+1.0]}"
      z="@{[$ArapucaIn_z]}"/>

     <subtraction name="ArapucaDoubleWalls">
      <first  ref="ArapucaOut"/>
      <second ref="ArapucaDoubleIn"/>
      <position name="posArapucaDoubleSub" x="0" y="0" z="0" unit="cm"/>
      </subtraction>

    <box name="ArapucaDoubleAcceptanceWindow" lunit="cm"
      x="@{[$ArapucaOut_y-0.02]}"
      y="@{[$ArapucaAcceptanceWindow_x]}"
      z="@{[$ArapucaAcceptanceWindow_z]}"/>

    <box name="ArapucaCathodeAcceptanceWindow" lunit="cm"
      x="@{[$ArapucaAcceptanceWindow_y]}"
      y="@{[$ArapucaAcceptanceWindow_x]}"
      z="@{[$ArapucaAcceptanceWindow_z]}"/>

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
    <volume name="volExternalActive">
      <materialref ref="LAr"/>
      <solidref ref="ExternalActive"/>
      <auxiliary auxtype="SensDet" auxvalue="SimEnergyDeposit"/>
      <auxiliary auxtype="StepLimit" auxunit="cm" auxvalue="0.5208*cm"/>
      <auxiliary auxtype="Efield" auxunit="V/cm" auxvalue="0*V/cm"/>
      <colorref ref="green"/>
    </volume>
EOF
#including single sided arapucas over the cathode while there is only the top volume
#if double sided, use
#    <volume name="volArapucaDouble_$i\-$j\-$p">
#      <materialref ref="G10" />
#      <solidref ref="ArapucaDoubleWalls" />
#    </volume>
#    <volume name="volOpDetSensitive_ArapucaDouble_$i\-$j\-$p">
#      <materialref ref="LAr"/>
#      <solidref ref="ArapucaDoubleAcceptanceWindow"/>
#    </volume>
if ($pdsconfig == 0){  #4-pi PDS converage
for($i=0 ; $i<$nCRM_x/2 ; $i++){ #arapucas over the cathode
for($j=0 ; $j<$nCRM_z/2 ; $j++){
for($p=0 ; $p<4 ; $p++){
  print CRYO <<EOF;
    <volume name="volArapucaDouble_$i\-$j\-$p">
      <materialref ref="G10" />
      <solidref ref="ArapucaWalls" />
    </volume>
    <volume name="volOpDetSensitive_ArapucaDouble_$i\-$j\-$p">
      <materialref ref="LAr"/>
      <solidref ref="ArapucaCathodeAcceptanceWindow"/>
    </volume>
EOF
}
}
}
}

if ($nCRM_x==8){ #arapucas on the laterals
for($j=0 ; $j<$nCRM_z/2 ; $j++){
for($p=0 ; $p<8 ; $p++){
  print CRYO <<EOF;
    <volume name="volArapucaLat_$j\-$p">
      <materialref ref="G10" />
      <solidref ref="ArapucaWalls" />
    </volume>
    <volume name="volOpDetSensitive_ArapucaLat_$j\-$p">
      <materialref ref="LAr"/>
      <solidref ref="ArapucaAcceptanceWindow"/>
    </volume>
EOF
}
}
}

if ($pdsconfig == 1){  #Membrane PDS converage
if ($nCRM_x==8) { #arapucas on the laterals
for($j=0 ; $j<$nCRM_z/2 ; $j++){
for($p=8 ; $p<18 ; $p++){
  print CRYO <<EOF;
    <volume name="volArapucaLat_$j\-$p">
      <materialref ref="G10" />
      <solidref ref="ArapucaWalls" />
    </volume>
    <volume name="volOpDetSensitive_ArapucaLat_$j\-$p">
      <materialref ref="LAr"/>
      <solidref ref="ArapucaAcceptanceWindow"/>
    </volume>
EOF
}
}
}
}

      print CRYO <<EOF;

    <volume name="volCryostat">
      <materialref ref="LAr" />
      <solidref ref="Cryostat" />
      <physvol>
        <volumeref ref="volGaseousArgon"/>
        <position name="posGaseousArgon" unit="cm" x="@{[$Argon_x/2-$HeightGaseousAr/2]}" y="0" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volSteelShell"/>
        <position name="posSteelShell" unit="cm" x="0" y="0" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volExternalActive"/>
        <position name="posExternalActive" unit="cm" x="-$xLArBuffer/2" y="0" z="0"/>
      </physvol>
EOF

if ($tpc_on==1) # place TPC inside croysotat offsetting each pair of CRMs by borderCRP
{
  $posX =  $Argon_x/2 - $HeightGaseousAr - 0.5*($driftTPCActive + $ReadoutPlane); 
  $idx = 0;
  my $posZ = -0.5*$Argon_z + $zLArBuffer + 0.5*$lengthCRM;
  for(my $ii=0;$ii<$nCRM_z;$ii++)
  {
    if( $ii % 2 == 0 ){
	$posZ += $borderCRP;
	if( $ii>0 ){
	    $posZ += $borderCRP;
	}
    }
    my $posY = -0.5*$Argon_y + $yLArBuffer + 0.5*$widthCRM;
    for(my $jj=0;$jj<$nCRM_x;$jj++)
    {
	if( $jj % 2 == 0 ){
	    $posY += $borderCRP;
	    if( $jj>0 ){
		$posY += $borderCRP;
	    }
	}
	print CRYO <<EOF;
      <physvol>
        <volumeref ref="volTPC"/>
	<position name="posTPC\-$idx" unit="cm"
           x="$posX" y="$posY" z="$posZ"/>
      </physvol>
EOF
       $idx++;
       $posY += $widthCRM;
    }

    $posZ += $lengthCRM;
  }
}

#The +50 in the x positions must depend on some other parameter
  if ( $FieldCage_switch eq "on" ) {
    for ( $i=0; $i<$NFieldShapers; $i=$i+1 ) {
    $dist=$i*$FieldShaperSeparation;
$posX =  $Argon_x/2 - $HeightGaseousAr - 0.5*($driftTPCActive + $ReadoutPlane); 
	if ($pdsconfig==0){
		if ($dist>250){
	print CRYO <<EOF;
  <physvol>
     <volumeref ref="volFieldShaperSlim"/>
     <position name="posFieldShaper$i" unit="cm"  x="@{[-$OriginXSet+50+($i-$NFieldShapers*0.5)*$FieldShaperSeparation]}" y="@{[-0.5*$FieldShaperShortTubeLength-$FieldShaperTorRad]}" z="0" />
     <rotation name="rotFS$i" unit="deg" x="0" y="0" z="90" />
  </physvol>
EOF
		}else{
	print CRYO <<EOF;
  <physvol>
     <volumeref ref="volFieldShaper"/>
     <position name="posFieldShaper$i" unit="cm"  x="@{[-$OriginXSet+50+($i-$NFieldShapers*0.5)*$FieldShaperSeparation]}" y="@{[-0.5*$FieldShaperShortTubeLength-$FieldShaperTorRad]}" z="0" />
     <rotation name="rotFS$i" unit="deg" x="0" y="0" z="90" />
  </physvol>
EOF
		}
	}else{
	print CRYO <<EOF;
  <physvol>
     <volumeref ref="volFieldShaperSlim"/>
     <position name="posFieldShaper$i" unit="cm"  x="@{[-$OriginXSet+50+($i-$NFieldShapers*0.5)*$FieldShaperSeparation]}" y="@{[-0.5*$FieldShaperShortTubeLength-$FieldShaperTorRad]}" z="0" />
     <rotation name="rotFS$i" unit="deg" x="0" y="0" z="90" />
  </physvol>
EOF
	}
    }
  }


$CathodePosX =-$OriginXSet+50+(-1-$NFieldShapers*0.5)*$FieldShaperSeparation;
$CathodePosY = 0;
$CathodePosZ = 0;
  if ( $Cathode_switch eq "on" )
  {
      print CRYO <<EOF;
  <physvol>
   <volumeref ref="volGroundGrid"/>
   <position name="posGroundGrid01" unit="cm" x="$CathodePosX" y="@{[-$CathodePosY]}" z="@{[$CathodePosZ]}"/>
   <rotation name="rotGG01" unit="deg" x="0" y="0" z="90" />
  </physvol>

EOF
  }

$FrameLenght_x= 2.0*$widthCRM_active;
$FrameLenght_z= 2.0*$lengthCRM_active;
  
if ($pdsconfig == 0) {  #4-pi PDS converage

#for placing the Arapucas over the cathode
  $FrameCenter_x=-1.5*$FrameLenght_x+(4-$nCRM_x/2)/2*$FrameLenght_x;
  $FrameCenter_y=$CathodePosX;
  $FrameCenter_z=-9.5*$FrameLenght_z+(20-$nCRM_z/2)/2*$FrameLenght_z;
for($i=0;$i<$nCRM_x/2;$i++){
for($j=0;$j<$nCRM_z/2;$j++){
  place_OpDetsCathode($FrameCenter_x, $FrameCenter_y, $FrameCenter_z, $i, $j);
  $FrameCenter_z+=$FrameLenght_z;
}
  $FrameCenter_x+=$FrameLenght_x;
  $FrameCenter_z=-9.5*$FrameLenght_z+(20-$nCRM_z/2)/2*$FrameLenght_z;
}
}

if ($pdsconfig == 0) {  #4-pi PDS converage
#for placing the Arapucas on laterals
if ($nCRM_x==8) {
  $FrameCenter_y=$posZplane[0]; #anode position
  $FrameCenter_z=-19*$FrameLenght_z/2+(40-$nCRM_z)/2*$FrameLenght_z/2;
for($j=0;$j<$nCRM_z/2;$j++){#nCRM will give the collumn number (1 collumn per frame)
  place_OpDetsLateral($FrameCenter_y, $FrameCenter_z, $j);
  $FrameCenter_z+=$FrameLenght_z;
}
}

} else {  #membrane only PDS converage

if($pdsconfig == 1){
if ($nCRM_x==8) {
  $FrameCenter_y=$posZplane[0]; #anode position
  $FrameCenter_z=-19*$FrameLenght_z/2+(40-$nCRM_z)/2*$FrameLenght_z/2;
for($j=0;$j<$nCRM_z/2;$j++){#nCRM will give the collumn number (1 collumn per frame)
  place_OpDetsMembOnly($FrameCenter_y, $FrameCenter_z, $j);
  $FrameCenter_z+=$FrameLenght_z;
}
}
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
#++++++++++++++++++++++++++++++++++++ place_OpDets +++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sub place_OpDetsCathode()
{

    $FrameCenter_x = $_[0];
    $FrameCenter_y = $_[1];
    $FrameCenter_z = $_[2];
    $Frame_x = $_[3];
    $Frame_z = $_[4];

#Placing Arapucas over the Cathode 
#If there are both top and bottom volumes --> use double-sided:
#     <physvol>
#       <volumeref ref="volOpDetSensitive_ArapucaDouble_$Frame_x\-$Frame_z\-$ara"/>
#       <position name="posOpArapucaDouble$ara-Frame\-$Frame_x\-$Frame_z" unit="cm" 
#         x="@{[$Ara_X]}"
#	 y="@{[$Ara_Y]}" 
#	 z="@{[$Ara_Z]}"/>
#     </physvol>
#else
for ($ara = 0; $ara<4; $ara++)
{
             # All Arapuca centers will have the same Y coordinate
             # X and Z coordinates are defined with respect to the center of the current Frame

 	     $Ara_Y = $FrameCenter_x+$list_posx_bot[$ara]; #GEOMETRY IS ROTATED: X--> Y AND Y--> X
             $Ara_X = $FrameCenter_y + $FrameToArapucaSpace;
 	     $Ara_Z = $FrameCenter_z+$list_posz_bot[$ara];

	print CRYO <<EOF;
     <physvol>
       <volumeref ref="volArapucaDouble_$Frame_x\-$Frame_z\-$ara"/>
       <position name="posArapucaDouble$ara-Frame\-$Frame_x\-$Frame_z" unit="cm" 
         x="@{[$Ara_X]}"
	 y="@{[$Ara_Y]}" 
	 z="@{[$Ara_Z]}"/>
       <rotation name="rPlus90AboutXPlus90AboutZ" unit="deg" x="90" y="0" z="90"/>
     </physvol>
     <physvol>
       <volumeref ref="volOpDetSensitive_ArapucaDouble_$Frame_x\-$Frame_z\-$ara"/>
       <position name="posOpArapucaDouble$ara-Frame\-$Frame_x\-$Frame_z" unit="cm" 
         x="@{[$Ara_X+0.5*$ArapucaOut_y-0.5*$ArapucaAcceptanceWindow_y-0.01]}"
	 y="@{[$Ara_Y]}" 
	 z="@{[$Ara_Z]}"/>
     </physvol>
EOF

}#end Ara for-loop

}


sub place_OpDetsLateral()
{

    $FrameCenter_y = $_[0];
    $FrameCenter_z = $_[1];
    $Lat_z = $_[2];

#Placing Arapucas on the laterals if nCRM_x=8 -- Single Sided
for ($ara = 0; $ara<8; $ara++)
{
             # Arapucas on laterals
             # All Arapuca centers on a given collumn will have the same Z coordinate
             # X coordinates are on the left and right laterals
             # Y coordinates are defined with respect to the cathode position
             # There are two collumns per frame on each side.

             if ($ara<4) {$Ara_Y = -0.5*$Argon_y + $FrameToArapucaSpaceLat;
                         $Ara_YSens = ($Ara_Y+0.5*$ArapucaOut_y-0.5*$ArapucaAcceptanceWindow_y-0.01);
                         $rot= "rIdentity"; }
             else {$Ara_Y = 0.5*$Argon_y - $FrameToArapucaSpaceLat;
                         $Ara_YSens = ($Ara_Y-0.5*$ArapucaOut_y+0.5*$ArapucaAcceptanceWindow_y+0.01);
                         $rot = "rPlus180AboutX";} #GEOMETRY IS ROTATED: X--> Y AND Y--> X
             if ($ara==0||$ara==4) {$Ara_X = $FrameCenter_y-40.0;} #first tile's center 40 cm bellow anode
             else {$Ara_X-=$VerticalPDdist;} #other tiles separated by VerticalPDdist
             $Ara_Z = $FrameCenter_z;

#        print "lateral arapucas: $Ara_X, $Ara_Y, $Ara_Z \n";
        
	print CRYO <<EOF;
     <physvol>
       <volumeref ref="volArapucaLat_$Lat_z\-$ara"/>
       <position name="posArapuca$ara-Lat\-$Lat_z" unit="cm" 
         x="@{[$Ara_X]}"
	 y="@{[$Ara_Y]}" 
	 z="@{[$Ara_Z]}"/>
       <rotationref ref="$rot"/>
     </physvol>
     <physvol>
       <volumeref ref="volOpDetSensitive_ArapucaLat_$Lat_z\-$ara"/>
       <position name="posOpArapuca$ara-Lat\-$Lat_z" unit="cm" 
         x="@{[$Ara_X]}"
	 y="@{[$Ara_YSens]}" 
	 z="@{[$Ara_Z]}"/>
     </physvol>
EOF
        
}#end Ara for-loop

}


sub place_OpDetsMembOnly()
{

    $FrameCenter_y = $_[0];
    $FrameCenter_z = $_[1];
    $Lat_z = $_[2];

#Placing Arapucas on the laterals if nCRM_x=8 -- Single Sided
for ($ara = 0; $ara<18; $ara++)
{
             # Arapucas on laterals
             # All Arapuca centers on a given collumn will have the same Z coordinate
             # X coordinates are on the left and right laterals
             # Y coordinates are defined with respect to the cathode position
             # There are two collumns per frame on each side.

             if($ara<9) {$Ara_Y = -0.5*$Argon_y + $FrameToArapucaSpaceLat;
                         $Ara_YSens = ($Ara_Y+0.5*$ArapucaOut_y-0.5*$ArapucaAcceptanceWindow_y-0.01);
                         $rot= "rIdentity"; }
             else {$Ara_Y = 0.5*$Argon_y - $FrameToArapucaSpaceLat;
                         $Ara_YSens = ($Ara_Y-0.5*$ArapucaOut_y+0.5*$ArapucaAcceptanceWindow_y+0.01);
                         $rot = "rPlus180AboutX";} #GEOMETRY IS ROTATED: X--> Y AND Y--> X
             if($ara==0||$ara==9) {$Ara_X = $FrameCenter_y-$ArapucaOut_x/2;} #first tile's center right below anode
             else {$Ara_X-=$ArapucaOut_x - $FrameToArapucaSpace;} #other tiles separated by minimal distance + buffer
             $Ara_Z = $FrameCenter_z;

#        print "lateral arapucas: $Ara_X, $Ara_Y, $Ara_Z \n";
        
	print CRYO <<EOF;
     <physvol>
       <volumeref ref="volArapucaLat_$Lat_z\-$ara"/>
       <position name="posArapuca$ara-Lat\-$Lat_z" unit="cm" 
         x="@{[$Ara_X]}"
	 y="@{[$Ara_Y]}" 
	 z="@{[$Ara_Z]}"/>
       <rotationref ref="$rot"/>
     </physvol>
     <physvol>
       <volumeref ref="volOpDetSensitive_ArapucaLat_$Lat_z\-$ara"/>
       <position name="posOpArapuca$ara-Lat\-$Lat_z" unit="cm" 
         x="@{[$Ara_X]}"
	 y="@{[$Ara_YSens]}" 
	 z="@{[$Ara_Z]}"/>
     </physvol>
EOF
        
}#end Ara for-loop



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
      x="$DetEncX"
      y="$DetEncY"
      z="$DetEncZ"/>

</solids>
EOF


# Detector enclosure structure
    print ENCL <<EOF;
<structure>
    <volume name="volFoamPadding">
      <materialref ref="fibrous_glass"/>
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
      x="@{[$DetEncX+2*$RockThickness]}" 
      y="@{[$DetEncY+2*$RockThickness]}" 
      z="@{[$DetEncZ+2*$RockThickness]}"/>
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


print "Some of the principal parameters for this TPC geometry (unit cm unless noted otherwise)\n";
print " CRM active area       : $widthCRM_active x $lengthCRM_active\n";
print " CRM total area        : $widthCRM x $lengthCRM\n";
print " Wire pitch in U, Y, Z : $wirePitchU, $wirePitchY, $wirePitchZ\n";
print " TPC active volume  : $driftTPCActive x $widthTPCActive x $lengthTPCActive\n";
print " Argon volume       : ($Argon_x, $Argon_y, $Argon_z) \n"; 
print " Argon buffer       : ($xLArBuffer, $yLArBuffer, $zLArBuffer) \n"; 
print " Detector enclosure : $DetEncX x $DetEncY x $DetEncZ\n";
print " TPC Origin         : ($OriginXSet, $OriginYSet, $OriginZSet) \n";
print " Field Cage         : $FieldCage_switch \n";
print " Cathode            : $Cathode_switch \n";
print " Workspace          : $workspace \n";
print " Wires              : $wires_on \n";

# run the sub routines that generate the fragments
if ( $FieldCage_switch eq "on" ) {  gen_FieldCage();	}
#if ( $Cathode_switch eq "on" ) {  gen_Cathode();	} #Cathode for now has the same geometry as the Ground Grid

gen_Extend();    # generates the GDML color extension for the refactored geometry 
gen_Define(); 	 # generates definitions at beginning of GDML
gen_Materials(); # generates materials to be used
gen_TPC();       # generate TPC for a given unit CRM
gen_Cryostat();  # 
gen_Enclosure(); # 
gen_World();	 # places the enclosure among DUSEL Rock
write_fragments(); # writes the XML input for make_gdml.pl
		   # which zips together the final GDML
print "--- done\n\n\n";
exit;
