#!/usr/bin/perl

use Math::Trig;
use XML::LibXML;
use Getopt::Long;

$outfile					=	"dunend.gdml";
#Specific H,W,L of DUNE
$CryostatLength		=	538.5;
$CryostatRadius		=	312.4/2;
$CryostatWidth 	        =       2*$CryostatRadius;

#detector variables
$ArgonLength		=	510.6;
$ArgonRadius			=  279.4/2;


$TPCWidth					= 179.4;
$TPCHeight				=	179.4;
$TPCLength				=	399.9;

$CathodeThickness = 0.1;
$PlanesSpacing = 0.3;

$cos30=0.86602540378;
$sin30=0.5;
$tg30=0.5773502692;


$TPCWireThickness=0.015;
$TPC_WireSpacing=0.3;

#$TPCActiveWidth					= $TPCWidth - $CathodeThickness - 2*$PlanesSpacing - 2*$TPCWireThickness;
#$TPCActiveHeight				=	$TPCHeight -0.5;
#$TPCActiveLength				=	$TPCLength - 0.5;
$TPCActiveWidth					= $TPCWidth;
$TPCActiveHeight				=	$TPCHeigth;
$TPCActiveLength				=	$TPCLength;


$TPC_NLongWires=148;
$TPC_NTopWires=518;
$LongWiresSpread=89.1;
$LongWiresPadding=($LongWiresSpread - ($TPC_WireSpacing/$sin30*($TPC_NLongWires - 1)))/2;
$TBWiresPadding=$TPC_WireSpacing/$sin30 - $LongWiresPadding;
$TPC_NBottomWires=518;
$TPC_TotNWiresZ=$TPC_NTopWires+$TPC_NLongWires;
$WiresPadding=($TPCLength - ($TPC_WireSpacing/$sin30*($TPC_TotNWiresZ - 1)))/2;
$TPC_LongWireLength=2*$TPCHeight -0.1;	
$TPC_VertWireLength=$TPCHeight;	
$TPC_NVertWires = 1333;
$VertWiresPadding=($TPCLength - ($TPC_WireSpacing*($TPC_NVertWires - 1)))/2;

$Pi=3.14159;

gen_header();
gen_materials();
gen_dune();
exit;






sub gen_dune()
{
	# Set up the output file.
	$DUNE = ">>" . $outfile;
	open(DUNE) or die("Could not open file $DUNE for writing");

	print DUNE <<EOF;
 <solids>
    <box name="World" lunit="cm" 
      x="2*$CryostatWidth+100" 
      y="4*$CryostatRadius+100" 
      z="2*$CryostatLength + 100"/>
     <box name="DetEnclosure" lunit="cm" 
      x="$CryostatWidth+2"
      y="2*$CryostatRadius+2"
      z="$CryostatLength+2"/>
    <box name="TPCActive" lunit="cm" 
      x="$TPCActiveWidth"
      y="$TPCActiveHeigth"
      z="$TPCActiveLength"/>
      <tube name="Cryostat" lunit="cm" aunit="deg" 
      rmax="$CryostatRadius" 
      z="$CryostatLength" 
      deltaphi="360"/>
<tube name="CryostatInt" lunit="cm" aunit="deg" 
      rmax="$CryostatRadius - 0.01" 
      z="$CryostatLength - 0.01" 
      deltaphi="360"/>

      <tube name="ArgonInterior" lunit="cm" aunit="deg" 
      rmax="$ArgonRadius"
      z="$ArgonLength"
      deltaphi="360"/>
      <subtraction name="SteelShell">
      <first ref="CryostatInt"/>
      <second ref="ArgonInterior"/>
      </subtraction>
      <box name="TPC" lunit="cm" 
      x="$TPCWidth" 
      y="$TPCHeight" 
      z="$TPCLength"/>
      <box name="TPCPlane" lunit="cm" 
      x="0.1" 
      y="$TPCHeight" 
      z="$TPCLength"/>
      <tube name="TPCLongWire"
      rmax="0.5*$TPCWireThickness"
      z="$TPC_LongWireLength"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>

      <tube name="TPCVertWire"
      rmax="0.5*$TPCWireThickness"
      z="$TPC_VertWireLength"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>


EOF
	#lengths of the top wires: the 0.03 is needed due to a small approximation on cos30
	$countertop=1;
	@topzpos;
	for ( $i= $WiresPadding; $i <= $TPCLength- $LongWiresSpread; $i+=$TPC_WireSpacing/$sin30 ) {
		$wire_zpos=$i;
		print DUNE <<EOF;
	<tube name="TPCWireTop$countertop"
	rmax="0.5*$TPCWireThickness"
	z="$wire_zpos/$cos30 -0.03"
	deltaphi="360"
	aunit="deg"
	lunit="cm"/>

EOF
		push(@topzpos,$wire_zpos);
		$countertop++;
	}
	#lengths of the bottom wires. The 0.015 is needed due to very timy approximation effects in the calculation of $i
	#the 0.03 is needed due to a small approximation on cos30

	$counterbottom=1;
	@bottomzpos;
	for ( $i= $LongWiresSpread + $TBWiresPadding; $i <= ($TPCLength - $WiresPadding + 0.015); $i+=$TPC_WireSpacing/$sin30 ) {
		$wire_zpos=$i;
		print DUNE <<EOF;
	<tube name="TPCWireBottom$counterbottom"
	rmax="0.5*$TPCWireThickness"
	z="($TPCLength - $wire_zpos)/$cos30 -0.03"
	deltaphi="360"
	aunit="deg"
	lunit="cm"/>

EOF

		push(@bottomzpos,$TPCLength - $wire_zpos);
		$counterbottom++;
	}


	print DUNE <<EOF;
      <box name="Cathode" lunit="cm"
      x="$CathodeThickness"
      y="$TPCHeight"
      z="$TPCLength"/>
      </solids>
      <structure>
      <volume name="volTPCLongWire">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCLongWire" />
      </volume>
      <volume name="volTPCVertWire">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCVertWire" />
      </volume>
      <volume name="volSteelShell">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="SteelShell" />
      </volume>


EOF
	#top wires

	for($i = 1; $i < $countertop; $i+=1){
		print DUNE<<EOF;
<volume name="volTPCWireTop$i">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCWireTop$i" />
      </volume>
EOF
	}
	#bottom wires

	for($i = 1; $i <$counterbottom; $i+=1){

		print DUNE<<EOF;
      <volume name="volTPCWireBottom$i">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCWireBottom$i" />
      </volume>
EOF

	}

	print DUNE<<EOF;

    <volume name="volTPCPlane">
    <materialref ref="LAr" />
    <solidref ref="TPCPlane" />
EOF

	#place long wires
	#The 0.015 is needed due to very timy approximation effects in the calculation of $i
	$count=1;
	for ( $i= (-$LongWiresSpread + 2*$LongWiresPadding)/2 ; $i <=  ($LongWiresSpread - 2*$LongWiresPadding)/2 + 0.015 ; $i+=$TPC_WireSpacing/$sin30 ) {
		$wire_zpos=$i;
		print DUNE <<EOF;
	 <physvol>
	 <volumeref ref="volTPCLongWire"/>
	 <position name="posTPCLWire$count" unit="cm" x="0" y="0" z="$wire_zpos "/>
	 <rotation name="rTPCLWire$count" unit="deg" x="30" y="0" z="0"/>
	 </physvol>
EOF
		$count++;
	}
	$longwiresplaced = $count - 1;
	if($longwiresplaced != $TPC_NLongWires){
		print "number of long wires placed = $longwiresplaced  differs from the nominal one = $TPC_NLongWires\n";
	}
	#place top wires
	$count=1;
	for ( $i= (-$TPCLength + 2*$WiresPadding)/2 ; $i <= $TPCLength/2 - $LongWiresSpread; $i+=$TPC_WireSpacing/$sin30 ) {
		$wire_zpos=$i - ($topzpos[$count-1])/2;
		$wire_ypos= $TPCHeight/2 - ($topzpos[$count-1])/2*$tg30;

		print DUNE <<EOF;
	 <physvol>
	 <volumeref ref="volTPCWireTop$count"/>
	 <position name="posTPCTWire$count" unit="cm" x="0" y="$wire_ypos" z="$wire_zpos"/>
	 <rotation name="rTPCTWire$count" unit="deg" x="30" y="0" z="0"/>
	 </physvol>
EOF
		$count++;
	}
	$topwiresplaced = $count -1;
	if($topwiresplaced != $TPC_NTopWires){
		print "number of top wires placed = $topwiresplaced  differs from the nominal one = $TPC_NTopWires\n";
	}

	#place bottom wires
	#The 0.015 is needed due to very timy approximation effects in the calculation of $i
	$count=1;
	for ( $i= -$TPCLength/2 + $LongWiresSpread + $TBWiresPadding ; $i <=  ($TPCLength - 2*$WiresPadding)/2 ; $i+=$TPC_WireSpacing/$sin30 ) {
		$wire_zpos=$i + ($bottomzpos[$count-1])/2;
		$wire_ypos= -$TPCHeight/2 + ($bottomzpos[$count-1])/2*$tg30;
		print DUNE <<EOF;
	 <physvol>
	 <volumeref ref="volTPCWireBottom$count"/>
	 <position name="posTPCBWire$count" unit="cm" x="0" y="$wire_ypos" z="$wire_zpos"/>
	 <rotation name="rTPCBWire$count" unit="deg" x="30" y="0" z="0"/>
	 </physvol>
EOF
		$count++;
	}

	$bottomwiresplaced = $count - 1;
	if($bottomwiresplaced != $TPC_NBottomWires){
		print "number of bottom wires placed = $bottomwiresplaced  differs from the nominal one = $TPC_NBottomWires\n";
	}


	print DUNE <<EOF;

</volume>
<volume name="volTPCPlaneVert">
    <materialref ref="LAr" />
    <solidref ref="TPCPlane" />

EOF
	#place vert wires
	#The 0.015 is needed due to very timy approximation effects in the calculation of $i
	$count=1;
	for ( $i= (-$TPCLength + 2*$VertWiresPadding)/2 ; $i <= ($TPCLength - 2*$VertWiresPadding)/2 + 0.015 ; $i+=$TPC_WireSpacing ) {
		$wire_zpos=$i;

		print DUNE <<EOF;
	 <physvol>
		 <volumeref ref="volTPCVertWire"/>
		 <position name="posTPCVWire$count" unit="cm" x="0" y="0" z="$wire_zpos "/>
		 <rotation name="rTPCVWire$count" unit="deg" x="90" y="0" z="0"/>
		 </physvol>
EOF
		$count++;
	}
	$vertwiresplaced = $count - 1;
	if($vertwiresplaced != $TPC_NVertWires){
		print "number of vertical wires placed = $vertwiresplaced  differs from the nominal one = $TPC_NVertWires\n";
	}


	print DUNE <<EOF;
 </volume>
	 <volume name="volCathode">
	 <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
	 <solidref ref="Cathode" />
	 </volume>
	 <volume name="volTPCActive">
	 <materialref ref="LAr" />
	 <solidref ref="TPC" />
	 <physvol>
	 <volumeref ref="volTPCPlane"/>
	 <position name="posTPCPlane1" unit="cm" x="(-$TPCWidth/2)+0.051+2*$PlanesSpacing" y="0" z="0"/>
	 </physvol>
	 <physvol>
	 <volumeref ref="volTPCPlane"/>
	 <position name="posTPCPlane2" unit="cm" x="(-$TPCWidth/2)+0.051 + $PlanesSpacing" y="0" z="0"/>
	 <rotation name="rTPCPlane2" unit="deg" x="0" y="180" z="0"/>
	 </physvol>
	 <physvol>
	<volumeref ref="volTPCPlaneVert"/>
	 <position name="posTPCPlaneVert" unit="cm" x="(-$TPCWidth/2)+0.051" y="0" z="0"/>
	 </physvol>
	 <physvol>
	 <volumeref ref="volCathode"/>
	 <position name="posCathode" unit="cm" x="($TPCWidth/2)-$CathodeThickness/2" y="0" z="0"/>
	 </physvol>
	 </volume>
         <volume name="volTPC">
	<materialref ref="LAr" />
	<solidref ref="TPC" />
	<physvol>
	<volumeref ref="volTPCActive"/>
	<position name="posTPCActive" unit="cm" x="0" y="0" z="0"/>
	</physvol>
	 </volume>
	<volume name="volCryostat">
	<materialref ref="LAr" />
	<solidref ref="Cryostat" />
	<physvol>
	<volumeref ref="volSteelShell"/>
	<position name="posSteelShell" unit="cm" x="0" y="0" z="0"/>
	</physvol>
	 <physvol>
	 <volumeref ref="volTPC"/>
	<position name="posTPC" unit="cm" x="0" y="0" z="0"/>
	</physvol>
	 </volume>
	 <volume name="volDetEnclosure">
	<materialref ref="Air"/>
	<solidref ref="DetEnclosure"/>
	<physvol>
	<volumeref ref="volCryostat"/>
	 <position name="posCryostat" unit="cm" x="0" y="0" z="0"/>
	 </physvol>
	 </volume>
	 <volume name="volWorld" >
	 <materialref ref="Air"/>
	 <solidref ref="World"/>
	 <physvol>
	 <volumeref ref="volDetEnclosure"/>
	 <position name="posDetEnclosure" unit="cm" x="0" y="0" z="0"/>
	 </physvol>
	 </volume>
	 </structure>
	 <setup name="Default" version="1.0">
	 <world ref="volWorld" />
	 </setup>
	 </gdml>
EOF
	close(DUNE);
}


#generates necessary gd/xml header
sub gen_header() 
{

	$DUNE = ">" . $outfile;
	open(DUNE) or die("Could not open file $DUNE for writing");
	print DUNE <<EOF;
	<?xml version="1.0" encoding="UTF-8" ?>
		<gdml xmlns:gdml="http://cern.ch/2001/Schemas/GDML"
		xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
		xsi:noNamespaceSchemaLocation="GDMLSchema/gdml.xsd">
EOF
}


sub gen_materials() 
{

	$DUNE = ">>" . $outfile;
	open(DUNE) or die("Could not open file $DUNE for writing");
	print DUNE <<EOF;
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
		<D value="2.82" unit="g/cc"/>
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

		<material formula=" " name="Air">
		<D value="0.001205" unit="g/cc"/>
		<fraction n="0.781154" ref="nitrogen"/>
		<fraction n="0.209476" ref="oxygen"/>
		<fraction n="0.00934" ref="argon"/>
		</material>

		<material formula=" " name="G10">
		<D value="1.7" unit="g/cc"/>
		<fraction n="0.2805" ref="silicon"/>
		<fraction n="0.3954" ref="oxygen"/>
		<fraction n="0.2990" ref="carbon"/>
		<fraction n="0.0251" ref="hydrogen"/>
		</material>

		<material formula=" " name="Granite">
		<D value="2.7" unit="g/cc"/>
		<fraction n="0.438" ref="oxygen"/>
		<fraction n="0.257" ref="silicon"/>
		<fraction n="0.222" ref="sodium"/>
		<fraction n="0.049" ref="aluminum"/>
		<fraction n="0.019" ref="iron"/>
		<fraction n="0.015" ref="potassium"/>
		</material>

		<material formula=" " name="ShotRock">
		<D value="2.7*0.6" unit="g/cc"/>
		<fraction n="0.438" ref="oxygen"/>
		<fraction n="0.257" ref="silicon"/>
		<fraction n="0.222" ref="sodium"/>
		<fraction n="0.049" ref="aluminum"/>
		<fraction n="0.019" ref="iron"/>
		<fraction n="0.015" ref="potassium"/>
		</material>

		<material formula=" " name="Dirt">
		<D value="1.7" unit="g/cc"/>
		<fraction n="0.438" ref="oxygen"/>
		<fraction n="0.257" ref="silicon"/>
		<fraction n="0.222" ref="sodium"/>
		<fraction n="0.049" ref="aluminum"/>
		<fraction n="0.019" ref="iron"/>
		<fraction n="0.015" ref="potassium"/>
		</material>

		<material formula=" " name="Concrete">
		<D value="2.3" unit="g/cc"/>
		<fraction n="0.530" ref="oxygen"/>
		<fraction n="0.335" ref="silicon"/>
		<fraction n="0.060" ref="calcium"/>
		<fraction n="0.015" ref="sodium"/>
		<fraction n="0.020" ref="iron"/>
		<fraction n="0.040" ref="aluminum"/>
		</material>

		<material formula="H2O" name="Water">
		<D value="1.0" unit="g/cc"/>
		<fraction n="0.1119" ref="hydrogen"/>
		<fraction n="0.8881" ref="oxygen"/>
		</material>

		<material formula="Ti" name="Titanium">
		<D value="4.506" unit="g/cc"/>
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
}
