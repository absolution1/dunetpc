#!/usr/bin/perl

use Math::Trig;
use XML::LibXML;
use Getopt::Long;

#Specific H,W,L of LBNE
$CryostatWidth=1900;
$CryostatHeight=1700;
$CryostatLength=8900;

#detector variables
$SteelThickness=0.5*2.54;
$ArgonWidth=$CryostatWidth-2*$SteelThickness;
$ArgonHeight=$CryostatHeight-2*$SteelThickness;
$ArgonLength=$CryostatLength-2*$SteelThickness;
$RockThickness=2000;
$ConcretePadding=50;
$GlassFoamPadding=100;
$TotalPadding=$ConcretePadding+$GlassFoamPadding+$SteelThickness;
$CavernWidth=$ArgonWidth+2*$TotalPadding;
$CavernHeight=$ArgonHeight+$TotalPadding;
$CavernLength=$ArgonLength+2*$TotalPadding;
$TPCWidth=$CryostatWidth-100;
$TPCHeight=0.75*$CryostatHeight;
$TPCLength=0.75*$CryostatLength;
$TPCWireThickness=0.015;
$Pi=3.14159;

gen_header();
gen_materials();
gen_lbne();
exit;






sub gen_lbne()
{
  # Set up the output file.
  $LBNE = "lbne.gdml";
  $LBNE = ">>" . $LBNE;
  open(LBNE) or die("Could not open file $LBNE for writing");

  print LBNE <<EOF;
 <solids>
    <box name="World" lunit="cm" 
      x="$CryostatWidth+2*($RockThickness+$TotalPadding)+$TPCWidth" 
      y="$CryostatHeight+2*($RockThickness+$TotalPadding)" 
      z="$CryostatLength+2*($RockThickness+$TotalPadding)+$TPCLength"/>
    <box name="LowerRock" lunit="cm" 
      x="$CryostatWidth+2*($RockThickness+$TotalPadding)" 
      y="$CryostatHeight+$RockThickness+$TotalPadding" 
      z="$CryostatLength+2*($RockThickness+$TotalPadding)"/>
    <box name="RockBottom" lunit="cm" 
      x="$CryostatWidth+2*$TotalPadding" 
      y="$RockThickness" 
      z="$CryostatLength+2*$TotalPadding"/>
    <box name="LowerCavern" lunit="cm" 
      x="$CryostatWidth+2*$TotalPadding"
      y="$CryostatHeight+$RockThickness+$TotalPadding+2"
      z="$CryostatLength+2*$TotalPadding"/>
    <subtraction name="LowerRockWithCavern">
      <first ref="LowerRock"/>
      <second ref="LowerCavern"/>
    </subtraction>
    <box name="UpperRock" lunit="cm" 
      x="$CryostatWidth+2*($RockThickness+$TotalPadding)" 
      y="$RockThickness" 
      z="$CryostatLength+2*($RockThickness+$TotalPadding)"/>
    <box name="UpperCavern" lunit="cm" 
      x="$CryostatWidth+2*$TotalPadding+1000"
      y="$RockThickness+20"
      z="$CryostatLength+2*$TotalPadding+1000"/>
    <box name="RockTop" lunit="cm" 
      x="$CryostatWidth+2*$TotalPadding+1000" 
      y="$RockThickness-500" 
      z="$CryostatLength+2*$TotalPadding+1000"/>
    <subtraction name="UpperRockWithCavern">
      <first ref="UpperRock"/>
      <second ref="UpperCavern"/>
    </subtraction>
    <box name="DetEnclosure" lunit="cm" 
      x="$CryostatWidth+2*$TotalPadding"
      y="$CryostatHeight+2*$TotalPadding"
      z="$CryostatLength+2*$TotalPadding"/>
    <box name="Concrete" lunit="cm" 
      x="$CryostatWidth+2*$TotalPadding"
      y="$CryostatHeight+$TotalPadding"
      z="$CryostatLength+2*$TotalPadding"/>
    <box name="ConcreteBottom" lunit="cm" 
      x="$CryostatWidth+2*($TotalPadding-$ConcretePadding)"
      y="$ConcretePadding"
      z="$CryostatLength+2*($TotalPadding-$ConcretePadding)"/>
    <box name="ConcreteCavern" lunit="cm" 
      x="$CryostatWidth+2*($TotalPadding-$ConcretePadding)"
      y="$CryostatHeight+$TotalPadding+2"
      z="$CryostatLength+2*($TotalPadding-$ConcretePadding)"/>
    <subtraction name="ConcreteWithCavern">
      <first ref="Concrete"/>
      <second ref="ConcreteCavern"/>
    </subtraction>
    <box name="Cryostat" lunit="cm" 
      x="$CryostatWidth" 
      y="$CryostatHeight" 
      z="$CryostatLength"/>
    <box name="ArgonInterior" lunit="cm" 
      x="$ArgonWidth"
      y="$ArgonHeight"
      z="$ArgonLength"/>
    <subtraction name="SteelShell">
      <first ref="Cryostat"/>
      <second ref="ArgonInterior"/>
    </subtraction>
    <box name="TPC" lunit="cm" 
      x="$TPCWidth" 
      y="$TPCHeight" 
      z="$TPCLength"/>
    <box name="TPCPlane" lunit="cm" 
      x="0.1" 
      y="0.9*$TPCHeight" 
      z="0.9*$TPCLength"/>
    <tube name="TPCWire"
      rmax="0.5*$TPCWireThickness"
      z="0.89*$TPCHeight"               
      deltaphi="2*$Pi"
      aunit="rad"
      lunit="cm"/>
  </solids>
  <structure>
    <volume name="volSteelShell">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="SteelShell" />
    </volume>
    <volume name="volTPCWire">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
      <solidref ref="TPCWire" />
    </volume>
    <volume name="volTPCPlane">
      <materialref ref="LAr"/>
      <solidref ref="TPCPlane"/>
EOF

    $count=1;
    for ( $i=-0.3*$TPCLength ; $i < 0.3*$TPCLength ; $i+=500 ) {
      $wire_zpos=$i;
      print LBNE <<EOF;
      <physvol>
        <volumeref ref="volTPCWire"/>
        <position name="posTPCWire$count" unit="cm" x="0" y="0" z="$wire_zpos"/>
        <rotation name="rTPCWire$count" unit="deg" x="60" y="0" z="0"/>
      </physvol>
EOF
      $count++;
    }
    print LBNE <<EOF;
    </volume>
    <volume name="volTPC">
      <materialref ref="LAr" />
      <solidref ref="TPC" />
      <physvol>
        <volumeref ref="volTPCPlane"/>
        <position name="posTPCPlane1" unit="cm" x="-0.45*$TPCWidth" y="0" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volTPCPlane"/>
        <position name="posTPCPlane2" unit="cm" x="-0.475*$TPCWidth" y="0" z="0"/>
        <rotation name="rTPCPlane2" unit="deg" x="0" y="180" z="0"/>
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
    <volume name="volConcreteWithCavern">
      <materialref ref="Concrete"/>
      <solidref ref="ConcreteWithCavern"/>
    </volume>
    <volume name="volConcreteBottom">
      <materialref ref="Concrete"/>
      <solidref ref="ConcreteBottom"/>
    </volume>
    <volume name="volDetEnclosure">
      <materialref ref="Air"/>
      <solidref ref="DetEnclosure"/>
      <physvol>
        <volumeref ref="volCryostat"/>
        <position name="posCryostat" unit="cm" x="0" y="0" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volConcreteWithCavern"/>
        <position name="posConcreteWithCavern" unit="cm" x="0" y="-0.5*($TotalPadding-$ConcretePadding)" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volConcreteBottom"/>
        <position name="posConcreteBottom" unit="cm" x="0" y="-0.5*($CryostatHeight+$TotalPadding)" z="0"/>
      </physvol>
    </volume>
    <volume name="volLowerRockWithCavern">
      <materialref ref="DUSEL_Rock"/>
      <solidref ref="LowerRockWithCavern"/>
    </volume>
    <volume name="volRockTop">
      <materialref ref="DUSEL_Rock"/>
      <solidref ref="RockTop"/>
    </volume>
    <volume name="volUpperRockWithCavern">
      <materialref ref="DUSEL_Rock"/>
      <solidref ref="UpperRockWithCavern"/>
    </volume>
    <volume name="volRockBottom">
      <materialref ref="DUSEL_Rock"/>
      <solidref ref="RockBottom"/>
    </volume>
    <volume name="volWorld" >
      <materialref ref="Air"/>
      <solidref ref="World"/>
      <physvol>
        <volumeref ref="volDetEnclosure"/>
        <position name="posDetEnclosure" unit="cm" x="0.5*$TPCWidth" y="0" z="0.5*$TPCLength"/>
      </physvol>
      <physvol>
        <volumeref ref="volLowerRockWithCavern"/>
        <position name="posLowerRockWithCavern" unit="cm" x="0.5*$TPCWidth" y="-0.5*$RockThickness" z="0.5*$TPCLength"/>
      </physvol>
      <physvol>
        <volumeref ref="volRockBottom"/>
        <position name="posRockBottom" unit="cm" x="0.5*$TPCWidth" y="-0.5*($RockThickness+$CryostatHeight+$TotalPadding)" z="0.5*$TPCLength"/>
      </physvol>
      <physvol>
        <volumeref ref="volUpperRockWithCavern"/>
        <position name="posLowerRockWithCavern" unit="cm" x="0.5*$TPCWidth" y="0.5*($RockThickness+$CryostatHeight)" z="0.5*$TPCLength"/>
      </physvol>
      <physvol>
        <volumeref ref="volRockTop"/>
        <position name="posRockTop" unit="cm" x="0.5*$TPCWidth" y="0.5*($RockThickness-500+$CryostatHeight)+500" z="0.5*$TPCLength"/>
      </physvol>
    </volume>
  </structure>
  <setup name="Default" version="1.0">
    <world ref="volWorld" />
  </setup>
</gdml>
EOF
   close(LBNE);
}


#generates necessary gd/xml header
sub gen_header() 
{

  $LBNE = "lbne.gdml";
  $LBNE = ">" . $LBNE;
  open(LBNE) or die("Could not open file $LBNE for writing");
  print LBNE <<EOF;
<?xml version="1.0" encoding="UTF-8" ?>
<gdml xmlns:gdml="http://cern.ch/2001/Schemas/GDML"
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:noNamespaceSchemaLocation="GDMLSchema/gdml.xsd">
EOF
}


sub gen_materials() 
{


  $LBNE = "lbne.gdml";
  $LBNE = ">>" . $LBNE;
  open(LBNE) or die("Could not open file $LBNE for writing");
  print LBNE <<EOF;
  <materials>
    <element name="videRef" formula="VACUUM" Z="1">  <atom value="1"/> </element>
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
    <element name="sodium" formula="Na" Z="11"> <atom value="22.99"/>    </element>
    <element name="argon" formula="Ar" Z="18"> <atom value="39.9480"/>  </element>

    <material name="Vacuum" formula="Vacuum">
      <D value="1.e-25" unit="g/cm3"/>
      <fraction n="1.0" ref="videRef"/>
    </material>
  
    <material name="STEEL_STAINLESS_Fe7Cr2Ni" formula="STEEL_STAINLESS_Fe7Cr2Ni">
      <D value="7.9300" unit="g/cm3"/>
      <fraction n="0.0010" ref="carbon"/>
      <fraction n="0.1800" ref="chromium"/>
      <fraction n="0.7298" ref="iron"/>
      <fraction n="0.0900" ref="nickel"/>
    </material>
  
    <material name="LAr" formula="LAr">
      <D value="1.40" unit="g/cm3"/>
      <fraction n="1.0000" ref="argon"/>
    </material>
  
    <material formula=" " name="Air">
      <D value="0.001205" unit="g/cc"/>
      <fraction n="0.78084" ref="nitrogen"/>
      <fraction n="0.209476" ref="oxygen"/>
      <fraction n="0.00934" ref="argon"/>
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

    <material formula=" " name="DUSEL_Rock">
      <D value="2.82" unit="g/cc"/>
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

    <material formula=" " name="ShotRock">
      <D value="2.7*0.6" unit="g/cc"/>
      <fraction n="0.438" ref="oxygen"/>
      <fraction n="0.257" ref="silicon"/>
      <fraction n="0.222" ref="sodium"/>
      <fraction n="0.049" ref="aluminum"/>
      <fraction n="0.019" ref="iron"/>
      <fraction n="0.015" ref="potassium"/>
    </material>


  </materials>
EOF
}
