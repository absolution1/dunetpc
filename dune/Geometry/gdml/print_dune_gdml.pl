#!/usr/bin/perl

    $RockMaterial="ShotRock";	# current options: "Granite", "Dirt", "ShotRock" 

    $DetectorWidth=1900;
    $DetectorHeight=1700;
    $DetectorLength=8900;
    $CavernWidth=$DetectorWidth+2*$TotalPadding;
    $CavernHeight=$DetectorHeight+$TotalPadding;
    $CavernLength=$DetectorLength+2*$TotalPadding;
    $RockThickness=2000;
    $ConcretePadding=50;
    $GlassFoamPadding=100;
    $TotalPadding=$ConcretePadding+$GlassFoamPadding;


$GDML = "dune.gdml";
$GDML = ">" . $GDML; 
open(GDML) or die("Could not open file $GDML for writing");

print GDML <<EOF;
<?xml version="1.0" encoding="UTF-8" ?>
<gdml xmlns:gdml="http://cern.ch/2001/Schemas/GDML"
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:noNamespaceSchemaLocation="GDMLSchema/gdml.xsd">
<define>
  <constant name="kInch"	value="2.54" />
  <constant name="kPi"	value="3.14159265358979" />
  <constant name="kSteelThickness"	  value="0.5*kInch" />
</define>

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
   <fraction n="0.1800" ref="chromium"/>
   <fraction n="0.7298" ref="iron"/>
   <fraction n="0.0900" ref="nickel"/>
  </material>

  <material name="LAr" formula="LAr">
   <D value="1.40" unit="g/cm3"/>
   <fraction n="1.0000" ref="argon"/>
  </material>

  <material formula=" " name="Air">
   <D value="0.001205" unit="g/cm3"/>
   <fraction n="0.78084" ref="nitrogen"/>
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

<solids>
   <box name="World" lunit="cm" 
     x="$CavernWidth+$RockThickness" 
     y="$CavernHeight+$RockThickness" 
     z="$CavernLength+$RockThickness"/>
   <box name="Rock" lunit="cm" 
     x="$CavernWidth+$RockThickness" 
     y="$CavernHeight+$RockThickness" 
     z="$CavernLength+$RockThickness"/>
   <box name="RockBottom" lunit="cm" 
     x="$CavernWidth" 
     y="$RockThickness" 
     z="$CavernLength"/>
   <box name="Cavern" lunit="cm" 
     x="$CavernWidth"
     y="$CavernHeight+2*$RockThickness+100"
     z="$CavernLength"/>
   <subtraction name="RockWithCavern">
     <first ref="Rock"/>
     <second ref="Cavern"/>
   </subtraction>
   <box name="Concrete" lunit="cm" 
     x="1900" 
     y="1700" 
     z="8900"/>
   <box name="ConcreteBottom" lunit="cm" 
     x="1850" 
     y="50" 
     z="8850"/>
   <box name="ConcreteCavern" lunit="cm" 
     x="1850"
     y="1800"
     z="8850"/>
   <subtraction name="ConcreteWithCavern">
     <first ref="Concrete"/>
     <second ref="ConcreteCavern"/>
   </subtraction>
   <box name="SteelBox" lunit="cm" 
     x="1900-2*(150)" 
     y="1600" 
     z="8900-2*(150)"/>
   <box name="ArgonInterior" lunit="cm" 
     x="1900-2*150-kSteelThickness"
     y="1600-kSteelThickness"
     z="1700-2*150-kSteelThickness"/>
   <subtraction name="CryostatSteelShell">
     <first ref="SteelBox"/>
     <second ref="ArgonInterior"/>
   </subtraction>
</solids>

<structure>
  <volume name="volArgonInterior">
    <materialref ref="LAr" />
    <solidref ref="ArgonInterior" />
  </volume>
  <volume name="volCryostatSteelShell">
    <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni" />
    <solidref ref="CryostatSteelShell" />
  </volume>
  <volume name="volConcreteWithCavern">
    <materialref ref="Concrete" />
    <solidref ref="ConcreteWithCavern" />
  </volume>
  <volume name="volConcreteBottom">
    <materialref ref="Concrete" />
    <solidref ref="ConcreteBottom" />
  </volume>
  <volume name="volRockWithCavern">
    <materialref ref="$RockMaterial" />
    <solidref ref="RockWithCavern" />
  </volume>
  <volume name="volRockBottom">
    <materialref ref="$RockMaterial" />
    <solidref ref="RockBottom" />
  </volume>
  <volume name="volWorld" >
    <materialref ref="Air"/>
    <solidref ref="World"/>
      <physvol>
        <volumeref ref="volCryostatSteelShell"/>
        <position name="posCryostatSteelShell" unit="cm" x="0" y="250+(1550/2)" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volArgonInterior"/>
        <position name="posArgonInterior" unit="cm" x="0" y="250+(1550/2)" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volConcreteWithCavern"/>
        <position name="posConcreteWithCavern" unit="cm" x="0" y="150+(1700/2)" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volConcreteBottom"/>
        <position name="posConcreteWithCavern" unit="cm" x="0" y="150+(50/2)" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volRockWithCavern"/>
        <position name="posRockWithCavern" unit="cm" x="0" y="0" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volRockBottom"/>
        <position name="posRockBottom" unit="cm" x="0" y="-850" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volRockWithCavern"/>
        <position name="posRockWithCavern" unit="cm" x="0" y="0" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volRockBottom"/>
        <position name="posRockBottom" unit="cm" x="0" y="-850" z="0"/>
      </physvol>
  </volume>
</structure>

<setup name="Default" version="1.0">
  <world ref="volWorld" />
</setup>

</gdml>

EOF
