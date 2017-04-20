#!/bin/sh

# On Alex Himmel's request -- turn all of the photon detector Acrylic to LAr
# in the DUNE 10kt *_nowires.gdml files, in order to make photon library generation
# easier.

# some assumptions -- the photon detectors are the only acrylic in these files.
# this is the case for the 10kt geometries as of April 2017, but we'll skip doing
# this for ProtoDUNE's and 35t -- they have acrylic cosmic-ray taggers, and their
# photon detector volumes are already LAr in the *_nowires.gdml files for them.

# to do later -- de-acrylify only those Acrylic lines right after the optical
# detector names

tmpfile=deacrylifytmpgdml.txt

for filename in dune10kt*_nowires.gdml
do
  echo Processing $filename
  touch $tmpfile
  rm $tmpfile
  sed -e 's/ref="Acrylic"/ref="LAr"/' $filename > $tmpfile
  mv $tmpfile $filename
done

