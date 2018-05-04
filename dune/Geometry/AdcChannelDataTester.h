// AdcChannelDataTester.h
//
// David Adams
// April 2018
//
// Class with functions to fill AdcChannelData for test purposes.
//
// Wire strums put data starting at tick i for channel or wire i.
// Where needed, the geometry is obtained from the geometry service.

#ifndef AdcChannelDataTester_H
#define AdcChannelDataTester_H

#include "dune/DuneInterface/AdcChannelData.h"

class AdcChannelDataTester {

public:

  using Index = unsigned int;

  // Data used to fill each channel data object.
  Index run =1;
  Index subrun =0;
  Index event =0;
  Index channel = 0;
  float pedestal =700.0;
  Index nsam = 500;
  std::vector<short> rawwf ={1200};
  std::vector<float> samwf ={50.0};

  // Fill ADC data with passed channel number and the data associated with this tester.
  // If nneded, the raw and prepared signals are resized to nsam.
  // The waveforms offset by isam0 are added to the raw and prepared signals in ADC data.
  int fill(AdcChannelData& acd, Index icha, Index isam0);

  // Strum channels 0 to ncha.
  // If ncha ==0 , then # channels is taken from the geometry.
  int strumChannels(AdcChannelDataMap& acds, Index ncha =0);

  // Strum all wires in TPC itpc in crystat icry.
  int strumTpcWires(AdcChannelDataMap& acds, Index itpc, Index icry =0);

  // Strum all wires in a cryostat.
  int strumCryoWires(AdcChannelDataMap& acds, Index icry =0);

  // Strum all wires in the detector.
  int strumDetectorWires(AdcChannelDataMap& acds);

};

#endif
