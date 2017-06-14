#ifndef UnpackFragment_h
#define UnpackFragment_h

#include "artdaq-core/Data/Fragment.hh"

/*

  Utility functions for getting information out of the online data format

  Currently assumed:
  16 RCEs (Reconfigurable Cluster Elements) running board readers
  1  ART-DAQ Fragment per RCE
  Each RCE sends fragments that are milli-slices for 128 channels (i.e. 1 FEB [Front End Board])
  A milli-slice is some integer multiple, N,  of drift times long
  A milli-slice contains N micro-slices, each 1 drift time long
  A micro-slice contains M nano-slices, where M is the number of TDCs in a drift window
  A nano-slice is 1 time tick for 128 channels (i.e. has 128 ADC values, corresponding to ADC at that TDC)

  Future plans:
  -Probably want to read in the geometry upon intialisation from a text file or database  

*/



namespace UnpackFragment{

  //Each channel's ADC values can be uniquely addressed by: [FragmentID, nSlice Sample]

  //Which FragmentID would contain this channel
  unsigned int getFragIDForChan(unsigned int channel);

  //Which Nano-slice sample number would contain this channel
  unsigned int getNanoSliceSampleForChan(unsigned int channel);
  
}



#endif //UnpackFragment_h
