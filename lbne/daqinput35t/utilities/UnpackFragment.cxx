#include "UnpackFragment.h"


//Which FragmentID would contain this channel
unsigned int UnpackFragment::getFragIDForChan(unsigned int channel){

  /*
    For now assume 
    
    channels   fragment ID
    [0:127]    0
    [128:255]  1

  */

  return channel / 128;

}

//Which Nano-slice group would contain this channel
unsigned int UnpackFragment::getNanoSliceGroupForChan(unsigned int channel){

  /*
    For now assume
    
    channels  group
    [0:31]    0
    [32:63]   1
    [64:95]   2
    [96:127]  3
    
  */
  
  unsigned int chanNumInFrag = channel % 128;
  return chanNumInFrag / 32;
  
}

//Which Nano-slice sample number would contain this channel
unsigned int UnpackFragment::getNanoSliceSampleForChan(unsigned int channel){

  /*
    For now assume
    
    channels  group sample
    [0:31]    0     [0:31]
    [32:63]   1     [0:31]     
    [64:95]   2     [0:31]
    [96:127]  3     [0:31]
    
  */
  
  unsigned int chanNumInFrag = channel % 128;
  return chanNumInFrag % 32;

}

