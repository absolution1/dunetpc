#include "UnpackFragment.h"

// Format correct as of 17 Mar 15 (M Wallbank)

//Which FragmentID would contain this channel
unsigned int UnpackFragment::getFragIDForChan(unsigned int channel){

  /*
    channels    fragment ID
    [0:127]     100
    [128:255]   101
    [256:383]   102
    [384:511]   103
    [512:639]   104
    [640:767]   105
    [768:895]   106
    [896:1023]  107
    [1024:1151] 108
    [1152:1279] 109
    [1280:1407] 110
    [1408:1535] 111
    [1536:1663] 112
    [1664:1791] 113
    [1792:1919] 114
    [1920:2047] 115
  */

  return (channel / 128) + 100;

}

//Which Nano-slice sample number would contain this channel
unsigned int UnpackFragment::getNanoSliceSampleForChan(unsigned int channel){

  /*   
    channels    sample
    [0:127]     0-127
    [128:255]   0-127
    [256:383]   0-127
    [384:511]   0-127
    [512:639]   0-127
    [640:767]   0-127
    [768:895]   0-127
    [896:1023]  0-127
    [1024:1151] 0-127
    [1152:1279] 0-127
    [1280:1407] 0-127
    [1408:1535] 0-127
    [1536:1663] 0-127
    [1664:1791] 0-127
    [1792:1919] 0-127
    [1920:2047] 0-127
  */
  
  unsigned int chanNumInFrag = channel % 128;
  return chanNumInFrag;

}

